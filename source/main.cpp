// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: Copyright 2019-2022 Heal Research

#include <chrono>
#include <cmath>
#include <cstdlib>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <memory>
#include <thread>
#include <taskflow/taskflow.hpp>
#if TF_MINOR_VERSION > 2
#include <taskflow/algorithm/reduce.hpp>
#endif
#include "operon/algorithms/nsga2.hpp"
#include "operon/core/format.hpp"
#include "operon/core/version.hpp"
#include "operon/core/problem.hpp"
#include "operon/interpreter/interpreter.hpp"
#include "operon/operators/creator.hpp"
#include "operon/operators/crossover.hpp"
#include "operon/operators/evaluator.hpp"
#include "operon/operators/generator.hpp"
#include "operon/operators/initializer.hpp"
#include "operon/operators/mutation.hpp"
#include "operon/operators/non_dominated_sorter.hpp"
#include "operon/operators/reinserter.hpp"
#include "operon/operators/selector.hpp"

#include "util.hpp"
#include "operator_factory.hpp"

// atomic
#include "atomic.hpp"

auto main(int argc, char** argv) -> int // NOLINT
{
    auto opts = Operon::InitOptions("operon_gp", "Genetic programming symbolic regression");
    auto result = Operon::ParseOptions(std::move(opts), argc, argv);

    // parse and set default values
    Operon::GeneticAlgorithmConfig config{};
    config.Generations = result["generations"].as<size_t>();
    config.PopulationSize = result["population-size"].as<size_t>();
    config.PoolSize = result["pool-size"].as<size_t>();
    config.Epsilon = result["epsilon"].as<Operon::Scalar>();
    config.Evaluations = result["evaluations"].as<size_t>();
    config.Iterations = result["iterations"].as<size_t>();
    config.CrossoverProbability = result["crossover-probability"].as<Operon::Scalar>();
    config.MutationProbability = result["mutation-probability"].as<Operon::Scalar>();
    config.TimeLimit = result["timelimit"].as<size_t>();
    config.Seed = std::random_device {}();

    // parse remaining config options
    Operon::Range training_range;
    Operon::Range test_range;
    std::unique_ptr<Operon::Dataset> dataset;
    std::string target;
    bool show_pset = false;
    auto threads = std::thread::hardware_concurrency();
    Operon::NodeType pset_cfg = Operon::PrimitiveSet::Arithmetic;

    auto max_length = result["maxlength"].as<size_t>();
    auto max_depth = result["maxdepth"].as<size_t>();
    auto crossover_internal_probability = result["crossover-internal-probability"].as<Operon::Scalar>();

    auto symbolic = result["symbolic"].as<bool>();

    try {
        for (const auto& kv : result.arguments()) {
            const auto& key = kv.key();
            const auto& value = kv.value();

            if (key == "dataset") {
                dataset = std::make_unique<Operon::Dataset>(value, true);
                ENSURE(!dataset->IsView());
            }
            if (key == "seed") {
                config.Seed = kv.as<size_t>();
            }
            if (key == "train") {
                training_range = Operon::ParseRange(value);
            }
            if (key == "test") {
                test_range = Operon::ParseRange(value);
            }
            if (key == "target") {
                target = value;
            }
            if (key == "maxlength") {
                max_length = kv.as<size_t>();
            }
            if (key == "maxdepth") {
                max_depth = kv.as<size_t>();
            }
            if (key == "enable-symbols") {
                auto mask = Operon::ParsePrimitiveSetConfig(value);
                pset_cfg |= mask;
            }
            if (key == "disable-symbols") {
                auto mask = ~Operon::ParsePrimitiveSetConfig(value);
                pset_cfg &= mask;
            }
            if (key == "threads") {
                threads = static_cast<decltype(threads)>(kv.as<size_t>());
            }
            if (key == "show-primitives") {
                show_pset = true;
            }
        }

        if (show_pset) {
            Operon::PrintPrimitives(pset_cfg);
            return 0;
        }
        if (auto res = dataset->GetVariable(target); !res.has_value()) {
            fmt::print(stderr, "error: target variable {} does not exist in the dataset.", target);
            return EXIT_FAILURE;
        }
        if (result.count("train") == 0) {
            training_range = Operon::Range{ 0, 2 * dataset->Rows() / 3 }; // by default use 66% of the data as training
        }
        if (result.count("test") == 0) {
            // if no test range is specified, we try to infer a reasonable range based on the training_range
            if (training_range.Start() > 0) {
                test_range = Operon::Range{ 0, training_range.Start() };
            } else if (training_range.End() < dataset->Rows()) {
                test_range = Operon::Range{ training_range.End(), dataset->Rows() };
            } else {
                test_range = Operon::Range{ 0, 1};
            }
        }
        // validate training range
        if (training_range.Start() >= dataset->Rows() || training_range.End() > dataset->Rows()) {
            fmt::print(stderr, "error: the training range {}:{} exceeds the available data range ({} rows)\n", training_range.Start(), training_range.End(), dataset->Rows());
            return EXIT_FAILURE;
        }

        if (training_range.Start() > training_range.End()) {
            fmt::print(stderr, "error: invalid training range {}:{}\n", training_range.Start(), training_range.End());
            return EXIT_FAILURE;
        }
        fmt::print("training range = {}:{}, test range = {}:{}\n", training_range.Start(), training_range.End(), test_range.Start(), test_range.End());

        std::vector<Operon::Variable> inputs;
        if (result.count("inputs") == 0) {
            auto variables = dataset->Variables();
            std::copy_if(variables.begin(), variables.end(), std::back_inserter(inputs), [&](auto const& var) { return var.Name != target; });
        } else {
            auto str = result["inputs"].as<std::string>();
            auto tokens = Operon::Split(str, ',');

            for (auto const& tok : tokens) {
                if (auto res = dataset->GetVariable(tok); res.has_value()) {
                    inputs.push_back(res.value());
                } else {
                    fmt::print(stderr, "error: variable {} does not exist in the dataset.", tok);
                    return EXIT_FAILURE;
                }
            }
        }

        auto problem = Operon::Problem(*dataset).Inputs(inputs).Target(target).TrainingRange(training_range).TestRange(test_range);
        problem.GetPrimitiveSet().SetConfig(pset_cfg);

        std::unique_ptr<Operon::CreatorBase> creator;
        creator = ParseCreator(result["tree-creator"].as<std::string>(), problem.GetPrimitiveSet(), problem.InputVariables());

        auto [amin, amax] = problem.GetPrimitiveSet().FunctionArityLimits();
        Operon::UniformTreeInitializer tree_init(*creator);
        tree_init.ParameterizeDistribution(amin+1, max_length);
        tree_init.SetMinDepth(1);
        tree_init.SetMaxDepth(1000); // NOLINT

        std::unique_ptr<Operon::CoefficientInitializerBase> coeff_init;
        std::unique_ptr<Operon::MutatorBase> one_point;
        if (symbolic) {
            using dist = std::uniform_int_distribution<int>;
            coeff_init = std::make_unique<Operon::CoefficientInitializer<dist>>();
            int constexpr range{5};
            dynamic_cast<Operon::CoefficientInitializer<dist>*>(coeff_init.get())->ParameterizeDistribution(-range, +range);
            one_point = std::make_unique<Operon::OnePointMutation<dist>>();
            dynamic_cast<Operon::OnePointMutation<dist>*>(one_point.get())->ParameterizeDistribution(-range, +range);
        } else {
            using dist = std::normal_distribution<Operon::Scalar>;
            coeff_init = std::make_unique<Operon::CoefficientInitializer<dist>>();
            dynamic_cast<Operon::NormalCoefficientInitializer*>(coeff_init.get())->ParameterizeDistribution(Operon::Scalar{0}, Operon::Scalar{1});
            one_point = std::make_unique<Operon::OnePointMutation<dist>>();
            dynamic_cast<Operon::OnePointMutation<dist>*>(one_point.get())->ParameterizeDistribution(Operon::Scalar{0}, Operon::Scalar{1});
        }

        auto const& [error, scale] = Operon::ParseErrorMetric(result["error-metric"].as<std::string>());
        Operon::Interpreter interpreter;

        // TODO: make these paths into cli options
        auto coord_path = result["coordinates"].as<std::string>();
        auto cluster_size = result["cluster-size"].as<size_t>();
        auto cutoff_radius = result["cutoff-radius"].as<Operon::Scalar>();


        atomic::summation_function sum(interpreter);
        sum.load_data(coord_path, static_cast<int>(cluster_size), cutoff_radius);

        Operon::RandomGenerator random(config.Seed);
        if (result["shuffle"].as<bool>()) {
            problem.GetDataset().Shuffle(random);
            auto idx = problem.GetDataset().GetValues("index");
            sum.set_index(idx.begin(), idx.end());
        }
        if (result["standardize"].as<bool>()) {
            problem.StandardizeData(problem.TrainingRange());
        }

        interpreter.GetDispatchTable().RegisterCallable(atomic::summation_function::hash, sum);

        Operon::Node sum_node(Operon::NodeType::Dynamic, atomic::summation_function::hash);
        sum_node.Arity = atomic::summation_function::arity;

        auto& pset = problem.GetPrimitiveSet();
        pset.AddPrimitive(sum_node, 1, 1, 1);
        //pset.SetMinMaxArity(static_cast<Operon::Hash>(Operon::NodeType::Div), 1, 1);

        Operon::Evaluator error_eval(problem, interpreter, *error, scale);
        error_eval.SetLocalOptimizationIterations(config.Iterations);
        error_eval.SetBudget(config.Evaluations);

        Operon::LengthEvaluator length_eval(problem, max_length);

        Operon::MultiEvaluator evaluator(problem);
        evaluator.SetBudget(config.Evaluations);
        evaluator.Add(error_eval);
        evaluator.Add(length_eval);

        EXPECT(problem.TestRange().Size() > 0);

        Operon::CrowdedComparison comp;

        auto female_selector = Operon::ParseSelector(result["female-selector"].as<std::string>(), comp);
        auto male_selector = Operon::ParseSelector(result["male-selector"].as<std::string>(), comp);


        Operon::SubtreeCrossover crossover{ crossover_internal_probability, max_depth, max_length };
        Operon::MultiMutation mutator{};

        Operon::ChangeVariableMutation change_var { problem.InputVariables() };
        Operon::ChangeFunctionMutation change_func { problem.GetPrimitiveSet() };
        Operon::ReplaceSubtreeMutation replace_subtree { *creator, *coeff_init, max_depth, max_length };
        Operon::InsertSubtreeMutation insert_subtree { *creator, *coeff_init, max_depth, max_length };
        Operon::RemoveSubtreeMutation remove_subtree { problem.GetPrimitiveSet() };
        Operon::DiscretePointMutation discrete_point;
        for (auto v : Operon::Math::Constants) {
            discrete_point.Add(static_cast<Operon::Scalar>(v), 1);
        }
        mutator.Add(*one_point, 1.0);
        mutator.Add(change_var, 1.0);
        mutator.Add(change_func, 1.0);
        mutator.Add(replace_subtree, 1.0);
        mutator.Add(insert_subtree, 1.0);
        mutator.Add(remove_subtree, 1.0);
        mutator.Add(discrete_point, 1.0);

        auto generator = Operon::ParseGenerator(result["offspring-generator"].as<std::string>(), evaluator, crossover, mutator, *female_selector, *male_selector);
        auto reinserter = Operon::ParseReinserter(result["reinserter"].as<std::string>(), comp);

        tf::Executor executor(threads);

        auto t0 = std::chrono::high_resolution_clock::now();
        Operon::RankIntersectSorter sorter;
        Operon::NSGA2 gp { problem, config, tree_init, *coeff_init, *generator, *reinserter, sorter };

        auto target_values = problem.TargetValues();
        auto target_train = target_values.subspan(training_range.Start(), training_range.Size());
        auto target_test = target_values.subspan(test_range.Start(), test_range.Size());

        // some boilerplate for reporting results
        const size_t idx { 0 };
        auto get_best = [&](Operon::Span<Operon::Individual const> pop) -> Operon::Individual {
            return *std::min_element(pop.begin(), pop.end(), [&](auto const& lhs, auto const& rhs) { return lhs[idx] < rhs[idx]; });
        };

        Operon::Individual best(1);
        //auto const& pop = gp.Parents();

        auto get_size = [](Operon::Individual const& ind) { return sizeof(ind) + sizeof(ind.Genotype) + sizeof(Operon::Node) * ind.Genotype.Nodes().capacity(); };

        tf::Executor exe(threads);

        auto report = [&]() {
            auto const& pop = gp.Parents();
            auto const& off = gp.Offspring();

            best = get_best(pop);

            Operon::Vector<Operon::Scalar> estimated_train;
            Operon::Vector<Operon::Scalar> estimated_test;

            tf::Taskflow taskflow;

            auto eval_train = taskflow.emplace([&]() {
                estimated_train = interpreter.Evaluate<Operon::Scalar>(best.Genotype, problem.GetDataset(), training_range);
            });

            auto eval_test = taskflow.emplace([&]() {
                estimated_test = interpreter.Evaluate<Operon::Scalar>(best.Genotype, problem.GetDataset(), test_range);
            });

            // scale values
            Operon::Scalar a{};
            Operon::Scalar b{};
            auto lst = taskflow.emplace([&]() {
                auto [a_, b_] = Operon::FitLeastSquares(estimated_train, target_train);
                a = static_cast<Operon::Scalar>(a_);
                b = static_cast<Operon::Scalar>(b_);
                // add scaling terms to the tree
                auto& nodes = best.Genotype.Nodes();
                auto const sz = nodes.size();
                if (std::abs(a - Operon::Scalar{1}) > std::numeric_limits<Operon::Scalar>::epsilon()) {
                    nodes.emplace_back(Operon::Node::Constant(a));
                    nodes.emplace_back(Operon::Node(Operon::NodeType::Mul));
                }
                if (std::abs(b) > std::numeric_limits<Operon::Scalar>::epsilon()) {
                    nodes.emplace_back(Operon::Node::Constant(b));
                    nodes.emplace_back(Operon::Node(Operon::NodeType::Add));
                }
                if (nodes.size() > sz) {
                    best.Genotype.UpdateNodes();
                }
            });

            double r2_train{};
            double r2_test{};
            double nmse_train{};
            double nmse_test{};
            double mae_train{};
            double mae_test{};

            auto scale_train = taskflow.emplace([&]() {
                Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>> estimated(estimated_train.data(), ssize(estimated_train));
                estimated = estimated * a + b;
            });

            auto scale_test = taskflow.emplace([&]() {
                Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>> estimated(estimated_test.data(), ssize(estimated_test));
                estimated = estimated * a + b;
            });

            auto calc_stats = taskflow.emplace([&]() {
                // negate the R2 because this is an internal fitness measure (minimization) which we here repurpose
                r2_train = -Operon::R2{}(estimated_train, target_train);
                r2_test = -Operon::R2{}(estimated_test, target_test);

                nmse_train = Operon::NMSE{}(estimated_train, target_train);
                nmse_test = Operon::NMSE{}(estimated_test, target_test);

                mae_train = Operon::MAE{}(estimated_train, target_train);
                mae_test = Operon::MAE{}(estimated_test, target_test);
            });

            double avg_length = 0;
            double avg_quality = 0;
            double total_memory = 0;

            auto calculate_length = taskflow.transform_reduce(pop.begin(), pop.end(), avg_length, std::plus<double>{}, [](auto const& ind) { return ind.Genotype.Length(); });
            auto calculate_quality = taskflow.transform_reduce(pop.begin(), pop.end(), avg_quality, std::plus<double>{}, [idx=idx](auto const& ind) { return ind[idx]; });
            auto calculate_pop_memory = taskflow.transform_reduce(pop.begin(), pop.end(), total_memory, std::plus{}, [get_size](auto const& ind) { return get_size(ind); });
            auto calculate_off_memory = taskflow.transform_reduce(off.begin(), off.end(), total_memory, std::plus{}, [get_size](auto const& ind) { return get_size(ind); });

            // define task graph
            lst.succeed(eval_train, eval_test);
            lst.precede(scale_train, scale_test);
            calc_stats.succeed(scale_train, scale_test);
            calc_stats.precede(calculate_length, calculate_quality, calculate_pop_memory, calculate_off_memory);

            exe.run(taskflow).wait();

            avg_length /= static_cast<double>(pop.size());
            avg_quality /= static_cast<double>(pop.size());

            auto t1 = std::chrono::high_resolution_clock::now();
            auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) / 1e6; // NOLINT

            using tup = std::tuple<std::string, double, std::string>;
            auto const* format = ":>#8.3g"; // see https://fmt.dev/latest/syntax.html
            std::array stats {
                tup{ "iteration", gp.Generation(), ":>" },
                tup{ "r2_tr", r2_train, format },
                tup{ "r2_te", r2_test, format },
                tup{ "mae_tr", mae_train, format },
                tup{ "mae_te", mae_test, format },
                tup{ "nmse_tr", nmse_train, format },
                tup{ "nmse_te", nmse_test, format },
                tup{ "avg_fit", avg_quality, format },
                tup{ "avg_len", avg_length, format },
                tup{ "eval_cnt", evaluator.CallCount , ":>" },
                tup{ "res_eval", evaluator.ResidualEvaluations, ":>" },
                tup{ "jac_eval", evaluator.JacobianEvaluations, ":>" },
                tup{ "seed", config.Seed, ":>" },
                tup{ "elapsed", elapsed, ":>"},
            };
            Operon::PrintStats({ stats.begin(), stats.end() }, gp.Generation() == 0);
        };

        gp.Run(executor, random, report);

        // write best values to file
        std::ofstream f(fmt::format("./results/{}.csv", config.Seed));
        f << "y_pred,y\n";
        auto y = problem.GetDataset().GetValues("energy");
        auto y_pred = interpreter.Evaluate<Operon::Scalar>(best.Genotype, problem.GetDataset(), {0UL, problem.GetDataset().Rows()});
        for (int i = 0UL; i < y.size(); ++i) {
            f << fmt::format("{},{}\n", y_pred[i], y[i]);
        }
        fmt::print("{}\n", Operon::InfixFormatter::Format(best.Genotype, problem.GetDataset(), 6)); // NOLINT
    } catch (std::exception& e) {
        fmt::print(stderr, "error: {}\n", e.what());
        return EXIT_FAILURE;
    }

    return 0;
}

