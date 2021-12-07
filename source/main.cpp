// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: Copyright 2019-2021 Heal Research

#include <chrono>
#include <cstdlib>
#include <charconv>
#include <memory>
#include <thread>

#include <cxxopts.hpp>
#include <fmt/core.h>

#include <operon/algorithms/gp.hpp>
#include <operon/algorithms/nsga2.hpp>
#include <operon/core/format.hpp>
#include <operon/core/metrics.hpp>
#include <operon/operators/creator.hpp>
#include <operon/operators/crossover.hpp>
#include <operon/operators/evaluator.hpp>
#include <operon/operators/generator.hpp>
#include <operon/operators/initializer.hpp>
#include <operon/operators/mutation.hpp>
#include <operon/operators/reinserter.hpp>
#include <operon/operators/selector.hpp>

#include "lib.hpp"

using namespace Operon;

namespace util {
// parse a range specified as stard:end
auto parse_range(const std::string& range) -> Range
{
    auto pos = static_cast<int64_t>(range.find_first_of(':'));
    auto range_first = std::string(range.begin(), range.begin() + pos);
    auto range_last = std::string(range.begin() + pos + 1, range.end());
    size_t begin = 0;
    size_t end = 0;
    if (auto [p, ec] = std::from_chars(range_first.data(), range_first.data() + range_first.size(), begin); ec != std::errc()) {
        throw std::runtime_error(fmt::format("Could not parse training range from argument \"{}\"", range));
    }
    if (auto [p, ec] = std::from_chars(range_last.data(), range_last.data() + range_last.size(), end); ec != std::errc()) {
        throw std::runtime_error(fmt::format("Could not parse training range from argument \"{}\"", range));
    }
    return { begin, end };
}

// parses a double value
auto parse_double(const std::string& str)
{
    char* end;
    double val = std::strtod(str.data(), &end);
    return std::make_pair(val, std::isspace(*end) || end == str.data() + str.size());
}

// splits a string into substrings separated by delimiter
auto split(const std::string& s, char delimiter) -> std::vector<std::string>
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream token_stream(s);
    while (std::getline(token_stream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// formats a duration as dd:hh:mm:ss.ms
auto format_duration(std::chrono::duration<double> d) -> std::string
{
    auto h = std::chrono::duration_cast<std::chrono::hours>(d);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(d - h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(d - h - m);
    auto l = std::chrono::duration_cast<std::chrono::milliseconds>(d - h - m - s);
    return fmt::format("{:#02d}:{:#02d}:{:#02d}.{:#03d}", h.count(), m.count(), s.count(), l.count());
}

auto format_bytes(size_t bytes) -> std::string
{
    constexpr char sizes[] = " KMGT";
    auto p = static_cast<size_t>(std::floor(std::log2(bytes) / std::log2(1024)));
    return fmt::format("{:.2f} {}b", (double)bytes / std::pow(1024, p), sizes[p]);
}

static const std::unordered_map<std::string, NodeType> primitives {
    { "add",      NodeType::Add },
    { "mul",      NodeType::Mul },
    { "sub",      NodeType::Sub },
    { "div",      NodeType::Div },
    { "aq",       NodeType::Aq },
    { "pow",      NodeType::Pow },
    { "cbrt",     NodeType::Cbrt },
    { "cos",      NodeType::Cos },
    { "exp",      NodeType::Exp },
    { "log",      NodeType::Log },
    { "sin",      NodeType::Sin },
    { "sqrt",     NodeType::Sqrt },
    { "square",   NodeType::Square },
    { "tan",      NodeType::Tan },
    { "tanh",     NodeType::Tanh },
    { "dyn",      NodeType::Dynamic },
    { "constant", NodeType::Constant },
    { "variable", NodeType::Variable }
};

auto parse_primitive_set_config(const std::string& options) -> PrimitiveSetConfig
{
    auto config = static_cast<PrimitiveSetConfig>(0);
    for (auto& s : split(options, ',')) {
        if (auto it = primitives.find(s); it != primitives.end()) {
            config |= it->second;
        } else {
            fmt::print("Unrecognized symbol {}\n", s);
            std::exit(1);
        }
    }
    return config;
}
} // namespace util


auto main(int argc, char** argv) -> int
{
    cxxopts::Options opts("operon_cli", "C++ large-scale genetic programming");
    std::string symbols = "add, sub, mul, div, exp, log, square, sqrt, cbrt, sin, cos, tan, "
                          "asin, acos, atan, sinh, cosh, tanh, abs, aq, ceil, floor, fmin, fmax, "
                          "log1p";

    opts.add_options()
        ("dataset", "Dataset file name (csv) (required)", cxxopts::value<std::string>())
        ("shuffle", "Shuffle the input data", cxxopts::value<bool>()->default_value("false"))
        ("standardize", "Standardize the training partition (zero mean, unit variance)", cxxopts::value<bool>()->default_value("false"))
        ("train", "Training range specified as start:end (required)", cxxopts::value<std::string>())
        ("test", "Test range specified as start:end", cxxopts::value<std::string>())
        ("target", "Name of the target variable (required)", cxxopts::value<std::string>())
        ("inputs", "Comma-separated list of input variables", cxxopts::value<std::string>())
        ("error-metric", "The error metric used for calculating fitness", cxxopts::value<std::string>()->default_value("r2"))
        ("population-size", "Population size", cxxopts::value<size_t>()->default_value("1000"))
        ("pool-size", "Recombination pool size (how many generated offspring per generation)", cxxopts::value<size_t>()->default_value("1000"))
        ("seed", "Random number seed", cxxopts::value<Operon::RandomGenerator::result_type>()->default_value("0"))
        ("generations", "Number of generations", cxxopts::value<size_t>()->default_value("1000"))
        ("evaluations", "Evaluation budget", cxxopts::value<size_t>()->default_value("1000000"))
        ("iterations", "Local optimization iterations", cxxopts::value<size_t>()->default_value("0"))
        ("selection-pressure", "Selection pressure", cxxopts::value<size_t>()->default_value("100"))
        ("maxlength", "Maximum length", cxxopts::value<size_t>()->default_value("50"))
        ("maxdepth", "Maximum depth", cxxopts::value<size_t>()->default_value("10"))
        ("crossover-probability", "The probability to apply crossover", cxxopts::value<Operon::Scalar>()->default_value("1.0"))
        ("crossover-internal-probability", "Crossover bias towards swapping function nodes", cxxopts::value<Operon::Scalar>()->default_value("0.9"))
        ("mutation-probability", "The probability to apply mutation", cxxopts::value<Operon::Scalar>()->default_value("0.25"))
        ("tree-creator", "Tree creator operator to initialize the population with.", cxxopts::value<std::string>())
        ("female-selector", "Female selection operator, with optional parameters separated by : (eg, --selector tournament:5)", cxxopts::value<std::string>())
        ("male-selector", "Male selection operator, with optional parameters separated by : (eg, --selector tournament:5)", cxxopts::value<std::string>())
        ("offspring-generator", "OffspringGenerator operator, with optional parameters separated by : (eg --offspring-generator brood:10:10)", cxxopts::value<std::string>())
        ("reinserter", "Reinsertion operator merging offspring in the recombination pool back into the population", cxxopts::value<std::string>())
        ("enable-symbols", "Comma-separated list of enabled symbols ("+symbols+")", cxxopts::value<std::string>())
        ("disable-symbols", "Comma-separated list of disabled symbols ("+symbols+")", cxxopts::value<std::string>())
        ("show-primitives", "Display the primitive set used by the algorithm")
        ("threads", "Number of threads to use for parallelism", cxxopts::value<size_t>()->default_value("0"))
        ("timelimit", "Time limit after which the algorithm will terminate", cxxopts::value<size_t>()->default_value(fmt::format("{}", std::numeric_limits<size_t>::max())))
        ("debug", "Debug mode (more information displayed)")
        ("help", "Print help")
        ("version", "Print version and program information");

    auto result = opts.parse(argc, argv);
    if (result.arguments().empty() || result.count("help") > 0) {
        fmt::print("{}\n", opts.help());
        exit(EXIT_SUCCESS);
    }

    if (result.count("version") > 0) {
        fmt::print("{}\n", 1234);
        exit(EXIT_SUCCESS);
    }

    // parse and set default values
    GeneticAlgorithmConfig config{};
    config.Generations = result["generations"].as<size_t>();
    config.PopulationSize = result["population-size"].as<size_t>();
    config.PoolSize = result["pool-size"].as<size_t>();
    config.Evaluations = result["evaluations"].as<size_t>();
    config.Iterations = result["iterations"].as<size_t>();
    config.CrossoverProbability = result["crossover-probability"].as<Operon::Scalar>();
    config.MutationProbability = result["mutation-probability"].as<Operon::Scalar>();
    config.TimeLimit = result["timelimit"].as<size_t>();
    config.Seed = std::random_device {}();

    // parse remaining config options
    Range training_range;
    Range test_range;
    std::unique_ptr<Dataset> dataset;
    std::string target;
    bool showPrimitiveSet = false;
    auto threads = std::thread::hardware_concurrency();
    NodeType primitiveSetConfig = PrimitiveSet::Arithmetic;

    auto maxLength = result["maxlength"].as<size_t>();
    auto maxDepth = result["maxdepth"].as<size_t>();
    auto crossoverInternalProbability = result["crossover-internal-probability"].as<Operon::Scalar>();

    try {
        for (auto kv : result.arguments()) {
            auto& key = kv.key();
            auto& value = kv.value();

            if (key == "dataset") {
                dataset = std::make_unique<Dataset>(value, true);
                ENSURE(!dataset->IsView());
            }
            if (key == "seed") {
                config.Seed = kv.as<size_t>();
            }
            if (key == "train") {
                training_range = util::parse_range(value);
            }
            if (key == "test") {
                test_range = util::parse_range(value);
            }
            if (key == "target") {
                target = value;
            }
            if (key == "maxlength") {
                maxLength = kv.as<size_t>();
            }
            if (key == "maxdepth") {
                maxDepth = kv.as<size_t>();
            }
            if (key == "enable-symbols") {
                auto mask = util::parse_primitive_set_config(value);
                primitiveSetConfig |= mask;
            }
            if (key == "disable-symbols") {
                auto mask = ~util::parse_primitive_set_config(value);
                primitiveSetConfig &= mask;
            }
            if (key == "threads") {
                threads = static_cast<decltype(threads)>(kv.as<size_t>());
            }
            if (key == "show-primitives") {
                showPrimitiveSet = true;
            }
        }

        if (showPrimitiveSet) {
            PrimitiveSet tmpSet;
            tmpSet.SetConfig(primitiveSetConfig);
            fmt::print("Built-in primitives:\n");
            fmt::print("{:<8}\t{:<50}\t{:>7}\t\t{:>9}\n",
                "Symbol",
                "Description",
                "Enabled",
                "Frequency");
            for (size_t i = 0; i < NodeTypes::Count; ++i) {
                auto type = static_cast<NodeType>(1u << i);
                auto hash = Node(type).HashValue;
                auto enabled = tmpSet.Contains(hash) && tmpSet.IsEnabled(hash);
                auto freq = enabled ? tmpSet.Frequency(hash) : 0;
                Node node(type);
                fmt::print("{:<8}\t{:<50}\t{:>7}\t\t{:>9}\n",
                    node.Name(),
                    node.Desc(),
                    enabled,
                    freq ? std::to_string(freq) : "-");
            }
            return 0;
        }

        if (!dataset) {
            fmt::print(stderr, "{}\n{}\n", "Error: no dataset given.", opts.help());
            exit(EXIT_FAILURE);
        }
        if (result.count("target") == 0) {
            fmt::print(stderr,
                "{}\n{}\n",
                "Error: no target variable given.",
                opts.help());
            exit(EXIT_FAILURE);
        } else if (auto res = dataset->GetVariable(target); !res.has_value()) {
            fmt::print(stderr,
                "Target variable {} does not exist in the dataset.",
                target);
            exit(EXIT_FAILURE);
        }
        if (result.count("train") == 0) {
            training_range = {
                0,
                2 * dataset->Rows()
                    / 3
            }; // by default use 66% of the data as training
        }
        if (result.count("test") == 0) {
            // if no test range is specified, we try to infer a reasonable range
            // based on the trainingRange
            if (training_range.Start() > 0) {
                test_range = { 0, training_range.Start() };
            } else if (training_range.End() < dataset->Rows()) {
                test_range = { training_range.End(), dataset->Rows() };
            } else {
                test_range = { 0, 0 };
            }
        }
        // validate training range
        if (training_range.Start() >= dataset->Rows()
            || training_range.End() > dataset->Rows()) {
            fmt::print(stderr,
                "The training range {}:{} exceeds the available data range "
                "({} rows)\n",
                training_range.Start(),
                training_range.End(),
                dataset->Rows());
            exit(EXIT_FAILURE);
        }

        if (training_range.Start() > training_range.End()) {
            fmt::print(stderr,
                "Invalid training range {}:{}\n",
                training_range.Start(),
                training_range.End());
            exit(EXIT_FAILURE);
        }

        std::vector<Variable> inputs;
        if (result.count("inputs") == 0) {
            auto variables = dataset->Variables();
            std::copy_if(variables.begin(),
                variables.end(),
                std::back_inserter(inputs),
                [&](auto const& var) { return var.Name != target; });
        } else {
            auto str = result["inputs"].as<std::string>();
            auto tokens = util::split(str, ',');

            for (auto const& tok : tokens) {
                if (auto res = dataset->GetVariable(tok); res.has_value()) {
                    inputs.push_back(res.value());
                } else {
                    fmt::print(
                        stderr, "Variable {} does not exist in the dataset.", tok);
                    exit(EXIT_FAILURE);
                }
            }
        }

        auto problem = Problem(*dataset)
                           .Inputs(inputs)
                           .Target(target)
                           .TrainingRange(training_range)
                           .TestRange(test_range);
        problem.GetPrimitiveSet().SetConfig(primitiveSetConfig);

        using Reinserter = ReinserterBase;
        using OffspringGenerator = OffspringGeneratorBase;

        std::unique_ptr<Operon::CreatorBase> creator;

        if (result.count("tree-creator") == 0) {
            creator = std::make_unique<BalancedTreeCreator>(problem.GetPrimitiveSet(), problem.InputVariables(), 0.0);
        } else {
            auto value = result["tree-creator"].as<std::string>();

            if (value == "ptc2") {
                creator = std::make_unique<ProbabilisticTreeCreator>(problem.GetPrimitiveSet(), problem.InputVariables());
            } else if (value == "grow") {
                creator = std::make_unique<GrowTreeCreator>(problem.GetPrimitiveSet(), problem.InputVariables());
            } else {
                auto tokens = util::split(value, ':');
                double irregularity_bias = 0.0;
                if (tokens.size() > 1) {
                    if (auto [val, ok] = util::parse_double(tokens[1]); ok) {
                        irregularity_bias = val;
                    } else {
                        fmt::print(stderr,
                            "{}\n{}\n",
                            "Error: could not parse BTC bias argument.",
                            opts.help());
                        exit(EXIT_FAILURE);
                    }
                }
                creator = std::make_unique<BalancedTreeCreator>(problem.GetPrimitiveSet(), problem.InputVariables(), irregularity_bias);
            }
        }

        const size_t max_depth = 1000;
        auto [amin, amax] = problem.GetPrimitiveSet().FunctionArityLimits();
        UniformTreeInitializer initializer(*creator);
        initializer.ParameterizeDistribution(amin + 1, maxLength);
        initializer.SetMinDepth(1);
        initializer.SetMaxDepth(max_depth);

        UniformCoefficientInitializer coeff_init;
        initializer.ParameterizeDistribution(0, 1);

        auto crossover = SubtreeCrossover { crossoverInternalProbability, maxDepth, maxLength };
        auto mutator = MultiMutation {};
        auto one_point = OnePointMutation<std::uniform_real_distribution<Operon::Scalar>> {};
        one_point.ParameterizeDistribution(-2, 2);
        auto change_var = ChangeVariableMutation { problem.InputVariables() };
        auto change_func = ChangeFunctionMutation { problem.GetPrimitiveSet() };
        auto replace_subtree = ReplaceSubtreeMutation { *creator, maxDepth, maxLength };
        auto insert_subtree = InsertSubtreeMutation { *creator, maxDepth, maxLength, problem.GetPrimitiveSet() };
        auto remove_subtree = RemoveSubtreeMutation { problem.GetPrimitiveSet() };
        mutator.Add(one_point, 1.0);
        mutator.Add(change_var, 1.0);
        mutator.Add(change_func, 1.0);
        mutator.Add(replace_subtree, 1.0);
        mutator.Add(insert_subtree, 1.0);
        mutator.Add(remove_subtree, 1.0);


        constexpr size_t cluster_size = 32;
        //static const std::string coord_path { "./data/poet_run/data/DFT_Cu/NPT_FCC_1400K/F_coord_direct.csv" };
        static const std::string coord_path { "./data/poet_run/data/DFT_Cu/NPT_FCC_1400K/F_coord.data" };
        static const std::string energy_path { "./data/poet_run/data/DFT_Cu/NPT_FCC_1400K/E0.data" };
        atomic::summation sum(coord_path, energy_path, cluster_size);
        Node n(NodeType::Dynamic, sum.hash);
        n.Arity = 0;
        n.Length = 0;

        auto& pset = problem.GetPrimitiveSet();
        pset.AddPrimitive(n, /*frequency*/1, /*min arity*/0, /*max arity*/0);

        // initialize an interpreter using a default dispatch table
        // (the dispatch table is by default initialized with the common
        // operations, but they can be overriden)
        Interpreter::DT ft;
        ft.RegisterCallable(sum.hash, sum);
        Interpreter interpreter(ft);

        std::unique_ptr<EvaluatorBase> error_evaluator;
        auto error_metric = result["error-metric"].as<std::string>();
        if (error_metric == "r2") {
            error_evaluator = std::make_unique<RSquaredEvaluator>(problem, interpreter);
        } else if (error_metric == "nmse") {
            error_evaluator = std::make_unique<NormalizedMeanSquaredErrorEvaluator>(problem, interpreter);
        } else if (error_metric == "mse") {
            error_evaluator = std::make_unique<MeanSquaredErrorEvaluator>(problem, interpreter);
        } else if (error_metric == "rmse") {
            error_evaluator = std::make_unique<RootMeanSquaredErrorEvaluator>(problem, interpreter);
        } else if (error_metric == "mae") {
            error_evaluator = std::make_unique<MeanAbsoluteErrorEvaluator>(problem, interpreter);
        } else if (error_metric == "l2") {
            error_evaluator = std::make_unique<L2NormEvaluator>(problem, interpreter);
        } else {
            throw std::runtime_error(
                fmt::format("Unknown metric {}\n", error_metric));
        }
        error_evaluator->SetLocalOptimizationIterations(config.Iterations);

        UserDefinedEvaluator length_evaluator(
            problem,
            [](Operon::RandomGenerator& /*unused*/, Individual& ind) {
                return Operon::Vector<Operon::Scalar> {
                    static_cast<Operon::Scalar>(ind.Genotype.Length())
                };
            });
        UserDefinedEvaluator shape_evaluator(
            problem,
            [](Operon::RandomGenerator& /*unused*/, Individual& ind) {
                return Operon::Vector<Operon::Scalar> {
                    static_cast<Operon::Scalar>(ind.Genotype.VisitationLength())
                };
            });

        std::unique_ptr<MultiEvaluator> evaluator(new MultiEvaluator(problem));
        evaluator->SetBudget(config.Evaluations);
        evaluator->Add(*error_evaluator);
        evaluator->Add(length_evaluator);
        evaluator->Add(shape_evaluator);

        EXPECT(problem.TrainingRange().Size() > 0);

        // auto comp = [](auto const& lhs, auto const& rhs) { return lhs[0] <
        // rhs[0]; };
        CrowdedComparison comp;

        auto parse_selector = [&](const std::string& name) -> SelectorBase* {
            constexpr size_t default_tournament_size{5};

            if (result.count(name) == 0) {
                auto *sel = new TournamentSelector(comp);
                sel->SetTournamentSize(default_tournament_size);
                return sel;
            }
            auto value = result[name].as<std::string>();
            auto tokens = util::split(value, ':');
            if (tokens[0] == "tournament") {
                size_t t_size = default_tournament_size;
                if (tokens.size() > 1) {
                    if (auto [p, ec] = std::from_chars(tokens[1].data(),
                            tokens[1].data() + tokens[1].size(),
                            t_size);
                        ec != std::errc()) {
                        fmt::print(stderr,
                            "{}\n{}\n",
                            "Error: could not parse tournament size argument.",
                            opts.help());
                        exit(EXIT_FAILURE);
                    }
                }
                auto *sel = new TournamentSelector(comp);
                sel->SetTournamentSize(t_size);
                return sel;
            }
            if (tokens[0] == "proportional") {
                auto *sel = new ProportionalSelector(comp);
                sel->SetObjIndex(0);
                return sel;
            }
            if (tokens[0] == "rank") {
                size_t t_size = default_tournament_size;
                if (tokens.size() > 1) {
                    if (auto [p, ec] = std::from_chars(tokens[1].data(),
                            tokens[1].data() + tokens[1].size(),
                            t_size);
                        ec != std::errc()) {
                        fmt::print(stderr,
                            "{}\n{}\n",
                            "Error: could not parse tournament size argument.",
                            opts.help());
                        exit(EXIT_FAILURE);
                    }
                }
                auto *sel = new RankTournamentSelector(comp);
                sel->SetTournamentSize(t_size);
                return sel;
            }
            if (tokens[0] == "random") {
                return new RandomSelector();
            }
            auto *sel = new TournamentSelector(comp);
            sel->SetTournamentSize(default_tournament_size);
            return sel;
        };

        std::unique_ptr<SelectorBase> female_selector;
        std::unique_ptr<SelectorBase> male_selector;

        female_selector.reset(parse_selector("female-selector"));
        male_selector.reset(parse_selector("male-selector"));

        std::unique_ptr<OffspringGenerator> generator;
        if (result.count("offspring-generator") == 0) {
                generator = std::make_unique<BasicOffspringGenerator>(*evaluator, crossover, mutator, *female_selector, *male_selector);
        } else {
            auto value = result["offspring-generator"].as<std::string>();
            auto tokens = util::split(value, ':');
            if (tokens[0] == "basic") {
                generator = std::make_unique<BasicOffspringGenerator>(*evaluator, crossover, mutator, *female_selector, *male_selector);
            } else if (tokens[0] == "brood") {
                size_t brood_size = 10;
                if (tokens.size() > 1) {
                    if (auto [p, ec] = std::from_chars(tokens[1].data(),
                            tokens[1].data() + tokens[1].size(),
                            brood_size);
                        ec != std::errc()) {
                        fmt::print(stderr,
                            "{}\n{}\n",
                            "Error: could not parse brood size argument.",
                            opts.help());
                        exit(EXIT_FAILURE);
                    }
                }
                generator = std::make_unique<BroodOffspringGenerator>(*evaluator, crossover, mutator, *female_selector, *male_selector);
                auto *ptr = dynamic_cast<BroodOffspringGenerator*>(generator.get());
                ptr->BroodSize(brood_size);
            } else if (tokens[0] == "poly") {
                size_t brood_size = 10;
                if (tokens.size() > 1) {
                    if (auto [p, ec] = std::from_chars(tokens[1].data(),
                            tokens[1].data() + tokens[1].size(),
                            brood_size);
                        ec != std::errc()) {
                        fmt::print(stderr,
                            "{}\n{}\n",
                            "Error: could not parse brood size argument.",
                            opts.help());
                        exit(EXIT_FAILURE);
                    }
                }
                generator = std::make_unique<PolygenicOffspringGenerator>(*evaluator, crossover, mutator, *female_selector, *male_selector);
                auto *ptr = dynamic_cast<PolygenicOffspringGenerator*>(generator.get());
                ptr->PolygenicSize(brood_size);
            } else if (tokens[0] == "os") {
                size_t selection_pressure = 100;
                double comparison_factor = 1.0;
                if (tokens.size() > 1) {
                    if (auto [p, ec] = std::from_chars(tokens[1].data(),
                            tokens[1].data() + tokens[1].size(),
                            selection_pressure);
                        ec != std::errc()) {
                        fmt::print(
                            stderr,
                            "{}\n{}\n",
                            "Error: could not parse maximum selection pressure argument.",
                            opts.help());
                        exit(EXIT_FAILURE);
                    }
                }
                if (tokens.size() > 2) {
                    if (auto [val, ok] = util::parse_double(tokens[2]); ok) {
                        comparison_factor = val;
                    } else {
                        fmt::print(stderr,
                            "{}\n{}\n",
                            "Error: could not parse comparison factor argument.",
                            opts.help());
                        exit(EXIT_FAILURE);
                    }
                }
                generator = std::make_unique<OffspringSelectionGenerator>(*evaluator, crossover, mutator, *female_selector, *male_selector);
                auto *ptr = dynamic_cast<OffspringSelectionGenerator*>(generator.get());
                ptr->MaxSelectionPressure(selection_pressure);
                ptr->ComparisonFactor(comparison_factor);
            }
        }
        std::unique_ptr<Reinserter> reinserter;
        if (result.count("reinserter") == 0) {
            reinserter = std::make_unique<ReplaceWorstReinserter>(comp);
        } else {
            auto value = result["reinserter"].as<std::string>();
            if (value == "keep-best") {
                reinserter = std::make_unique<KeepBestReinserter>(comp);
            } else if (value == "replace-worst") {
                reinserter = std::make_unique<ReplaceWorstReinserter>(comp);
            }
        }

        Operon::RandomGenerator random(config.Seed);
        if (result["shuffle"].as<bool>()) {
            problem.GetDataset().Shuffle(random);
        }
        if (result["standardize"].as<bool>()) {
            problem.StandardizeData(problem.TrainingRange());
        }

        tf::Executor executor(threads);

        auto t0 = std::chrono::high_resolution_clock::now();
        RankSorter sorter;

        NSGA2 gp { problem, config, initializer, coeff_init, *generator, *reinserter, sorter };

        auto target_values = problem.TargetValues();
        auto target_train = target_values.subspan(training_range.Start(), training_range.Size());
        auto target_test = target_values.subspan(test_range.Start(), test_range.Size());

        // some boilerplate for reporting results
        const size_t idx { 0 };
        auto get_best = [&](const Operon::Span<const Individual> pop) {
            const auto *min_elem = std::min_element(pop.begin(),
                pop.end(),
                [&](auto const& lhs, auto const& rhs) { return lhs[idx] < rhs[idx]; });
            ENSURE(min_elem->Genotype.Length() > 0);
            ENSURE(min_elem->Fitness.size() > 1);
            return *min_elem;
        };

        Individual best;
        // auto const& pop = gp.Parents();

        auto get_size = [](const Individual& ind) {
            return sizeof(ind) + sizeof(ind.Genotype)
                + sizeof(Node) * ind.Genotype.Nodes().capacity();
        };

        tf::Executor exe(threads);

        auto report = [&]() {
            auto pop = gp.Parents();
            auto off = gp.Offspring();


            auto contains_sum = [&](auto const& tree) {
                auto const& nodes = tree.Nodes();
                bool res = std::any_of(nodes.begin(), nodes.end(), [&](auto const& node) { return node.HashValue == sum.hash && node.Arity == 0; });
                //if (res) {
                //    fmt::print("tree: {}\n", InfixFormatter::Format(tree, problem.GetDataset(), 20));
                //}
                return res;
            };

            bool has_sum = std::any_of(pop.begin(), pop.end(), [&](auto const& ind) { return contains_sum(ind.Genotype); });

            if (!has_sum) {
                throw std::runtime_error("the population does not contain summation symbols\n");
            }

            auto best_front= std::vector<Individual>(gp.Best().begin(), gp.Best().end());
            std::sort(best_front.begin(),
                best_front.end(),
                [&](auto const& lhs, auto const& rhs) {
                    EXPECT(lhs.Fitness.size() > 1);
                    EXPECT(rhs.Fitness.size() > 1);
                    return lhs[idx] < rhs[idx];
                });

            best = get_best(gp.Best());

            Operon::Vector<Operon::Scalar> estimated_train;
            Operon::Vector<Operon::Scalar> estimated_test;

            tf::Taskflow taskflow;

            auto eval_train = taskflow.emplace(
                [&]() {
                    estimated_train = interpreter.Evaluate<Operon::Scalar>(
                        best.Genotype, problem.GetDataset(), training_range);
                });

            auto eval_test = taskflow.emplace(
                [&]() {
                    estimated_test = interpreter.Evaluate<Operon::Scalar>(
                        best.Genotype, problem.GetDataset(), test_range);
                });

            // scale values
            Operon::Scalar a = 1;
            Operon::Scalar b = 0;
            auto linear_scaling = taskflow.emplace(
                [&]() {
                    auto stats = bivariate::accumulate<double>(estimated_train.data(),
                        target_train.data(),
                        estimated_train.size());
                    a = static_cast<Operon::Scalar>(stats.covariance
                        / stats.variance_x);
                    b = static_cast<Operon::Scalar>(stats.mean_y - a * stats.mean_x);
                });

            double r2_train;
            double r2_test;
            double nmse_train;
            double nmse_test;
            double mae_train;
            double mae_test;

            auto scale_train = taskflow.emplace(
                [&]() {
                    Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>> estimated(
                        estimated_train.data(), estimated_train.size());
                    estimated = estimated * a + b;
                });

            auto scale_test = taskflow.emplace(
                [&]() {
                    Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>> estimated(
                        estimated_test.data(), estimated_test.size());
                    estimated = estimated * a + b;
                });

            auto calc_stats = taskflow.emplace(
                [&]() {
                    r2_train = RSquared<Operon::Scalar>(estimated_train, target_train);
                    r2_test = RSquared<Operon::Scalar>(estimated_test, target_test);

                    nmse_train = NormalizedMeanSquaredError<Operon::Scalar>(
                        estimated_train, target_train);
                    nmse_test = NormalizedMeanSquaredError<Operon::Scalar>(
                        estimated_test, target_test);

                    mae_train = MeanAbsoluteError<Operon::Scalar>(estimated_train,
                        target_train);
                    mae_test = MeanAbsoluteError<Operon::Scalar>(estimated_test, target_test);
                });

            double avg_length = 0;
            double avg_quality = 0;
            double total_memory = 0;

            EXPECT(std::all_of(pop.begin(),
                pop.end(),
                [](auto const& ind) { return ind.Genotype.Length() > 0; }));

            auto calculate_length = taskflow.transform_reduce(
                pop.begin(),
                pop.end(),
                avg_length,
                std::plus<double> {},
                [](auto const& ind) { return ind.Genotype.Length(); });
            auto calculate_quality = taskflow.transform_reduce(
                pop.begin(),
                pop.end(),
                avg_quality,
                std::plus<double> {},
                [idx = idx](auto const& ind) { return ind[idx]; });
            auto calculate_pop_memory = taskflow.transform_reduce(
                pop.begin(),
                pop.end(),
                total_memory,
                std::plus {},
                [get_size](auto const& ind) { return get_size(ind); });
            auto calculate_off_memory = taskflow.transform_reduce(
                off.begin(),
                off.end(),
                total_memory,
                std::plus {},
                [get_size](auto const& ind) { return get_size(ind); });

            // define task graph
            linear_scaling.succeed(eval_train, eval_test);
            linear_scaling.precede(scale_train, scale_test);
            calc_stats.succeed(scale_train, scale_test);
            calc_stats.precede(calculate_length, calculate_quality, calculate_pop_memory, calculate_off_memory);

            exe.run(taskflow).wait();

            avg_length /= static_cast<double>(pop.size());
            avg_quality /= static_cast<double>(pop.size());

            auto t1 = std::chrono::high_resolution_clock::now();
            auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1- t0).count()) / 1e6;

            fmt::print("{:.4f}\t{}\t", elapsed, gp.Generation());
            fmt::print("{:.4f}\t{:.4f}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}\t",
                r2_train,
                r2_test,
                mae_train,
                mae_test,
                nmse_train,
                nmse_test);
            fmt::print("{:.4g}\t{:.1f}\t{:.3f}\t{:.3f}\t{}\t{}\t{}\t",
                avg_quality,
                avg_length,
                /* shape */ 0.0,
                /* diversity */ 0.0,
                evaluator->FitnessEvaluations(),
                evaluator->LocalEvaluations(),
                evaluator->TotalEvaluations());
            fmt::print("{}\t{}\n", total_memory, config.Seed);
        };

        gp.Run(executor, random, report);
        fmt::print("{}\n", InfixFormatter::Format(best.Genotype, problem.GetDataset(), 20));
    } catch (std::exception& e) {
        fmt::print("{}\n", e.what());
        std::exit(EXIT_FAILURE);
    }

    return 0;
}
