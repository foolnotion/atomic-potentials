#ifndef ATOMIC_SUM_HPP
#define ATOMIC_SUM_HPP

#include <Eigen/Core>
#include <operon/interpreter/interpreter.hpp>
#include <string>
#include <vector>

#include <operon/core/dataset.hpp>
#include <operon/core/node.hpp>
#include <operon/core/types.hpp>
#include <operon/hash/hash.hpp>
#include <operon/core/format.hpp>

#include <operon/ceres/jet.h>


using matrix = Eigen::Array<Operon::Scalar, -1, -1>;

namespace atomic {

auto parse_coordinates(std::string const& coord_path, size_t cluster_size) -> std::vector<matrix>;
auto parse_energy(std::string const& energy_path) -> std::vector<Operon::Scalar>;

Operon::Hash constexpr sum_hash{1238717263861283123UL};
Operon::Hash constexpr smooth_hash{8127387129837126261UL};
Operon::Hash constexpr input_hash{9605156509256513176UL};

struct smooth { // NOLINT
    static constexpr Operon::Hash hash{smooth_hash};

    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, size_t row) -> void
    {
        auto& res = m[index];
        auto& arg = m[index-1];

        using S = typename std::remove_reference_t<decltype(res)>::Scalar; // NOLINT

        auto r2 = arg * arg;
        auto r2in = S{inner_radius * inner_radius};
        auto r2out = S{outer_radius * outer_radius};

        res = ((S{2} * r2 - S{3} * r2in + r2out) * (r2out - r2).pow(S{2}) * ceres::pow(r2out - r2in, S{-3}));
    }

    static constexpr size_t arity{1};

    // inner and outer cut-off radii, in Angstrom. Taken from Hernandez et al. 2019.
    static constexpr double inner_radius{3};
    static constexpr double outer_radius{5};
};

struct sum { // NOLINT
    static constexpr Operon::Hash hash{sum_hash};
    static constexpr size_t arity{1};

    sum(Operon::Interpreter& interpreter, std::string const&, size_t);

    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, size_t row) -> void
    {
        auto& res = m[index];
        auto idx = static_cast<Eigen::Index>(index);

        // if this summation node is nested within another summation node, then behave like a passthrough
        if (std::any_of(nodes.begin() + idx + 1, nodes.end(), [](auto const& n) { return n.HashValue == sum_hash; })) {
            res = m[index-1];
            return;
        }

        // create a subtree skipping all previous summation nodes
        std::vector<Operon::Node> subtree;
        subtree.reserve(nodes[index].Length);

        for (auto i = idx - nodes[idx].Length; i < idx; ++i) {
            if (nodes[i].HashValue == sum_hash) { continue; }
            subtree.push_back(nodes[i]);
        }
        auto tree = Operon::Tree(subtree).UpdateNodes();
        std::vector<Operon::Scalar> result(data_[0].Rows());

        using S = typename std::remove_reference_t<decltype(res)>::Scalar; // NOLINT

        for (auto r = row; r < std::min(data_.size(), row + res.size()); ++r) {
            interpreter_.get().Evaluate(tree, data_[r], Operon::Range(0, data_[r].Rows()), Operon::Span<Operon::Scalar> { result.data(), result.size() }); 
            res(r - row) = S{Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>>(result.data(), static_cast<Eigen::Index>(result.size())).sum()};
        }
    }

    private:
    std::vector<Operon::Dataset> data_;
    std::reference_wrapper<Operon::Interpreter> interpreter_;
    size_t cluster_size_;
};
} // namespace atomic

#endif
