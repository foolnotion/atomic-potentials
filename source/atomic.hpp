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

struct smoothing_function { // NOLINT
    static constexpr Operon::Hash hash{smooth_hash};

    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, size_t row) -> void
    {
        auto& res = m[index];
        auto& arg = m[index-1];

        // act as passthrough if not under a summation
        if (std::none_of(nodes.begin() + index + 1, nodes.end(), [](auto const& n) { return n.HashValue == sum_hash; })) {
            res = arg;
            return;
        }

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

struct summation_function { // NOLINT
    static constexpr Operon::Hash hash{sum_hash};
    static constexpr size_t arity{1};

    summation_function(Operon::Interpreter& interpreter, std::string const&, size_t);

    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, size_t row) -> void
    {
        EXPECT(row < data_.size());

        auto& res = m[index];
        auto idx = static_cast<Eigen::Index>(index);

        // if this summation node is nested within another summation node, then behave like a passthrough
        if (std::any_of(nodes.begin() + idx + 1, nodes.end(), [](auto const& n) { return n.HashValue == sum_hash; })) {
            res = m[index-1];
            return;
        }

        std::vector<Operon::Scalar> result(data_[0].Rows());
        Operon::Range range(0, result.size());
        Operon::Span<Operon::Scalar> span{result.data(), result.size()};

        // create a subtree skipping all previous summation nodes
        std::vector<Operon::Node> subtree;
        subtree.reserve(nodes[index].Length);
        auto it = nodes.begin() + idx;
        std::copy_if(it - nodes[idx].Length, it, std::back_inserter(subtree), [](auto const& n) { return n.HashValue != atomic::summation_function::hash; });
        auto tree = Operon::Tree(subtree).UpdateNodes();

        using S = typename std::remove_reference_t<decltype(res)>::Scalar; // NOLINT

        auto p = std::min(res.size(), static_cast<int64_t>(data_.size() - row));
        for (auto r = row; r < row + p; ++r) {
            interpreter_.get().Evaluate(tree, data_[r], range, span); 
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
