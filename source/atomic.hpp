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
Operon::Hash constexpr input_hash{9605156509256513176UL};

struct summation_function { // NOLINT

    static constexpr Operon::Hash hash{sum_hash};
    static constexpr size_t arity{1};

    summation_function(Operon::Interpreter& interpreter, Operon::Range range);

    auto load_data(std::string const&, int, int) -> void;

    // inner and outer cut-off radii, in Angstrom. Taken from Hernandez et al. 2019.
    static constexpr double inner_radius{3};
    static constexpr double outer_radius{5};

    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, size_t row) -> void
    {
        EXPECT(row < data_.size());

        auto& res = m[index];
        auto idx = static_cast<Eigen::Index>(index);

        // if this summation node is nested within another summation node, then behave like a passthrough
        if (std::any_of(nodes.begin() + idx + 1, nodes.end(), [](auto const& n) { return n.HashValue == atomic::summation_function::hash; })) {
            res = m[index-1];
            return;
        }

        // create a subtree skipping all previous summation nodes
        std::vector<Operon::Node> subtree;
        subtree.reserve(nodes[index].Length);
        auto it = nodes.begin() + idx;
        std::copy_if(it - nodes[idx].Length, it, std::back_inserter(subtree), [](auto const& n) { return n.HashValue != atomic::summation_function::hash; });
        auto tree = Operon::Tree(subtree).UpdateNodes();


        using S = typename std::remove_reference_t<decltype(res)>::Scalar; // NOLINT
        auto const r2in = inner_radius * inner_radius;
        auto const r2out = outer_radius * outer_radius;

        Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>> x(buf_.data(), buf_.size());

        auto remaining_rows = std::min(res.size(), Eigen::Index(range_.End() - row));
        for (auto r = 0; r < remaining_rows; ++r) {
            ENSURE(row + r < data_.size());

            interpreter_.get().Evaluate<Operon::Scalar>(tree, data_[row + r], {0UL, static_cast<size_t>(buf_.size())}, Operon::Span<Operon::Scalar>{buf_.data(), static_cast<size_t>(buf_.size())});

            //auto r2 = x*x;
            //res(r) = S{((2 * r2  - 3 * r2in + r2out) * (r2out - r2).pow(2) * std::pow(r2out - r2in, -3)).sum()};
            res(r) = S{x.sum()};
            //res(r) = S{std::reduce(x.begin(), x.end(), Operon::Scalar{0}, std::plus{})};
        }
    }

    private:
    std::reference_wrapper<Operon::Interpreter> interpreter_;
    Operon::Range range_;
    std::vector<Operon::Dataset> data_;
    std::vector<Operon::Scalar> buf_;

};
} // namespace atomic

#endif
