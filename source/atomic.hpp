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
#include <vstat/vstat.hpp>


using matrix = Eigen::Array<Operon::Scalar, -1, -1>;

namespace atomic {

auto parse_coordinates(std::string const& coord_path, size_t cluster_size) -> std::vector<matrix>;
auto parse_energy(std::string const& energy_path) -> std::vector<Operon::Scalar>;

Operon::Hash constexpr sum_hash{1238717263861283123UL};
Operon::Hash constexpr input_hash{9605156509256513176UL};

namespace vu = vstat::univariate;

struct summation_function { // NOLINT
    static constexpr Operon::Hash hash{sum_hash};
    static constexpr size_t arity{1};

    explicit summation_function(Operon::Interpreter& interpreter)
        : interpreter_(interpreter) {}

    auto load_data(std::string const&, int, double) -> void;

    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, Operon::Range range) -> void
    {
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
        auto const r2in = static_cast<Operon::Scalar>(inner_radius_ * inner_radius_);
        auto const r2out = static_cast<Operon::Scalar>(outer_radius_ * outer_radius_);

        // smoothing function - see Hernandez 2019 https://www.nature.com/articles/s41524-019-0249-1/
        auto f = [&](auto r)  {
            auto r2 = r * r;
            return (2 * r2 - 3 * r2in + r2out) * (r2out - r2) * (r2out - r2) * std::pow(r2out - r2in, -3);
        };

        // allocate enough memory for all the evaluations (to avoid reallocating inside the loop)
        std::vector<Operon::Scalar> buf(std::max_element(data_.begin(), data_.end(), [](auto const& a, auto const& b) { return a.Rows() < b.Rows(); })->Rows());

        for (auto row = 0; row < range.Size(); ++row) {
            ENSURE(range.Start() + row < data_.size());
            auto const& ds = data_[index_[range.Start() + row]];
            Operon::Range rg{0UL, ds.Rows()};

            interpreter_.get().Evaluate<Operon::Scalar>(tree, ds, rg, { buf.data(), ds.Rows() });
            Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>> x(buf.data(), ds.Rows());
            //res(row) = S{x.unaryExpr([&](auto r) { return r * f(r); }).sum()};
            //res(row) = S{static_cast<Operon::Scalar>(vstat::univariate::accumulate<Operon::Scalar>(buf_.data(), buf_.size()).sum)};
            res(row) = S{static_cast<Operon::Scalar>(vu::accumulate<Operon::Scalar>(buf.data(), ds.Rows(), [&](auto r) { return r * f(r); }).sum)};
        }
    }

    template<typename T>
    requires std::is_integral_v<T>
    void set_index(std::vector<T> const& index) {
        if constexpr (std::is_same_v<T, size_t>) {
            index_ = index;
        } else {
            std::copy(index.begin(), index.end(), index_.begin());
        }
    }

    template<typename Iter>
    void set_index(Iter begin, Iter end) {
        auto it = index_.begin();
        while (begin != end) {
            *it++ = static_cast<size_t>(*begin++);
        }
    }

    [[nodiscard]] auto index() const -> Operon::Span<size_t const> { return {index_.data(), index_.size()}; }

    void apply_permutation(std::vector<Eigen::Index> const& indices) {
        decltype(data_) tmp; tmp.reserve(data_.size());
        std::transform(indices.begin(), indices.end(), std::back_inserter(tmp), [&](auto i) { return std::move(data_[indices[i]]); });
        std::swap(data_, tmp);
    }

    private:
    std::reference_wrapper<Operon::Interpreter> interpreter_;
    std::vector<size_t> index_;
    std::vector<Operon::Dataset> data_;

    // inner and outer cut-off radii, in Angstrom. Taken from Hernandez et al. 2019.
    double inner_radius_{3};
    double outer_radius_{5};
};
} // namespace atomic

#endif
