#ifndef ATOMIC_SUM_HPP
#define ATOMIC_SUM_HPP

#include <Eigen/Core>
#include <operon/interpreter/interpreter.hpp>
#include <ranges>
#include <string>
#include <vector>

#include <operon/core/dataset.hpp>
#include <operon/core/node.hpp>
#include <operon/core/types.hpp>
#include <operon/hash/hash.hpp>
#include <operon/formatter/formatter.hpp>

#include <operon/ceres/jet.h>
#include <vstat/vstat.hpp>


using matrix = Eigen::Array<Operon::Scalar, -1, -1>;

namespace atomic {

auto parse_coordinates(std::string const& coord_path, size_t cluster_size) -> std::vector<matrix>;
auto parse_energy(std::string const& energy_path) -> std::vector<Operon::Scalar>;

Operon::Hash constexpr sum_hash{1238717263861283123UL};
Operon::Hash constexpr pow_hash{9605156509256513176UL};

namespace vu = vstat::univariate;

template<typename Dist>
struct power_function {
    // a power function which draws a random exponent from a user-specified distribution
    static constexpr Operon::Hash hash{pow_hash};
    static constexpr size_t arity{1};

    power_function()
        : seed_(hash)
    {
    }

    power_function(power_function&&) noexcept = default;

    power_function(power_function const& other) {
        parameterize_distribution(other.param_);
        seed_.store(other.seed_.load());
    };

    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, Operon::Range /* unused */) -> void
    {
        Dist dist(param_);
        Operon::RandomGenerator random(seed_);
        using S = typename std::remove_reference_t<decltype(m[index])>::Scalar; // NOLINT
        auto exponent = S(static_cast<int64_t>(dist(random)));
        m[index] = m[index-1].pow(exponent);
        seed_.store(random());
    }

    template <typename... Args>
    auto parameterize_distribution(Args... args) const -> void
    {
        param_ = typename Dist::param_type { std::forward<Args&&>(args)... };
    }

    private:
    mutable std::atomic_ulong seed_{0};
    mutable typename Dist::param_type param_;
};

struct summation_function { // NOLINT
    static constexpr Operon::Hash hash{sum_hash};
    static constexpr size_t arity{1};

    explicit summation_function(Operon::Interpreter& interpreter)
        : interpreter_(interpreter) {}

    auto load_data(std::string const&, int, Operon::Scalar) -> void;

    // this function needs to be thread-safe!
    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, Operon::Range range) -> void
    {
        auto& res = m[index];
        auto is_sum = [](auto const& n) { return n.HashValue == atomic::summation_function::hash; };

        // if this summation node is nested within another summation node, then behave like a passthrough
        for (auto parent = nodes[index].Parent; parent != 0; parent = nodes[parent].Parent) {
            if (is_sum(nodes[parent])) {
                res = m[index-1];
                return;
            }
        }

        // create a subtree skipping all previous summation nodes
        auto it = nodes.begin() + static_cast<int64_t>(index);
        std::vector<Operon::Node> subtree;
        subtree.reserve(it->Length);
        std::copy_if(it - it->Length, it, std::back_inserter(subtree), std::not_fn(is_sum));
        auto tree = Operon::Tree(subtree).UpdateNodes();

        using S = typename std::remove_reference_t<decltype(res)>::Scalar; // NOLINT

        // allocate enough memory for all the evaluations (to avoid reallocating inside the loop)
        auto el = std::ranges::max_element(data_, std::less{}, &Operon::Dataset::Rows);
        std::unique_ptr<Operon::Scalar> buf(new (default_alignment) Operon::Scalar[el->Rows()]); // NOLINT

        for (auto row = 0; row < range.Size(); ++row) {
            ENSURE(range.Start() + row < data_.size());
            auto const& ds = data_[index_[range.Start() + row]];
            interpreter_.get().Evaluate<Operon::Scalar>(tree, ds, { 0UL, ds.Rows() }, { buf.get(), ds.Rows() });
            Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>> x(buf.get(), static_cast<Eigen::Index>(ds.Rows()));
            x *= ds.Values().col(2); // multiply with smoothing function
            res(row) = S{static_cast<Operon::Scalar>(vu::accumulate<Operon::Scalar>(buf.get(), ds.Rows()).sum)};
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

    private:
    std::reference_wrapper<Operon::Interpreter> interpreter_;
    std::vector<size_t> index_;
    std::vector<Operon::Dataset> data_;

    // inner and outer cut-off radii, in Angstrom. Taken from Hernandez et al. 2019.
    Operon::Scalar inner_radius_{default_inner_radius};
    Operon::Scalar outer_radius_{default_outer_radius};

    static constexpr Operon::Scalar default_inner_radius{3};
    static constexpr Operon::Scalar default_outer_radius{5};
    static constexpr std::align_val_t default_alignment{32};
};
} // namespace atomic

#endif
