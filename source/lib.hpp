#ifndef ATOMIC_SUM_HPP
#define ATOMIC_SUM_HPP

#include <string>
#include <vector>
#include <Eigen/Core>

#include <operon/core/types.hpp>
#include <operon/hash/hash.hpp>
#include <operon/core/node.hpp>

using matrix = Eigen::Array<Operon::Scalar, -1, -1>;

namespace atomic {

auto parse_coordinates(std::string const& coord_path, size_t cluster_size) -> std::vector<matrix>;
auto parse_energy(std::string const& energy_path) -> std::vector<Operon::Scalar>;

struct summation {
    static constexpr std::array name{"summation"};
    Operon::Hash const hash = Operon::Hasher{}(reinterpret_cast<uint8_t const*>(name.data()), name.size()); // NOLINT
    static constexpr int arity = 0;
    static constexpr int length = 0;

    summation(std::string const& coord_path, std::string const& energy_path, size_t cluster_size);

    template<typename T, typename U>
    auto smooth(T& r, U r_out) {
        U r_in{0}; // inner radius
        auto r2 = r * r;
        auto r2_in = r_in * r_in;
        auto r2_out = r_out * r_out;
        auto r2_out_in = r2_out - r2_in;
        return (U{2} * r2 - U{3} * r2_in + r2_out) * (r2_out - r2) * r2_out_in * r2_out_in * r2_out_in;
    }

    // here, the typename T represents the Eigen type (typically, an Eigen::Array) used to hold
    // the intermediate node evaluation values
    // typename T::Scalar will represent the actual number type (a scalar or a dual)
    template<typename T>
    auto operator()(T& m, std::vector<Operon::Node> const& nodes, size_t index, size_t row) -> void {
        auto r = m.col(index);
        //auto radius = typename T::Scalar { nodes[index].Value };
        auto radius = nodes[index].Value;
        ENSURE(nodes[index].HashValue == hash);

        auto remaining_rows = static_cast<Eigen::Index>(cluster_energy_.size() - row);
        //fmt::print("m.rows() = {}, remaining_rows = {}\n", m.rows(), remaining_rows);
        for (size_t i = 0; i < std::min(m.rows(), remaining_rows); ++i) {
            auto const& d = pairwise_distances_[i];
            r(i) = typename T::Scalar{(d > radius).select(0, d).sum()};
            //r(i) = typename T::Scalar{d.sum()};
            //r(i) = typename T::Scalar{smooth(d, radius).sum()};
        }
        //std::cout << r << "\n";
    }

    private:
        std::vector<matrix> cluster_coordinates_;
        std::vector<matrix> pairwise_distances_;
        std::vector<Operon::Scalar> cluster_energy_;
        size_t cluster_size_;
};

} // namespace atomic

#endif
