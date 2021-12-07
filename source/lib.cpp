#include "lib.hpp"

#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <scn/scn.h>
#include <string>
#include <vector>

#include <fmt/format.h>

using std::getline;
using std::ifstream;
using std::string;
using std::vector;

constexpr size_t ncoord = 3;

namespace atomic {

auto parse_coordinates(string const& coord_path, size_t cluster_size) -> std::vector<matrix>
{
    std::ifstream f(coord_path);

    // count lines
    size_t nline = std::count(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>(), '\n');
    ENSURE(nline % cluster_size == 0);

    f.seekg(0); // rewind file stream
    string line;

    vector<matrix> coordinates;
    coordinates.reserve(nline / cluster_size);

    Eigen::Index row { 0 };
    vector<Operon::Scalar> values;
    while (std::getline(f, line)) {
        if (row == 0) {
            coordinates.emplace_back(cluster_size, ncoord);
        }
        values.clear();
        auto res = scn::scan_list(line, values, ',');
        ENSURE(res);
        coordinates.back().row(row++) = Eigen::Map<Eigen::Array<Operon::Scalar, -1, 1>>(values.data(), static_cast<Eigen::Index>(ncoord));
        if (row == cluster_size) {
            row = 0;
        }
    }

    return coordinates;
}

auto parse_energy(string const& energy_path) -> std::vector<Operon::Scalar>
{
    ifstream f(energy_path);
    vector<Operon::Scalar> energy;
    string line;
    Operon::Scalar e {};
    while (getline(f, line)) {
        scn::scan(line, "{}", e);
        energy.push_back(e);
    }
    return energy;
}

// a new function type to represent the summation of interatomic distances
summation::summation(string const& coord_path, string const& energy_path, size_t cluster_size)
    : cluster_coordinates_(parse_coordinates(coord_path, cluster_size))
    , cluster_energy_(parse_energy(energy_path))
    , cluster_size_(cluster_size)
{
    ENSURE(cluster_energy_.size() == cluster_coordinates_.size());

    for (auto const& m : cluster_coordinates_) {
        Operon::Scalar sum { 0 };
        matrix d(cluster_size, cluster_size);
        for (int i = 0; i < m.rows() - 1; ++i) {
            for (int j = i + 1; j < m.rows(); ++j) {
                d(i, j) = (m.row(i) - m.row(j)).matrix().norm();
            }
        }
        //std::cout << d << "\n";
        pairwise_distances_.push_back(d);
    }
}
} // namespace atomic
