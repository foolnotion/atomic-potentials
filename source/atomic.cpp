#include "atomic.hpp"

#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <scn/scn.h>
#include <scn/scan/list.h>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <operon/core/dataset.hpp>

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
    coordinates.reserve(nline);

    Eigen::Index row { 0 };
    vector<Operon::Scalar> values;
    while (std::getline(f, line)) {
        if (row == 0) {
            coordinates.emplace_back(cluster_size, ncoord);
        }
        values.clear();
        scn::scan_list_options<char> opts{ {','}, {} };
        auto res = scn::scan_list_ex(line, values, opts);
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
        (void) scn::scan(line, "{}", e);
        energy.push_back(e);
    }
    return energy;
}

summation_function::summation_function(Operon::Interpreter& interpreter, Operon::Range range)
    : interpreter_(interpreter)
    , range_(std::move(range))
{
}

auto summation_function::load_data(std::string const& path, int cluster_size, int nearest_neighbors) -> void {
    auto coordinates = parse_coordinates(path, cluster_size);

    data_.clear();
    for (auto const& c : coordinates) {
        matrix dist(cluster_size, cluster_size);
        dist.matrix().diagonal().fill(0);

        for(auto i = 0; i < c.rows() - 1; ++i) { // NOLINT
            for (auto j = i+1; j < c.rows(); ++j) { // NOLINT
                dist(i, j) = dist(j, i) = (c.row(i) - c.row(j)).matrix().norm();
            }
        }

        std::vector<Operon::Scalar> distances;
        distances.reserve(size_t{static_cast<size_t>(nearest_neighbors * cluster_size)});

        for (auto col : dist.colwise()) {
            std::stable_sort(col.begin(), col.end());
            std::copy_n(col.begin(), nearest_neighbors, std::back_inserter(distances)); 
        }

        Eigen::Array<Operon::Scalar, -1, 1> m = Eigen::Map<decltype(m)>(distances.data(), ssize(distances));

        data_.emplace_back(m);
        data_.back().SetVariableNames({ "r" });
    }
    buf_.resize(data_[0].Rows());
    fmt::print("data size: {}\n", data_.size());
}

} // namespace atomic
