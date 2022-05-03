#include "atomic.hpp"

#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <scn/scn.h>
#include <scn/scan/list.h>
#include <string>
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

summation_function::summation_function(Operon::Interpreter& interpreter, Operon::Range range, /* path to coordinates file */ std::string const& cpath, /*cluster_size*/ size_t size)
    : interpreter_(interpreter), cluster_size_(size)
    , range_(range)
{
    auto coordinates = parse_coordinates(cpath, size);

    auto s = static_cast<Eigen::Index>(size);
    for (auto const& c : coordinates) {
        matrix dist(s, s);
        dist.matrix().diagonal().fill(0);

        for(auto i = 0; i < c.rows() - 1; ++i) { // NOLINT
            for (auto j = i+1; j < c.rows(); ++j) { // NOLINT
                dist(i, j) = dist(j, i) = (c.row(i) - c.row(j)).matrix().norm();
            }
        }

        Eigen::Array<Operon::Scalar, -1, 1> m = dist.reshaped();
        data_.emplace_back(m);
        data_.back().SetVariableNames({ "r" });
    }
    fmt::print("data size: {}, input size: {}\n", data_.size(), data_.front().Rows());
}
} // namespace atomic
