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

using Eigen::Array;
using Eigen::Map;

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

// this function loads the snaphots of atomic positions into separate datasets (one per snapshot)
// only distances below the outer cutoff radius are added
auto summation_function::load_data(std::string const& path, int cluster_size, Operon::Scalar cutoff_radius) -> void {
    auto coordinates = parse_coordinates(path, cluster_size);
    outer_radius_ = cutoff_radius;

    auto const r2in = inner_radius_ * inner_radius_;
    auto const r2out = outer_radius_ * outer_radius_;

    // smoothing function - see Hernandez 2019 https://www.nature.com/articles/s41524-019-0249-1/
    auto f = [&](auto r)  {
        auto r2 = r * r;
        auto s = (2 * r2 - 3 * r2in + r2out) * (r2out - r2) * (r2out - r2) * std::pow(r2out - r2in, -3);
        return static_cast<Operon::Scalar>(s);
    };

    data_.clear();

    vector<Operon::Scalar> distances;
    distances.reserve(size_t{static_cast<size_t>(cluster_size * cluster_size)});

    for (auto const& c : coordinates) {
        distances.clear();

        for(auto i = 0; i < c.rows() - 1; ++i) { // NOLINT
            for (auto j = i+1; j < c.rows(); ++j) { // NOLINT
                auto d = (c.row(i) - c.row(j)).matrix().norm();
                if (d > cutoff_radius) { continue; }
                distances.push_back(d);
            }
        }

        // copy the distances into an Eigen::Array
        Array<Operon::Scalar, -1, 3> m(ssize(distances), 3);
        m.col(0) = Map<Array<Operon::Scalar, -1, 1>>(distances.data(), ssize(distances));
        m.col(1) = m.col(0).inverse();
        m.col(2) = m.col(0).unaryExpr(f);

        // construct a dataset from the array and add it to our data_
        data_.emplace_back(m);
        data_.back().SetVariableNames({ "r", "q", "s" });

        auto vars = data_.back().Variables();
        for (const auto& v : vars) {
            if (v.Name == "r") { ENSURE(v.Index == 0); }
            if (v.Name == "q") { ENSURE(v.Index == 1); }
            if (v.Name == "s") { ENSURE(v.Index == 2); }
        }
    }
    index_.resize(data_.size());
    std::iota(index_.begin(), index_.end(), 0UL);
}

} // namespace atomic
