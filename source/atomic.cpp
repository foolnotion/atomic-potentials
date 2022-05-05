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

// this function loads the snaphots of atomic positions into separate datasets (one per snapshot)
// only distances below the outer cutoff radius are added
auto summation_function::load_data(std::string const& path, int cluster_size, double cutoff_radius) -> void {
    auto coordinates = parse_coordinates(path, cluster_size);
    outer_radius_ = cutoff_radius;

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
        distances.reserve(size_t{static_cast<size_t>(cluster_size * cluster_size)});

        std::vector<Operon::Scalar> tmp;
        for (auto i = 0; i < dist.cols(); ++i) {
            auto col = dist.col(i);
            for (int j = 0; j < col.size(); ++j) {
                if (i != j && col(j) < cutoff_radius) {
                    tmp.push_back(col(j));
                }
            }
            if (!tmp.empty()) {
                //std::stable_sort(tmp.begin(), tmp.end());
                std::copy(tmp.begin(), tmp.end(), std::back_inserter(distances));
            }
            tmp.clear();
        }

        // copy the distances into an Eigen::Array
        Eigen::Array<Operon::Scalar, -1, 1> m = Eigen::Map<decltype(m)>(distances.data(), ssize(distances));

        // construct a dataset from the array and add it to our data_
        data_.emplace_back(m);
        data_.back().SetVariableNames({ "r" });
    }
    index_.resize(data_.size());
    std::iota(index_.begin(), index_.end(), 0UL);
    fmt::print("data size: {}, dist size: {}\n", data_.size(), data_.front().Rows());
}

} // namespace atomic
