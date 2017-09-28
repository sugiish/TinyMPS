// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_GRID_H_INCLUDED
#define MPS_GRID_H_INCLUDED

#include <iostream>
#include <unordered_map>
#include <vector>
#include <Eigen/Dense>

namespace tiny_mps {

// Searches neighbor particles.
class Grid {
public:
    // Used for describing containers of neighbor particles.
    using Neighbors = std::vector<int>;

    Grid(double grid_width, const Eigen::MatrixXd& coordinates, const Eigen::Matrix<bool, Eigen::Dynamic, 1>& valid_coordinates, int dimension);
    virtual ~Grid();
    void getNeighbors(int index, Neighbors& neighbors) const;
    inline int getSize() const { return size; }
    inline int getDimension() const { return dimension; }
    inline double getGridWidth() const { return grid_width; }

private:
    void setHash();
    inline void getGridHashBegin(int hash, int& begin, int& end) const {
        if (begin_hash.empty()) {
            begin = -1; end = -1; return;
        }
        if (begin_hash.find(hash) == begin_hash.end()) {
            begin = -1; end = -1; return;
        }
        begin = begin_hash.at(hash).first;
        end = begin_hash.at(hash).second;
    }
    inline void getMaxCoordinates(Eigen::Vector3d& answer) const {
        answer = (coordinates.array().rowwise() * valid_coordinates.cast<double>().transpose().array()).rowwise().maxCoeff();
    }
    inline void getMinCoordinates(Eigen::Vector3d& answer) const {
        answer = (coordinates.array().rowwise() * valid_coordinates.cast<double>().transpose().array()).rowwise().minCoeff();
    }
    inline int getGridNumberX() const { return grid_number[0]; }
    inline int getGridNumberY() const { return grid_number[1]; }
    inline int getGridNumberZ() const {
        if(dimension == 2) return 0;
        return grid_number[2];
    }
    inline int toHash(const Eigen::Vector3d& vec) const {
        int dx, dy, dz;
        toIndex(vec, dx, dy, dz);
        return toHash(dx, dy, dz);
    }
    inline int toHash(int index_x, int index_y) const {
        return index_x + index_y * grid_number[0];
    }
    inline int toHash(int index_x, int index_y, int index_z) const {
        if (dimension == 2) return toHash(index_x, index_y);
        return index_x + index_y * grid_number[0] + index_z * grid_number[1] * grid_number[0];
    }
    inline void toIndex(const Eigen::Vector3d& vec, int& dx, int& dy, int& dz) const {
        dx = std::ceil((vec(0) - lower_bounds(0)) / grid_width);
        dy = std::ceil((vec(1) - lower_bounds(1)) / grid_width);
        if (dimension == 3) dz = std::ceil((vec(2) - lower_bounds(2)) / grid_width);
        else dz = 0;
    }
    inline void toIndex(int hash, int& dx, int& dy) const {
        dy = hash / grid_number[0];
        dx = hash % grid_number[0];
    }
    inline void toIndex(int hash, int& dx, int& dy, int& dz) const {
        if (dimension == 2) {
             toIndex(hash, dx, dy);
             dz = 0;
        } else {
            dz = hash / (grid_number[1] * grid_number[0]);
            dy = (hash % (grid_number[1] * grid_number[0])) / grid_number[0];
            dx = (hash % (grid_number[1] * grid_number[0])) % grid_number[0];
        }
    }

    /// hash, index
    std::vector<std::pair<int, int> > grid_hash;
    /// index, begin(order) end(order)
    std::unordered_map<int, std::pair<int, int> > begin_hash;
    const Eigen::MatrixXd coordinates;
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> valid_coordinates;
    int dimension;
    int size;
    int initial_neighbors_size;
    Eigen::Vector3d higher_bounds;
    Eigen::Vector3d lower_bounds;
    double grid_width;
    int grid_number[3];
};

} // namespace tiny_mps
#endif //MPS_GRID_H_INCLUDED