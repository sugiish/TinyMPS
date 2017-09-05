#include "grid.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>

namespace tiny_mps {

Grid::Grid(double grid_width, const Eigen::MatrixXd& coordinates, const Eigen::Matrix<bool, Eigen::Dynamic, 1>& valid_coordinates, int dimension) 
	: grid_hash(coordinates.cols()), coordinates(coordinates), valid_coordinates(valid_coordinates) {
	this->grid_width = grid_width;
	this->dimension = dimension;
	size = coordinates.cols();
	initial_neighbors_size = (int)(std::pow(grid_width, dimension) * 2);
	resetHash();
}

Grid::~Grid() {
	begin_hash.clear();
}

void Grid::getNeighbors(int index, std::vector<int>& neighbors) {
	neighbors.clear();
	if(valid_coordinates(index) == 0) return;
	int ix, iy, iz;
	toIndex(coordinates.col(index), ix, iy, iz);
	int x_begin, x_end, y_begin, y_end, z_begin, z_end;
	x_begin = std::max(ix - 1, 0);
	y_begin = std::max(iy - 1, 0);
	z_begin = std::max(iz - 1, 0);
	x_end = std::min(ix + 1, getGridNumberX() - 1);
	y_end = std::min(iy + 1, getGridNumberY() - 1);
	z_end = std::min(iz + 1, getGridNumberZ() - 1);
	if(dimension == 2) {
		z_begin = 0;
		z_end = 0;
	}
	
	for (int i = z_begin; i <= z_end; i++) {
		for (int j = y_begin; j <= y_end; j++) {
			for (int k = x_begin; k <= x_end; k++) {
				int begin, end;
				getGridHashBegin(toHash(k, j, i), begin, end);
				if (begin == -1 || end == -1) continue;

				Eigen::Vector3d r_i = coordinates.col(index);
				for (int n = begin; n <= end; n++) {
					int j_particle = grid_hash[n].second;
					if(index == j_particle) continue;
					Eigen::Vector3d r_j = coordinates.col(j_particle);
					Eigen::Vector3d r_ji = r_j - r_i;
					if(r_ji.norm() > grid_width) continue;
					neighbors.push_back(j_particle);
				}
			}
		}
	}
}

void Grid::getGridHashBegin(int hash, int& begin, int& end) {
	if (begin_hash.empty()) {
		begin = -1; end = -1; return;
	}
	if (begin_hash.find(hash) == begin_hash.end()) {
		begin = -1; end = -1; return;
	}
	begin = begin_hash[hash].first;
	end = begin_hash[hash].second;
}

double Grid::sumNeighborScalars(int index, std::function<double(int, int)> interaction) {
	if(valid_coordinates(index) == 0) return 0;	
	int ix, iy, iz;
	toIndex(coordinates.col(index), ix, iy, iz);
	int x_begin, x_end, y_begin, y_end, z_begin, z_end;
	x_begin = std::max(ix - 1, 0);
	y_begin = std::max(iy - 1, 0);
	z_begin = std::max(iz - 1, 0);
	x_end = std::min(ix + 1, getGridNumberX() - 1);
	y_end = std::min(iy + 1, getGridNumberY() - 1);
	z_end = std::min(iz + 1, getGridNumberZ() - 1);
	if(dimension == 2) {
		z_begin = 0;
		z_end = 0;
	}
	
	double ans = 0;
	for (int i = z_begin; i <= z_end; i++) {
		for (int j = y_begin; j <= y_end; j++) {
			for (int k = x_begin; k <= x_end; k++) {
				int begin, end;
				getGridHashBegin(toHash(k, j, i), begin, end);
				if (begin == -1 || end == -1) continue;

				Eigen::Vector3d r_i = coordinates.col(index);
				for (int n = begin; n <= end; n++) {
					int j_particle = grid_hash[n].second;
					if (index == j_particle) continue;
					if (valid_coordinates(j_particle) == 0) continue;
					Eigen::Vector3d r_j = coordinates.col(j_particle);
					Eigen::Vector3d r_ji = r_j - r_i;
					if (r_ji.norm() > grid_width) continue;

					ans += interaction(index, j_particle);
				}
			}
		}
	}
	return ans;
}
void Grid::sumNeighborVectors(int index, std::function<void(int, int, Eigen::Vector3d&)> interaction, Eigen::Vector3d& output) {
	output << 0.0, 0.0, 0.0;
	if (valid_coordinates(index) == 0) return;
	int ix, iy, iz;
	toIndex(coordinates.col(index), ix, iy, iz);
	int x_begin, x_end, y_begin, y_end, z_begin, z_end;
	x_begin = std::max(ix - 1, 0);
	y_begin = std::max(iy - 1, 0);
	z_begin = std::max(iz - 1, 0);
	x_end = std::min(ix + 1, getGridNumberX() - 1);
	y_end = std::min(iy + 1, getGridNumberY() - 1);
	z_end = std::min(iz + 1, getGridNumberZ() - 1);
	if (dimension == 2) {
		z_begin = 0;
		z_end = 0;
	}

	for (int i = z_begin; i <= z_end; i++) {
		for (int j = y_begin; j <= y_end; j++) {
			for (int k = x_begin; k <= x_end; k++) {
				int begin, end;
				getGridHashBegin(toHash(k, j, i), begin, end);
				if (begin == -1 || end == -1) continue;

				Eigen::Vector3d r_i = coordinates.col(index);
				for (int n = begin; n <= end; n++) {
					int j_particle = grid_hash[n].second;
					if (index == j_particle) continue;
					Eigen::Vector3d r_j = coordinates.col(j_particle);
					Eigen::Vector3d r_ji = r_j - r_i;
					if (r_ji.norm() > grid_width) continue;

					Eigen::Vector3d tmp;
					interaction(index, j_particle, tmp);
					output += tmp;
				}
			}
		}
	}
}

void Grid::sumAllNeighborScalars(std::function<double(int, int)> interaction, Eigen::VectorXd& output) {
	output = Eigen::VectorXd::Zero(size);
	std::vector<int> v(initial_neighbors_size);
	for (int i_particle = 0; i_particle < size; i_particle++) {
		// output(i_particle) = sumNeighborScalars(i_particle, interaction);
		getNeighbors(i_particle, v);
		for(int j_particle : v) {
			output(i_particle) += interaction(i_particle, j_particle);
		}
	}
}

void Grid::sumAllNeighborVectors(std::function<void(int, int, Eigen::Vector3d&)> interaction, Eigen::Matrix3Xd& output) {
	output = Eigen::MatrixXd::Zero(3, size);
	for (int i_particle = 0; i_particle < size; i_particle++) {
		Eigen::Vector3d ans;
		sumNeighborVectors(i_particle, interaction, ans);
		output.col(i_particle) = ans;
	}
}

void Grid::resetHash() {
	if (!begin_hash.empty()) begin_hash.clear();
	if (size != coordinates.cols()) {
		size = coordinates.cols();
		grid_hash.resize(size);
		if (size == 0) return;
	}

	getMaxCoordinates(higher_bounds);
	getMinCoordinates(lower_bounds);
	Eigen::Vector3d diff = higher_bounds - lower_bounds;
	if (getDimension() == 2) diff(2) = 0;
	for (int i = 0; i < 3; i++) {
		grid_number[i] = std::ceil(diff(i) / grid_width) + 1;
	}
	for (int i = 0; i < size; i++) {
		grid_hash[i] = std::make_pair(toHash(coordinates.col(i)), i);
	}
	std::sort(grid_hash.begin(), grid_hash.end());

	int start_i = 0;
	int start_value = grid_hash[0].first;
	for (int i = 1; i < size; i++) {
		if (start_value != grid_hash[i].first) {
			begin_hash[start_value] = std::make_pair(start_i, i - 1);
			start_i = i;
			start_value = grid_hash[i].first;
		}

		if (i == size - 1) {
			begin_hash[start_value] = std::make_pair(start_i, i);
		}
	}
}

} // namespace tiny_mps