#ifndef MPS_GRID_H_INCLUDED
#define MPS_GRID_H_INCLUDED

#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace tiny_mps {

/**
 * Class for neighbor searching. 
 */
class Grid {
public:
	Grid(double grid_width, const Eigen::MatrixXd& coordinates, const Eigen::Matrix<bool, Eigen::Dynamic, 1>& valid_coordinates, int dimension);
	virtual ~Grid();

	void sumNeighborScalars(Eigen::VectorXd& output, std::function<double(int, int)> interaction);
	void sumNeighborVectors(Eigen::MatrixXd& output, std::function<void(int, int, const Eigen::Vector3d&)> interaction);
	void resetHash();	

	inline void getGridHash(std::vector<std::pair<int, int> >& ghash) const {
		ghash = grid_hash;
	}	
	inline int getSize() const { return size; }
	inline int getDimension() const { return dimension; }
	inline int getGridNumberX() const { return grid_number[0]; }
	inline int getGridNumberY() const { return grid_number[1]; }
	inline int getGridNumberZ() const {
		if(dimension == 2) return 0;
		return grid_number[2];
	}

	inline int toHash(const Eigen::Vector3d& coordinates) const {
		int dx, dy, dz;
		toIndex(coordinates, dx, dy, dz);
		return toHash(dx, dy, dz);
	}
	inline int toHash(int index_x, int index_y) const {
		return index_x + index_y * grid_number[0];
	}
	inline int toHash(int index_x, int index_y, int index_z) const {
		if (dimension == 2) return toHash(index_x, index_y);
		return index_x + index_y * grid_number[0] + index_z * grid_number[1] * grid_number[0];
	}

	inline void toIndex(const Eigen::Vector3d& coordinates, int& dx, int& dy, int& dz) const {
		dx = std::ceil((coordinates(0) - lower_bounds(0)) / grid_width);
		dy = std::ceil((coordinates(1) - lower_bounds(1)) / grid_width);
		if (dimension == 3) dz = std::ceil((coordinates(2) - lower_bounds(2)) / grid_width);
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

private:
	void getNeighbors(int hash, int& begin, int& end);
	inline void getMaxCoordinates(Eigen::Vector3d& answer) const {
		answer = coordinates.rowwise().maxCoeff();
	}
	inline void getMinCoordinates(Eigen::Vector3d& answer) const {
		answer = coordinates.rowwise().minCoeff();
	}

	const Eigen::MatrixXd& coordinates;
	const Eigen::Matrix<bool, Eigen::Dynamic, 1>& valid_coordinates;
	int dimension;
	int size;

	Eigen::Vector3d lower_bounds;
	Eigen::Vector3d higher_bounds;

	double grid_width;
	int grid_number[3];

	// hash, index
	std::vector<std::pair<int, int> > grid_hash;
	// index, begin(order) end(order)
	std::unordered_map<int, std::pair<int, int> > begin_hash;
};

} // namespace tiny_mps
#endif //MPS_GRID_H_INCLUDED