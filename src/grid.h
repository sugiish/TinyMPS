#ifndef MPS_GRID_H_INCLUDED
#define MPS_GRID_H_INCLUDED

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

/**
 * Class for neighbor searching. 
 */
class Grid
{
public:
	Grid(int dimension, const MatrixXd& coordinates, const Vector3d& lower_coordinate, const Vector3d& higher_coordinate);
	virtual ~Grid();

private:
	int dimension;

	const MatrixXd& coordinates;
	const Vector3d& lower_coordinate;
	const Vector3d& higher_coordinate;

	std::vector<int> hash;
	std::vector<int> index;

	double grid_size;
	int grid_number[3];

};


#endif //MPS_GRID_H_INCLUDED
