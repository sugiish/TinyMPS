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
	Grid(int particles_number, int dimension, const MatrixXd& coordinates, const Vector3d& lower_coordinate, const Vector3d& higher_coordinate);
	virtual ~Grid();

	
private:
	int particles_number;
	int dimension;

	const MatrixXd& coordinates;
	const Vector3d& lower_coordinate;
	const Vector3d& higher_coordinate;

	std::vector<long long> hash;
	std::vector<int> index;

	double grid_size;
	int grid_number[3];

	void resetHash();
	int getHashValue(const Vector3d& position);

};


#endif //MPS_GRID_H_INCLUDED
