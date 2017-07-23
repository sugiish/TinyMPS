#ifndef MPS_GRID_H_INCLUDED
#define MPS_GRID_H_INCLUDED

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

class Particles;

/**
 * Class for neighbor searching. 
 */
class Grid
{
public:
	Grid(const Particles& particles, double grid_width);
	virtual ~Grid();

private:
	const Particles& particles;

	VectorXd lower_bounds;
	VectorXd higher_bounds;

	std::vector<std::pair<int, int> > hash;
	std::vector<int> index;

	double grid_width;
	int grid_number[3];

	void resetHash();
	int getHashValue(const VectorXd& position) const;

};


#endif //MPS_GRID_H_INCLUDED
