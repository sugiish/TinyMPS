#ifndef MPS_GRID_H_INCLUDED
#define MPS_GRID_H_INCLUDED

#include <unordered_map>
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

	void getNeighbor(int hash, int& begin, int& end);
	
	inline void getGridHash(std::vector<std::pair<int, int> >& ghash)
	{
		ghash = grid_hash;
	}

	void resetHash();
	int getHashValue(const VectorXd& position) const;

	inline int indexToHashValue(int index_x, int index_y)
	{
		return index_x + index_y * grid_number[0];
	}

	inline int indexToHashValue(int index_x, int index_y, int index_z)
	{
		return index_x + index_y * grid_number[0] + index_z * grid_number[1] * grid_number[0];
	}


private:
	const Particles& particles;

	VectorXd lower_bounds;
	VectorXd higher_bounds;

	double grid_width;
	int grid_number[3];

	// hash, index
	std::vector<std::pair<int, int> > grid_hash;
	
	// index, begin(order) end(order)
	std::unordered_map<int, std::pair<int, int> > begin_hash;
};


#endif //MPS_GRID_H_INCLUDED
