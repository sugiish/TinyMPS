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
	Grid(double grid_width, const MatrixXd& coordinates, const VectorXi& valid_coordinates, int dimension);
	virtual ~Grid();

	void getNeighbor(int hash, int& begin, int& end);

	bool getNeighbors(int index, std::vector<int>& neighbors);
	
	inline void getGridHash(std::vector<std::pair<int, int> >& ghash) const
	{
		ghash = grid_hash;
	}

	void resetHash();
	int getHashValue(const VectorXd& position) const;

	inline int getSize() const
	{
		return size;
	}

	inline int getDimension() const
	{
		return dimension;
	}

	inline int getGridNumberX() const
	{
		return grid_number[0];
	}

	inline int getGridNumberY() const
	{
		return grid_number[1];
	}

	inline int getGridNumberZ() const
	{
		if(dimension == 2) return 0;
		return grid_number[2];
	}
	
	inline void getMaxPosition(VectorXd& answer) const
	{
		answer = coordinates.rowwise().maxCoeff();
	}

	inline void getMinPosition(VectorXd& answer) const
	{
		answer = coordinates.rowwise().minCoeff();
	}	

	inline int isValidCoordinate(int index)
	{
		return valid_coordinates(index);
	}

	inline int indexToHashValue(int index_x, int index_y) const
	{
		return index_x + index_y * grid_number[0];
	}

	inline int indexToHashValue(int index_x, int index_y, int index_z) const
	{
		return index_x + index_y * grid_number[0] + index_z * grid_number[1] * grid_number[0];
	}

	inline void hashValueToIndex(int hash, int& dx, int& dy)
	{
		dy = hash / grid_number[0];
		dx = hash % grid_number[0];
	}

	inline void hashValueToIndex(int hash, int& dx, int& dy, int& dz)
	{
		if(dimension == 2)
		{
			 hashValueToIndex(hash, dx, dy);
			 dz = 0;
		}
		else
		{
			dz = hash / (grid_number[1] * grid_number[0]);
			dy = (hash % (grid_number[1] * grid_number[0])) / grid_number[0];
			dx = (hash % (grid_number[1] * grid_number[0])) % grid_number[0];
		}
	}


private:
	const MatrixXd& coordinates;
	const VectorXi& valid_coordinates;
	int dimension;
	int size;

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
