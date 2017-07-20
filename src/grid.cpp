#include "grid.h"

#include <cmath>
#include <vector>

Grid::Grid(int particles_number, int dimension, const MatrixXd & coordinates, const Vector3d & lower_coordinate, const Vector3d & higher_coordinate)
	: particles_number(particles_number), dimension(dimension), coordinates(coordinates), lower_coordinate(lower_coordinate), higher_coordinate(higher_coordinate)
{
	//hash = std::vector(particles_number);
	//index = std::vector(particles_number);

	Vector3d diff = higher_coordinate - lower_coordinate;
	for(int i = 0; i < 3; i++)
	{
		if(i < dimension)grid_number[i] = std::ceil(diff(i) / grid_size);
		else grid_number[i] = 0;
	}
}

Grid::~Grid()
{

}

void
Grid::resetHash()
{


}

int
getHashValue(const Vector3d& position)
{
	int grid_index[3];
	for(int i = 0; i < 3; i++)
	{
		if(i < dimension)grid_index[i] = std::ceil(position(i) / grid_size);
		else grid_index[i] = 0;
	}

	return grid_index[0] + grid_index[1] * grid_number[0] + grid_index[2] * grid_number[1] * grid_number[0];

}