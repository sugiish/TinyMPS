#include "grid.h"

Grid::Grid(int dimension, const MatrixXd & coordinates, const Vector3d & lower_coordinate, const Vector3d & higher_coordinate)
	: dimension(dimension), coordinates(coordinates), lower_coordinate(lower_coordinate), higher_coordinate(higher_coordinate)
{

}

Grid::~Grid()
{

}