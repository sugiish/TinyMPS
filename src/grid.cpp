#include "grid.h"
#include "particles.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

Grid::Grid(const Particles& particles, double grid_width) : particles(particles), grid_width(grid_width)
{
	int dim = particles.getDimension();
	particles.getMaxPosition(higher_bounds);
	particles.getMinPosition(lower_bounds);

	VectorXd diff = higher_bounds - lower_bounds;
	for(int i = 0; i < 3; i++)
	{
		if(i < dim)grid_number[i] = std::ceil(diff(i) / grid_width);
		else grid_number[i] = 0;
	}

	resetHash();

}

Grid::~Grid()
{

}

void
Grid::resetHash()
{
	int pt_num = particles.getParticlesNumber();
	hash.resize(pt_num);
	
	for(int i = 0; i < pt_num; i++)
	{
		hash[i] = std::make_pair(getHashValue(particles.position.row(i)), i);
	}

	std::sort(hash.begin(), hash.end());
}

int
Grid::getHashValue(const VectorXd& position) const
{
	int grid_index[3];
	for(int i = 0; i < 3; i++)
	{
		if(i < particles.getDimension())grid_index[i] = std::ceil((position(i) - lower_bounds(i)) / grid_width);
		else grid_index[i] = 0;
	}

	return grid_index[0] + grid_index[1] * grid_number[0] + grid_index[2] * grid_number[1] * grid_number[0];
}