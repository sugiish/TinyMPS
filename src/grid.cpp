#include "grid.h"
#include "particles.h"

#include <algorithm>
#include <cmath>
#include <iostream>

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
	begin_hash.clear();
}

void
Grid::getNeighbor(int hash, int& begin, int& end)
{
	if(begin_hash.empty())
	{
		begin = -1; end = -1; return;
	}

	if(begin_hash.find(hash) == begin_hash.end())
	{
		begin = -1; end = -1; return;
	}

	begin = begin_hash[hash].first;
	end = begin_hash[hash].second;

}

void
Grid::resetHash()
{
	if(!begin_hash.empty()) begin_hash.clear();

	int pt_num = particles.getParticlesNumber();
	if(pt_num == 0) return;

	grid_hash.resize(pt_num);
	for(int i = 0; i < pt_num; i++)
	{
		grid_hash[i] = std::make_pair(getHashValue(particles.position.row(i)), i);
	}

	std::sort(grid_hash.begin(), grid_hash.end());

	int start_i = 0;
	int start_value = grid_hash[0].first;
	for(int i = 1; i < pt_num; i++)
	{
		if(start_value != grid_hash[i].first)
		{
			begin_hash[start_value] = std::make_pair(start_i, i - 1);
			start_i = i;
			start_value = grid_hash[i].first;
		}

		if(i == pt_num - 1)
		{
			begin_hash[start_value] = std::make_pair(start_i, i);
		}
	}
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