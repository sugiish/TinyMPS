#include "grid.h"
#include "particles.h"

#include <algorithm>
#include <cmath>
#include <iostream>

Grid::Grid(double grid_width, const MatrixXd& coordinates, const VectorXi& valid_coordinates, int dimension) : coordinates(coordinates), valid_coordinates(valid_coordinates)
{
	this->grid_width = grid_width;
	this->dimension = dimension;
	size = coordinates.size();

	getMaxPosition(higher_bounds);
	getMinPosition(lower_bounds);

	VectorXd diff = higher_bounds - lower_bounds;
	for(int i = 0; i < 3; i++)
	{
		if(i < getDimension())grid_number[i] = std::ceil(diff(i) / grid_width);
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

bool
Grid::getNeighbors(int index, std::vector<int>& neighbors)
{
	if(valid_coordinates(index) == 0) return false;

	if(dimension == 2)
	{
		int hash = getHashValue(coordinates.row(index));
		int ix, iy;
		hashValueToIndex(hash, ix, iy);
		int x_begin, x_end, y_begin, y_end;
		x_begin = std::max(ix - 1, 0);
		y_begin = std::max(iy - 1, 0);
		x_end = std::min(ix + 1, getGridNumberX() - 1);
		y_end = std::min(iy + 1, getGridNumberY() - 1);

		int total = 0;
		for(int i = y_begin; i <= y_end; i++)
		{
			for(int j = x_begin; j <= x_end; j++)
			{
				int begin, end;
				getNeighbor(indexToHashValue(j, i), begin, end);
				if(begin == -1 || end == -1) continue;
				total += end - begin;
			}
		}

		std::vector<int> v(total);
		for(int i = y_begin; i <= y_end; i++)
		{
			for(int j = x_begin; j <= x_end; j++)
			{
				int begin, end;
				getNeighbor(indexToHashValue(j, i), begin, end);
				if(begin == -1 || end == -1) continue;
				for(int n = begin; n <= end; n++)
				{
					int tmp = grid_hash[n].second;
					v.push_back(tmp);
				}
			}
		}

		neighbors = v;
	}
	else
	{
		int hash = getHashValue(coordinates.row(index));
		int ix, iy, iz;
		hashValueToIndex(hash, ix, iy, iz);
		int x_begin, x_end, y_begin, y_end, z_begin, z_end;
		x_begin = std::max(ix - 1, 0);
		y_begin = std::max(iy - 1, 0);
		z_begin = std::max(iz - 1, 0);
		x_end = std::min(ix + 1, getGridNumberX() - 1);
		y_end = std::min(iy + 1, getGridNumberY() - 1);
		z_end = std::min(iz + 1, getGridNumberZ() - 1);

		int total = 0;
		for(int i = z_begin; i <= z_end; i++)
		{
			for(int j = y_begin; j <= y_end; j++)
			{
				for(int k = x_begin; k <= x_end; k++)
				{
					int begin, end;
					getNeighbor(indexToHashValue(k, j, i), begin, end);
					if(begin == -1 || end == -1) continue;
					total += end - begin;
				}
			}
		}
		
		std::vector<int> v(total);
		for(int i = z_begin; i <= z_end; i++)
		{
			for(int j = y_begin; j <= y_end; j++)
			{
				for(int k = x_begin; k <= x_end; k++)
				{
					int begin, end;
					getNeighbor(indexToHashValue(k, j, i), begin, end);
					if(begin == -1 || end == -1) continue;
					for(int n = begin; n <= end; n++)
					{
						int tmp = grid_hash[n].second;
						v.push_back(tmp);
					}
				}
			}
		}

		neighbors = v;
	}
	return true;
}

void
Grid::resetHash()
{
	if(!begin_hash.empty()) begin_hash.clear();

	int pt_num = coordinates.size();
	if(pt_num == 0) return;

	grid_hash.resize(pt_num);
	for(int i = 0; i < pt_num; i++)
	{
		grid_hash[i] = std::make_pair(getHashValue(coordinates.row(i)), i);
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
		if(i < getDimension())grid_index[i] = std::ceil((position(i) - lower_bounds(i)) / grid_width);
		else grid_index[i] = 0;
	}

	return grid_index[0] + grid_index[1] * grid_number[0] + grid_index[2] * grid_number[1] * grid_number[0];
}