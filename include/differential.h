#ifndef MPS_DIFFERENTIAL_H_INCLUDED
#define MPS_DIFFERENTIAL_H_INCLUDED

#include <Eigen/Dense>

class Differential
{
public:
	inline double weightFunction(double r, double r_e)
	{
		return (0 <= r && r < r_e) ? 
			r / r_e - 1
			: 0;
	}

	inline void gradient(VectorXd& scalars, VectorXd& output, const Grid& grid, const Condition& condition)
	{
		VectorXd ans = VectorXd::Zero(grid.getDimension());
		for(int j = 0; j < grid.size; j++)
		{
			if(grid.isValidCoordinate(j) == 0) continue;
			std::vector neighbors;
			grid.getNeighbors(j, neighbors);

			for(int i : neighbors)
			{
				VectorXd r_ji= coordinates.row(j) - coordinates.row(i);
				ans = ans + (scalars(i) - scalars(j)) * r_ji * weightFunction(r_ji.distance(), condition.gradient_influence) 
				/ r_ji.squaredNorm(); ;
			}

		}
	}

};

	
#endif //MPS_DIFFERENTIAL_H_INCLUDED