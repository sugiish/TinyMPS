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

};

	
#endif //MPS_DIFFERENTIAL_H_INCLUDED