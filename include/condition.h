#ifndef MPS_CONDITION_H_INCLUDED
#define MPS_CONDITION_H_INCLUDED

#include <Eigen/Dense>

class Condition
{
public:
	double average_distnace;

	Vector3d gravity;
	double temperature;
	double head_pressure;
	double kinematic_viscosity;

	double cfl_number;
	double diffusion_number;

	double gradient_ratio;
	double laplacian_ratio;

};


#endif //MPS_CONDITION_H_INCLUDED