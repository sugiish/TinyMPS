#ifndef MPS_CONDITION_H_INCLUDED
#define MPS_CONDITION_H_INCLUDED

#include <Eigen/Dense>

class Condition
{
public:
	Vector3d gravity;
	double temperature;
	double head_pressure;

	double kinematic_viscosity;


};


#endif //MPS_CONDITION_H_INCLUDED