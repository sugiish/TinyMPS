#ifndef MPS_CONDITION_H_INCLUDED
#define MPS_CONDITION_H_INCLUDED

#include "reader.h"

#include <Eigen/Dense>

using namespace Eigen;

class Condition
{
public:
	Condition(Reader& reader){
		reader.getValue("average_distance",  average_distance);
		reader.getValue("dimension", dimension);
		if(dimension != 2 && dimension != 3)
		{
			std::cerr << "Error: " << dimension << "-dimension is not supported." << std::endl;
		}

		gravity = VectorXd::Zero(dimension);
		initial_particle_number_density = 0;

		double gx, gy, gz;
		reader.getValue("gravity_x", gx);
		reader.getValue("gravity_y", gy);
		reader.getValue("gravity_z", gz);
		gravity(0) = gx;
		gravity(1) = gy;
		if(dimension == 3) gravity(2) = gz;

		reader.getValue("temperature", temperature);
		reader.getValue("head_pressure", head_pressure);

		reader.getValue("viscosity_calculation", viscosity_calculation);
		reader.getValue("kinematic_viscosity", kinematic_viscosity);

		reader.getValue("courant_number", courant_number);
		reader.getValue("diffusion_number", diffusion_number);

		reader.getValue("gradient_influence", gradient_influence);
		reader.getValue("laplacian_influence", laplacian_influence);
		
	}

	virtual ~Condition(){}

	double average_distance;
	int dimension;
	double initial_particle_number_density;

	VectorXd gravity;
	double temperature;
	double head_pressure;

	bool viscosity_calculation;
	double kinematic_viscosity;

	double courant_number;
	double diffusion_number;

	double gradient_influence;
	double laplacian_influence;

};


#endif //MPS_CONDITION_H_INCLUDED