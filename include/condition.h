#ifndef MPS_CONDITION_H_INCLUDED
#define MPS_CONDITION_H_INCLUDED

#include "reader.h"

#include <Eigen/Dense>

using namespace Eigen;

class Condition
{
public:
	Condition(std::string path) : reader(path)
	{
		reader.getValue("average_distance",  average_distance);
		reader.getValue("dimension", dimension);
		if(dimension != 2 && dimension != 3)
		{
			std::cerr << "Error: " << dimension << "-dimension is not supported." << std::endl;
		}

		initial_particle_number_density = 0;

		double gx, gy, gz;
		reader.getValue("gravity_x", gx);
		reader.getValue("gravity_y", gy);
		reader.getValue("gravity_z", gz);
		gravity(0) = gx;
		gravity(1) = gy;
		if(dimension == 3) gravity(2) = gz;
		else gravity(2) = gz;

		reader.getValue("temperature", temperature);
		reader.getValue("head_pressure", head_pressure);

		reader.getValue("viscosity_calculation", viscosity_calculation);
		reader.getValue("kinematic_viscosity", kinematic_viscosity);

		reader.getValue("courant_number", courant_number);
		reader.getValue("diffusion_number", diffusion_number);

		reader.getValue("pnd_influence", pnd_influence);
		reader.getValue("gradient_influence", gradient_influence);
		reader.getValue("laplacian_influence", laplacian_influence);
	}

	virtual ~Condition(){}

	Reader reader;

	double average_distance;
	int dimension;
	double initial_particle_number_density;

	Vector3d gravity;
	double temperature;
	double head_pressure;

	bool viscosity_calculation;
	double kinematic_viscosity;

	double courant_number;
	double diffusion_number;

	double pnd_influence;
	double gradient_influence;
	double laplacian_influence;
private:

};


#endif //MPS_CONDITION_H_INCLUDED