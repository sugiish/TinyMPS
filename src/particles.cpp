#include "particles.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <Eigen/Dense>

#include "timer.h"


/**
 * Class for describing particles status.
 */
Particles::Particles(const string& path, Condition& condition): condition(condition)
{
	readGridFile(path, condition.dimension);
}

Particles::~Particles()
{
	
}

void
Particles::initialize(int particles_number)
{
	this->particles_number = particles_number;
	
	position = MatrixXd::Zero(3, particles_number);
	velocity = MatrixXd::Zero(3, particles_number);
	pressure = VectorXd::Zero(particles_number);

	temporary_position = MatrixXd::Zero(3, particles_number);
	temporary_velocity = MatrixXd::Zero(3, particles_number);

	particles_type = VectorXi::Zero(particles_number);
	particles_valid = VectorXi::Zero(particles_number);
}

int
Particles::readGridFile(const string& path, int dimension)
{
	ifstream ifs(path);
	
    if(ifs.fail())
    {
        std::cerr << "Error: in readGridFile()" << std::endl;
        std::cerr << "Failed to read files: " << path << std::endl;

		return 1;
    }

	string tmp_str;

	//Line 0: Start time
	getline(ifs, tmp_str);
	

	//Line 1: particles_number
	getline(ifs, tmp_str);

	{
		stringstream ss;
		int ptcl_num = 0;

		ss.str(tmp_str);
		ss >> ptcl_num;
		initialize(ptcl_num);
	}

	int i_counter = 0;
    while(getline(ifs, tmp_str))
    {
		stringstream ss;
		ss.str(tmp_str);
        
		ss >> particles_type(i_counter);

		for(int i_dim = 0; i_dim < 3; i_dim++)
		{
			double tmp;
			ss >> tmp;
			if(i_dim < dimension)position(i_dim, i_counter) = tmp;
			
		}

		for(int i_dim = 0; i_dim < 3; i_dim++)
		{
			double tmp;
			ss >> tmp;
			if(i_dim < dimension)velocity(i_dim, i_counter) = tmp;
		}

		ss >> pressure(i_counter);

		i_counter++;
    }

	return 0;
}



void
Particles::moveParticlesExplicitly(double delta_time, const VectorXd& force)
{
	temporary_velocity = velocity + delta_time * (force);
	temporary_position = position + delta_time * temporary_velocity;
}
