#include "particles.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <Eigen/Dense>

#include "timer.h"


/**
 * Class for describing particles status.
 */
Particles::Particles(int particles_number, int dimension)
{
	initialize(particles_number, dimension);

}

Particles::Particles(const string& path, int dimension)
{
	readGridFile(path, dimension);
}

Particles::~Particles()
{
	
}

void
Particles::initialize(int particles_number, int dimension)
{
	this->particles_number = particles_number;
	this->dimension = dimension;
	
	position = MatrixXd::Zero(particles_number, dimension);
	velocity = MatrixXd::Zero(particles_number, dimension);
	pressure = VectorXd::Zero(particles_number);

	temporal_position = MatrixXd::Zero(particles_number, dimension);
	temporal_velocity = MatrixXd::Zero(particles_number, dimension);

	particles_type = VectorXi::Zero(particles_number);
	ghost_particles = VectorXi::Zero(particles_number);
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
		initialize(ptcl_num, dimension);
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
			if(i_dim < dimension)position(i_counter, i_dim) = tmp;
			
		}

		for(int i_dim = 0; i_dim < 3; i_dim++)
		{
			double tmp;
			ss >> tmp;
			if(i_dim < dimension)velocity(i_counter, i_dim) = tmp;
		}

		ss >> pressure(i_counter);

		i_counter++;
    }

	return 0;
}