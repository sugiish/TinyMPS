#include "particles.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <Eigen/Dense>


/**
 * Class for describing particles status.
 */
Particles::Particles(int particles_number, int dimension)
{
	initialize(particles_number, dimension);

}

Particles::Particles(string path, int dimension)
{
	
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
Particles::readGridFile(string path, int dimension)
{
	int ptcl_num;
	string tmp_str;

	ifstream ifs(path);
	stringstream ss;

    if(ifs.fail())
    {
        std::cerr << "Error: in readGridFile()" << std::endl;
        std::cerr << "Failed to read files: " << path << std::endl;

		return 1;
    }
	//Line 0: Start time
	getline(ifs, tmp_str);

	//Line 1: particles_number
	getline(ifs, tmp_str);
	ss.clear();
	ss.str(tmp_str);
	ss >> ptcl_num;

	initialize(ptcl_num, dimension);

    while(getline(ifs, tmp_str))
    {
        std::cout << "[" << tmp_str << "]" << std::endl;
    }


	return 0;
}