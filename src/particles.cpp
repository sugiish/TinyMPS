#include "particles.h"

#include <Eigen/Dense>

Particles::Particles(int particles_number, int dimension)
{
	this->particles_number = particles_number;
	this->dimension = dimension;
	
	position = MatrixXd::Zero(particles_number, dimension);
	velocity = MatrixXd::Zero(particles_number, dimension);
	pressure = VectorXd::Zero(particles_number);

	temporal_position = MatrixXd::Zero(particles_number, dimension);
	temporal_velocity = MatrixXd::Zero(particles_number, dimension);

	type = VectorXi::Zero(particles_number);



}

Particles::~Particles()
{
	
}