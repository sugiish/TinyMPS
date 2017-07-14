#include "particles.h"

#include <Eigen/Dense>


/**
 * Class for describing particles status.
 */
Particles::Particles(int particles_number, int dimension)
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

Particles::~Particles()
{
	
}