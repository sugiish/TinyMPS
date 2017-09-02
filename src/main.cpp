#include <iostream>

#include <Eigen/Dense>

#include "condition.h"
#include "grid.h"
#include "particles.h"

/**
 * This is a sample code using mps library.
 */
int main() {
	tiny_mps::Condition condition("./input/input.data");
	tiny_mps::Particles particles("./input/dambreak.grid", condition);
	tiny_mps::Timer timer(condition);
	Eigen::Matrix<bool, Eigen::Dynamic, 1> valid = particles.particle_types.array() != tiny_mps::ParticleType::GHOST;
	double re = condition.pnd_influence * condition.average_distance;
	tiny_mps::Grid grid(re, particles.position, valid, condition.dimension);
	auto weight = [&particles, re](int i,int j){
		Eigen::Vector3d v = particles.position.col(j) - particles.position.col(i);
		double r = v.norm();
		if(r < re) return (re/r - 1.0);
		else return 0.0;
	};
	/*auto number = [&particles, re](int i,int j){
		Vector3d v = particles.position.col(j) - particles.position.col(i);
		double r = v.norm();
		if(r < re) return 1;
		else return 0;
	};*/
	particles.updateParticleNumberDensity(grid, weight);
	particles.setInitialParticleNumberDensity(condition.initial_pnd_index);

	while(timer.hasNextLoop())
	{
		particles.moveParticlesExplicitly(condition.gravity, timer);
		

		timer.update();
	}

	particles.writeVtkFile("out.vtk", "test");
}