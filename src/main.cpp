#include <iostream>

#include <Eigen/Dense>

#include "condition.h"
#include "grid.h"
#include "particles.h"

using namespace Eigen;
using namespace std;

/**
 * This is a sample code using mps library.
 */
int main()
{
	Condition condition("./input/input.data");
	Particles particles("./input/dambreak.grid", condition);
	Timer timer(condition.initial_time, condition.finish_time, condition.delta_time);
	
	double re = condition.pnd_influence * condition.average_distance;
	Grid grid(re, particles.position, particles.particles_valid, condition.dimension);
	auto weight = [&particles, re](int i,int j){
		Vector3d v = particles.position.col(j) - particles.position.col(i);
		double r = v.norm();
		if(r < re) return (re/r - 1.0);
		else return 0.0;
	};
	auto number = [&particles, re](int i,int j){
		Vector3d v = particles.position.col(j) - particles.position.col(i);
		double r = v.norm();
		if(r < re) return 1;
		else return 0;
	};
	particles.updateParticleNumberDensity(grid, weight);

	while(timer.hasNextLoop())
	{
		particles.moveParticlesExplicitly(condition.gravity, timer);
		

		timer.update();
	}

	particles.writeVtkFile("out.vtk", "test");
	
}