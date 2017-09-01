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
	particles.timer.initialize(condition.initial_time, condition.finish_time, condition.delta_time);

	
	while(particles.timer.hasNextLoop())
	{
		particles.moveParticlesExplicitly(condition.gravity);
		

		particles.timer.update();
		//std::cout << particles.timer.getCurrentTime() << std::endl;
	}

	particles.writeVtkFile("out.vtk", "test");
	
}