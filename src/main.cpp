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
	
	while(timer.hasNextLoop())
	{
		particles.moveParticlesExplicitly(condition.gravity, timer);
		

		timer.update();
		//std::cout << particles.timer.getCurrentTime() << std::endl;
	}

	particles.writeVtkFile("out.vtk", "test");
	
}