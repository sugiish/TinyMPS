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
	while(timer.hasNextLoop())
	{
		particles.moveParticlesExplicitly(condition.gravity, timer);
		timer.update();
	}
	particles.writeVtkFile("out.vtk", "test");
}