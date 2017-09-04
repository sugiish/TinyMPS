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
	Eigen::Matrix<bool, 1, Eigen::Dynamic> lap_valid = particles.particle_types.array() == tiny_mps::ParticleType::NORMAL;
	Eigen::Matrix<bool, 1, Eigen::Dynamic> pnd_valid = particles.particle_types.array() == tiny_mps::ParticleType::NORMAL;	
	tiny_mps::Grid lap_grid(particles.laplacian_pressure_weight_radius, particles.position, lap_valid, condition.dimension);
	tiny_mps::Grid pnd_grid(particles.pnd_weight_radius, particles.position, pnd_valid, condition.dimension);
	while(timer.hasNextLoop()) {
		std::cout << timer.getCurrentTime() << std::endl;
		timer.limitCurrentDeltaTime(particles.getMaxSpeed(), condition);
		lap_grid.resetHash();
		particles.calculateTemporaryVelocity(condition.gravity, lap_grid, timer, condition);
		particles.moveExplicitly();
		pnd_grid.resetHash();
		particles.updateParticleNumberDensity(pnd_grid);
		particles.checkSurfaceParticles(condition.surface_parameter);
		timer.update();
	}
	particles.writeVtkFile("out.vtk", "test");
}