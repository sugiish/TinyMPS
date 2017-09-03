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
	Eigen::Matrix<bool, 1, Eigen::Dynamic> valid = particles.particle_types.array() == tiny_mps::ParticleType::NORMAL;
	tiny_mps::Grid lap_grid(particles.laplacian_pressure_weight_radius, particles.position, valid, condition.dimension);
	while(timer.hasNextLoop()) {
		std::cout << timer.getCurrentTime() << std::endl;
		lap_grid.resetHash();
		particles.moveParticlesExplicitly(condition.gravity, lap_grid, timer, condition);
		particles.velocity = particles.temporary_velocity;
		timer.update();
	}
	particles.writeVtkFile("out.vtk", "test");
}