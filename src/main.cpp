// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include <iostream>
#include <Eigen/Dense>
#include "condition.h"
#include "grid.h"
#include "particles.h"

// Sample code using TinyMPS library.
int main() {
    tiny_mps::Condition condition("./input/input.data");
    tiny_mps::Particles particles("./input/dambreak.grid", condition);
    tiny_mps::Timer timer(condition);
    Eigen::Matrix<bool, 1, Eigen::Dynamic> lap_valid = particles.particle_types.array() == tiny_mps::ParticleType::NORMAL;
    // tiny_mps::Grid pnd_grid(particles.pnd_weight_radius, particles.position, pnd_valid, condition.dimension);
    particles.writeVtkFile("first.vtk", "test");
    while(timer.hasNextLoop()) {
        std::cout << timer.getCurrentTime() << std::endl;
        timer.limitCurrentDeltaTime(particles.getMaxSpeed(), condition);
        tiny_mps::Grid lap_grid(particles.laplacian_pressure_weight_radius, particles.position, lap_valid, condition.dimension);
        particles.calculateTemporaryVelocity(condition.gravity, lap_grid, timer, condition);
        particles.moveExplicitly();
        Eigen::Matrix<bool, 1, Eigen::Dynamic> pnd_valid = particles.particle_types.array() != tiny_mps::ParticleType::GHOST;
        tiny_mps::Grid pnd_grid(particles.pnd_weight_radius, particles.position, pnd_valid, condition.dimension);
        particles.updateParticleNumberDensity(pnd_grid);
        particles.checkSurfaceParticles(condition.surface_parameter);
        Eigen::Matrix<bool, 1, Eigen::Dynamic> solve_valid = (particles.boundary_types.array() == tiny_mps::BoundaryType::INNER)
            || (particles.boundary_types.array() == tiny_mps::BoundaryType::SURFACE);
        tiny_mps::Grid solve_grid(particles.laplacian_pressure_weight_radius, particles.position, solve_valid, condition.dimension);
        particles.solvePressurePoission(solve_grid, timer, condition);
        tiny_mps::Grid grad_grid(condition.gradient_influence * condition.average_distance, particles.position, solve_valid, condition.dimension);
        // particles.advectVelocity(solve_grid, timer, condition);
        timer.update();
    }
    particles.writeVtkFile("out.vtk", "test");
}