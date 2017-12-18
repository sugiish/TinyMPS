// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include <iostream>
#include <boost/format.hpp>
#include <Eigen/Core>
#include "condition.h"
#include "grid.h"
#include "particles.h"
#include "bubble_particles.h"

// Sample code using TinyMPS library.
int main(int argc, char* argv[]) {
  try {
    std::string output_path = "./output/";
    std::string input_data = "./input/nozzle.data";
    std::string input_grid = "./input/nozzle3.grid";
    if (argc >= 2) output_path = argv[1];
    // if (argc >= 3) input_data = argv[2];
    // if (argc >= 4) input_grid = argv[3];
    tiny_mps::Condition condition(input_data);
    if (argc >= 3) condition.inflow_velocity(1) = std::stod(argv[2]);
    my_mps::BubbleParticles particles(input_grid, condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.1, -2.1 * condition.average_distance, 0);
    // Eigen::Vector3d minpos(-0.1, -0.1, 0);
    Eigen::Vector3d maxpos(1.1, 2.1, 0);
    particles.initAverageGrid(minpos, maxpos);
    while(particles.nextLoop(output_path, timer)) {
      particles.moveInflowParticles(timer);
      particles.calculateTemporaryVelocity(condition.gravity, timer);
      particles.updateTemporaryPosition(timer);
      // particles.shiftParticles(2.1, 0.03);
      particles.giveCollisionRepulsionForce();
      particles.updateTemporaryPosition(timer);
      particles.calculateTemporaryParticleNumberDensity();
      particles.checkSurface2();
      particles.calculateModifiedParticleNumberDensity();
      particles.solvePressurePoissonDuan(timer);
      particles.correctVelocityDuan(timer);
      particles.updateTemporaryPosition(timer);
      particles.calculateAveragePressure();
      particles.updateAverageGrid(3.0e-4, timer);
      particles.calculateBubbles();
      particles.updateVelocityAndPosition();
      particles.removeOutsideParticles(minpos, maxpos);

      // particles.moveInflowParticles(timer);
      // // particles.shiftParticles(2.1, 0.03);
      // particles.giveCollisionRepulsionForce();
      // particles.calculateTemporaryVelocity(condition.gravity, timer);
      // particles.calculateTemporaryParticleNumberDensity();
      // particles.checkSurfaceParticles();
      // particles.solvePressurePoissonTanakaMasunaga(timer);
      // particles.correctVelocityTanakaMasunagaWithTensor(timer);
      // // particles.correctVelocityExplicitly(timer);
      // particles.updateTemporaryPosition(timer);
      // particles.updateVelocityAndPosition();
      // particles.removeOutsideParticles(minpos, maxpos);
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }
  return 1;
}