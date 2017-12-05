// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include <iostream>
#include <boost/format.hpp>
#include <Eigen/Core>
#include "condition.h"
#include "grid.h"
#include "particles.h"

// Sample code using TinyMPS library.
int main(int argc, char* argv[]) {
  try {
    std::string output_path = "./output/";
    std::string input_data = "./input/input.data";
    std::string input_grid = "./input/dam.grid";
    if (argc >= 2) output_path = argv[1];
    if (argc >= 3) input_data = argv[2];
    if (argc >= 4) input_grid = argv[3];
    output_path += "output_%1%.vtk";
    tiny_mps::Condition condition(input_data);
    tiny_mps::Particles particles(input_grid, condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.1, -0.1, 0);
    Eigen::Vector3d maxpos(1.1, 2.1, 0);
    while(particles.nextLoop(output_path, timer)) {
      particles.calculateTemporaryVelocity(condition.gravity, timer);
      particles.updateTemporaryPosition(timer);
      particles.giveCollisionRepulsionForce();
      particles.updateTemporaryPosition(timer);
      particles.calculateTemporaryParticleNumberDensity();
      particles.checkSurfaceParticles();
      particles.solvePressurePoisson(timer);
      particles.setZeroOnNegativePressure();
      particles.correctVelocity(timer);
      particles.updateTemporaryPosition(timer);
      particles.updateVelocityAndPosition();
      particles.removeOutsideParticles(minpos, maxpos);
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }
  return 1;
}