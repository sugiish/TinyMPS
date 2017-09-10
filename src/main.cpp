// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include <iostream>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "condition.h"
#include "grid.h"
#include "particles.h"

// Sample code using TinyMPS library.
int main() {
    tiny_mps::Condition condition("./input/input.data");
    tiny_mps::Particles particles("./input/dam.grid", condition);
    tiny_mps::Timer timer(condition);
    particles.saveInterval("./output/output_%04d.vtk", timer);
    while(timer.hasNextLoop()) {
        timer.limitCurrentDeltaTime(particles.getMaxSpeed(), condition);
        timer.printTimeInfo();
        std::cout << boost::format("Max velocity: %f") % particles.getMaxSpeed() << std::endl;
        if (particles.checkNeedlessCalculation()) {
            std::cerr << "Error: All particles have become ghost." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (timer.isUnderMinDeltaTime()) {
            std::cerr << "Error: Delta time has become so small." << std::endl;
            exit(EXIT_FAILURE);
        }
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.calculateTemporaryParticleNumberDensity(condition);
        particles.checkSurfaceParticles();
        particles.solvePressurePoission(timer);
        particles.correctVelocity(timer);
        particles.saveInterval("./output/output_%04d.vtk", timer);
        std::cout << std::endl;
        timer.update();
    }
    std::cout << "Succeed in simulation." << std::endl;
}