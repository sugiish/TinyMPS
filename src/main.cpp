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
    tiny_mps::Particles particles("./input/dam.grid", condition);
    tiny_mps::Timer timer(condition);
    particles.writeVtkFile("first.vtk", "test");
    while(timer.hasNextLoop()) {
        std::cout << timer.getCurrentTime() << std::endl;
        timer.limitCurrentDeltaTime(particles.getMaxSpeed(), condition);
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.calculateTemporaryParticleNumberDensity(condition);
        particles.checkSurfaceParticles();
        particles.solvePressurePoission(timer);
        particles.correctVelocity(timer);
        timer.update();
    }
    particles.writeVtkFile("out.vtk", "test");
}