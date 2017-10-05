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
    tiny_mps::Condition condition("./input/input_cavitation.data");
    tiny_mps::Particles particles("./input/input.grid", condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.02, -0.1, 0);
    Eigen::Vector3d maxpos(0.02, 0.05, 0);
    while(particles.nextLoop("./output/output_%1%.vtk", timer)) {
        particles.moveInflowParticles(timer);
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.updateTemporaryPosition(timer);
        particles.giveCollisionRepulsionForce();
        particles.updateTemporaryPosition(timer);
        particles.calculateTemporaryParticleNumberDensity();
        particles.checkSurfaceParticles();
        particles.solvePressurePoission(timer);
        particles.correctVelocity(timer);
        particles.updateTemporaryPosition(timer);
        particles.updateVelocityAndPosition();
        particles.removeOutsideParticles(minpos, maxpos);
    }
}