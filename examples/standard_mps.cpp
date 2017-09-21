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
    tiny_mps::Condition condition("./input/input_standard.data");
    tiny_mps::Particles particles("./input/dam.grid", condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.1, -0.1, 0);
    Eigen::Vector3d maxpos(0.8, 0.8, 0);
    while(particles.nextLoop("./output/output_%04d.vtk", timer)) {
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.giveCollisionRepulsion(timer);
        particles.calculateTemporaryParticleNumberDensity();
        particles.checkSurfaceParticles();
        particles.solvePressurePoission(timer);
        particles.correctVelocity(timer);
        particles.updateFromTemporary();
        particles.removeOutsideParticles(minpos, maxpos);
    }
}