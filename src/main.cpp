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
    while(particles.next("./output/output_%04d.vtk", timer, condition)) {
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.calculateTemporaryParticleNumberDensity(condition);
        particles.checkSurfaceParticles();
        particles.solvePressurePoission(timer);
        particles.correctVelocity(timer);
    }
}