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
    tiny_mps::Condition condition("../input/input.data");
    tiny_mps::Particles particles("../input/dam_tm.grid", condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.1, -0.1, 0);
    Eigen::Vector3d maxpos(1.1, 2.1, 0);
    while(particles.nextLoop("../output/output_%04d.vtk", timer, condition)) {
        particles.updateParticleNumberDensity(condition);
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.checkTanakaMasunagaSurfaceParticles(condition.tanaka_masunaga_beta);
        // particles.checkSurfaceParticles();
        particles.solveTanakaMasunagaPressurePoission(timer, condition);
        particles.correctTanakaMasunagaVelocity(timer, condition);
        particles.updateFromTemporary();
        particles.removeOutsideParticles(minpos, maxpos);
        particles.removeFastParticles(170);
    }
}