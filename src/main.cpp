// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include <iostream>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "condition.h"
#include "grid.h"
#include "particles.h"
/*
// Sample code using TinyMPS library.
int main() {
    tiny_mps::Condition condition("./input/input.data");
    tiny_mps::Particles particles("./input/dam_tm.grid", condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.1, -0.1, 0);
    Eigen::Vector3d maxpos(1.1, 1.6, 0);
    while(particles.nextLoop("./output/output_%04d.vtk", timer, condition)) {
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.giveCollisionRepulsion(0.85, 0.2, timer, condition);
        particles.calculateTemporaryParticleNumberDensity(condition);
        particles.checkSurfaceParticlesWithTanakaMasunaga(condition);
        particles.solvePressurePoission(timer);
        particles.correctVelocity(timer);
        particles.updateFromTemporary();
        particles.removeOutsideParticles(minpos, maxpos);
    }
}

int main() {
    tiny_mps::Condition condition("./input/input.data");
    tiny_mps::Particles particles("./input/dam_tm.grid", condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.1, -0.1, 0);
    Eigen::Vector3d maxpos(1.1, 1.6, 0);
    while(particles.nextLoop("./output/output_%04d.vtk", timer, condition)) {
        particles.giveCollisionRepulsion(0.85, 0.2, timer, condition);
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        particles.calculateTemporaryParticleNumberDensity(condition);
        particles.checkSurfaceParticlesWithTanakaMasunaga(condition);
        particles.solvePressurePoission(timer);
        particles.correctVelocity(timer);
        // particles.giveCollisionRepulsion(0.85, 0.2, timer, condition);
        particles.updateFromTemporary();
        particles.removeOutsideParticles(minpos, maxpos);
        particles.removeFastParticles(10);
    }
}*/

// Sample code using TinyMPS library.
int main() {
    tiny_mps::Condition condition("./input/input.data");
    tiny_mps::Particles particles("./input/hydrostatic.grid", condition);
    tiny_mps::Timer timer(condition);
    Eigen::Vector3d minpos(-0.1, -0.1, 0);
    Eigen::Vector3d maxpos(1.1, 1.6, 0);
    while(particles.nextLoop("./output/output_%04d.vtk", timer, condition)) {
        particles.giveCollisionRepulsion(0.85, 0.2, timer, condition);
        particles.updateParticleNumberDensity(condition);
        particles.checkSurfaceParticlesWithTanakaMasunaga(condition);
        // particles.checkTanakaMasunagaSurfaceParticles(condition.tanaka_masunaga_beta);
        particles.calculateTemporaryVelocity(condition.gravity, timer);
        // particles.saveInterval("./output/pre_cor_%04d.vtk", timer);
        // particles.calculateTemporaryParticleNumberDensity(condition);
        // particles.checkSurfaceParticles();
        particles.solveTanakaMasunagaPressurePoission(timer, condition);
        particles.correctVelocityExplicitly(timer);
        particles.updateFromTemporary();
        particles.removeOutsideParticles(minpos, maxpos);
        // particles.removeFastParticles(100);
    }
}