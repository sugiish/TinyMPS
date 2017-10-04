// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <stack>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "condition.h"
#include "grid.h"
#include "timer.h"

namespace tiny_mps {

// Particle type
enum ParticleType {
    NORMAL = 0,
    WALL = 2,
    DUMMY_WALL = 3,
    INFLOW = 4,
    DUMMY_INFLOW = 5,
    GHOST = -1
};

// Used for solving Poisson's equation.
enum BoundaryType {
    INNER = 0,
    SURFACE = 1,
    OTHERS = -1
};

// Holds data on particles and manipulates them.
class Particles {
public:
    Particles(const std::string& path, const Condition& condition);
    virtual ~Particles();
    int writeVtkFile(const std::string& path, const std::string& title) const;
    bool saveInterval(const std::string& path, const Timer& timer) const;
    bool nextLoop(const std::string& path, Timer& timer);
    bool checkNeedlessCalculation() const;
    int addParticle();
    void setGhostParticle(int index);
    void removeOutsideParticles(const Eigen::Vector3d& minpos, const Eigen::Vector3d& maxpos);
    void removeFastParticles(double max_speed);
    void calculateTemporaryParticleNumberDensity();
    void updateParticleNumberDensity();
    void updateParticleNumberDensity(const Grid& grid);
    void moveInflowParticles(const Timer& timer);
    void calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer);
    void calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer, Grid& grid);
    void updateTemporaryPosition(const Timer& timer);
    void solvePressurePoission(const Timer& timer, const Grid& grid);
    void solvePressurePoission(const Timer& timer);
    void solveTanakaMasunagaPressurePoission(const Timer& timer);
    void solvePressurePoissionWithTanakaMasunaga(const Timer& timer);
    void correctVelocity(const Timer& timer);
    void correctVelocity(const Timer& timer, const Grid& grid);
    void correctVelocityExplicitly(const Timer& timer);
    void correctTanakaMasunagaVelocity(const Timer& timer);
    void updateVelocityAndPosition();
    void checkSurfaceParticles();
    void checkSurfaceParticles(double surface_parameter);
    void checkSurfaceParticlesWithTanakaMasunaga();
    void checkTanakaMasunagaSurfaceParticles(double surface_parameter);
    void giveCollisionRepulsionForce();
    void giveCollisionRepulsionForce(double influence_ratio, double restitution_coefficient);
    void showParticlesInfo();
    inline double getMaxSpeed() const {
        Eigen::VectorXd moving = (particle_types.array() != ParticleType::GHOST).cast<double>().transpose();
        Eigen::VectorXd norms = velocity.colwise().norm();
        return (norms.array() * moving.array()).maxCoeff();
    }
    inline int getSize() const { return size; }
    inline int getDimension() const { return dimension; }

    Eigen::Matrix3Xd position;
    Eigen::Matrix3Xd velocity;
    Eigen::VectorXd pressure;
    Eigen::VectorXd particle_number_density;
    Eigen::Matrix3Xd temporary_position;
    Eigen::Matrix3Xd temporary_velocity;
    Eigen::Matrix3Xd correction_velocity;
    Eigen::VectorXi particle_types;
    Eigen::VectorXi boundary_types;
    Eigen::VectorXi neighbor_particles;

private:
    using VectorXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
    Particles(const Particles&);
    Particles& operator=(const Particles&);
    void initialize(int particles_number);
    int readGridFile(const std::string& path, const Condition& condition);
    void setInitialParticleNumberDensity(int index);
    void calculateLaplacianLambda(int index, const Condition& condition);
    void solveConjugateGradient(Eigen::SparseMatrix<double> p_mat, Eigen::VectorXd source);
    inline double weightFunction(const Eigen::Vector3d& vec, double influence_radius) const {
        double r = vec.norm();
        if(r < influence_radius) return (influence_radius / r - 1.0);
        else return 0.0;
    }

    const Condition& condition_;
    std::stack<int> ghost_stack;
    int size;
    int dimension;
    double initial_particle_number_density;
    double laplacian_lambda_pressure;
    double laplacian_lambda_viscosity;
    double initial_neighbor_particles;
    double inflow_stride;
};

} // namespace tiny_mps
#endif // MPS_PARTICLES_H_INCLUDED