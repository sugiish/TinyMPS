// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

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
    bool checkNeedlessCalculation() const;
    bool nextLoop(const std::string& path, Timer& timer, const Condition& condition);
    void calculateTemporaryParticleNumberDensity(const Condition& condition);
    void updateParticleNumberDensity(const Condition& condition);
    void updateParticleNumberDensity(const Grid& grid);
    void calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer);
    void calculateTemporaryVelocity(const Eigen::Vector3d& force, Grid& grid, const Timer& timer, const Condition& condition);
    void solvePressurePoission(const Grid& grid, const Timer& timer, const Condition& condition);
    void solvePressurePoission(const Timer& timer);
    void correctVelocity(const Timer& timer);
    void correctVelocity(const Grid& grid, const Timer& timer, const Condition& condition);
    void checkSurfaceParticles();
    void checkSurfaceParticles(double surface_parameter);
    int writeVtkFile(const std::string& path, const std::string& title) const;
    bool saveInterval(const std::string& path, const Timer& timer) const;
    inline double getMaxSpeed() const {
        Eigen::VectorXd moving = (particle_types.array() != ParticleType::GHOST).cast<double>().transpose();
        Eigen::VectorXd norms = velocity.colwise().norm();
        return (norms.array() * moving.array()).maxCoeff();
    }
    inline int getSize() const { return size; }
    inline int getDimension() const { return dimension; }

private:
    using VectorXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
    void initialize(int particles_number);
    int readGridFile(const std::string& path, int dimension);
    void setInitialParticleNumberDensity(int index);
    void calculateLaplacianLambda(int index, const Grid& grid);
    inline double weightFunction(const Eigen::Vector3d& vec, double influence_radius) const {
        double r = vec.norm();
        if(r < influence_radius) return (influence_radius / r - 1.0);
        else return 0.0;
    }

    const Condition& condition_;
    Eigen::Matrix3Xd position;
    Eigen::Matrix3Xd velocity;
    Eigen::VectorXd pressure;
    Eigen::VectorXd particle_number_density;
    Eigen::Matrix3Xd temporary_position;
    Eigen::Matrix3Xd temporary_velocity;
    Eigen::Matrix3Xd correction_velocity;
    Eigen::VectorXi particle_types;
    Eigen::VectorXi boundary_types;
    int size;
    int dimension;
    double initial_particle_number_density;
    double laplacian_lambda;
};

} // namespace tiny_mps
#endif // MPS_PARTICLES_H_INCLUDED