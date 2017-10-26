// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <stack>
#include <string>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
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
  Particles(int size, const Condition& condition);
  Particles(const std::string& path, const Condition& condition);
  Particles(const Particles& other);
  Particles& operator=(const Particles& other);
  virtual ~Particles();
  virtual void writeVtkFile(const std::string& path, const std::string& title) const;
  bool saveInterval(const std::string& path, const Timer& timer) const;
  bool nextLoop(const std::string& path, Timer& timer);
  bool checkNeedlessCalculation() const;
  virtual void extendStorage(int extra_size);
  int addParticle();
  virtual void setGhostParticle(int index);
  void removeOutsideParticles(const Eigen::Vector3d& minpos, const Eigen::Vector3d& maxpos);
  void removeFastParticles(double max_speed);
  void calculateTemporaryParticleNumberDensity();
  void updateParticleNumberDensity();
  void updateParticleNumberDensity(const Grid& grid);
  void updateVoxelRatio(int half_width, const Grid& grid);
  void moveInflowParticles(const Timer& timer);
  void calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer);
  void calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer, Grid& grid);
  void updateTemporaryPosition(const Timer& timer);
  void solvePressurePoission(const Timer& timer);
  void solvePressurePoissionOriginal(const Timer& timer, const Grid& grid);
  void solvePressurePoissionTanakaMasunaga(const Timer& timer);
  void correctVelocity(const Timer& timer);
  void correctVelocity(const Timer& timer, const Grid& grid);
  void correctVelocityExplicitly(const Timer& timer);
  void correctTanakaMasunagaVelocity(const Timer& timer);
  void correctVelocityWithTensor(const Timer& timer);
  void correctVelocityTanakaMasunagaWithTensor(const Timer& timer);
  void updateVelocityAndPosition();
  void checkSurfaceParticles();
  void checkSurfaceParticlesRemovingIsolated();
  void giveCollisionRepulsionForce();
  void giveCollisionRepulsionForce(double influence_ratio, double restitution_coefficient);
  void shiftParticles(double influence_ratio, double alpha);
  void showParticlesInfo();

  inline int getSize() const { return size; }
  inline int getDimension() const { return dimension; }
  inline double getMaxSpeed() const {
    Eigen::VectorXd moving = (particle_types.array() != ParticleType::GHOST).cast<double>().transpose();
    Eigen::VectorXd norms = velocity.colwise().norm();
    return (norms.array() * moving.array()).maxCoeff();
  }

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
  Eigen::VectorXd voxel_ratio;

 protected:
  virtual double weightForParticleNumberDensity(const Eigen::Vector3d& vec) const;
  virtual double weightForGradientPressure(const Eigen::Vector3d& vec) const;
  virtual double weightForLaplacianPressure(const Eigen::Vector3d& vec) const;
  virtual double weightForLaplacianViscosity(const Eigen::Vector3d& vec) const;

  const Condition& condition_;
  int size;
  const int dimension;
  std::stack<int> ghost_stack;
  double initial_particle_number_density;
  double laplacian_lambda_pressure;
  double laplacian_lambda_viscosity;
  double initial_neighbor_particles;
  double inflow_stride;

 private:
  using VectorXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
  void initialize(int particles_number);
  void readGridFile(const std::string& path, const Condition& condition);
  void setInitialParticleNumberDensity();
  void setLaplacianLambda();
  void solveConjugateGradient(Eigen::SparseMatrix<double> p_mat, Eigen::VectorXd source);

  static inline double weightStandard(const double distance, const double influence_radius) {
    if (distance < influence_radius) return (influence_radius / distance - 1.0);
    else return 0.0;
  }
  static inline double weightStandard(const Eigen::Vector3d& vec, const double influence_radius) {
    return weightStandard(vec.norm(), influence_radius);
  }
  static inline int weightCount(const double distance, const double influence_radius) {
    if (distance < influence_radius) return 1;
    else return 0;
  }
  static inline int weightCount(const Eigen::Vector3d& vec, const double influence_radius) {
    return weightCount(vec.norm(), influence_radius);
  }
};

} // namespace tiny_mps
#endif // MPS_PARTICLES_H_INCLUDED