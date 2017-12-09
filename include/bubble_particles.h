// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_BUBBLE_PARTICLES_H_INCLUDED
#define MPS_BUBBLE_PARTICLES_H_INCLUDED

#include "particles.h"
#include "timer.h"

namespace my_mps{

enum SurfaceLayer {
  INNER = 0,
  OUTER_SURFACE = 1,
  INNER_SURFACE = 2,
  OTHERS = -1
};

class BubbleParticles : public tiny_mps::Particles {
 public:
  BubbleParticles(const std::string& path, const tiny_mps::Condition& condition);
  BubbleParticles(const tiny_mps::Particles& other) = delete;
  BubbleParticles& operator=(const tiny_mps::Particles& other) = delete;
  virtual ~BubbleParticles() {};
  void writeVtkFile(const std::string& path, const std::string& title) const;
  void extendStorage(int extra_size);
  void setGhostParticle(int index);
  void calculateBubbles();
  void calculateModifiedParticleNumberDensity();
  void solvePressurePoisson(const tiny_mps::Timer& timer);
  void solvePressurePoissonDuan(const tiny_mps::Timer& timer);
  void checkSurface();
  void checkSurface2();
  void correctVelocityDuan(const tiny_mps::Timer& timer);

 private:
  Eigen::Matrix3Xd normal_vector;
  Eigen::VectorXd modified_pnd;
  Eigen::VectorXd bubble_radius;
  Eigen::VectorXd void_fraction;
  Eigen::VectorXi free_surface_type;
  double init_bubble_radius;
};

}

#endif // MPS_BUBBLE_PARTICLES_H_INCLUDED