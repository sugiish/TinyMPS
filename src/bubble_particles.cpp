// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include "bubble_particles.h"

namespace my_mps{

BubbleParticles::BubbleParticles(const std::string& path, const tiny_mps::Condition& condition)
    : Particles(path, condition) {
  bubble_radius = Eigen::VectorXd::Zero(getSize());
}

void BubbleParticles::writeVtkFile(const std::string& path, const std::string& title) const {
  int size = getSize();
  std::ofstream ofs(path);
  if(ofs.fail()) {
    std::cerr << "Error: in writeVtkFile() in particles.cpp." << std::endl;
    throw std::ios_base::failure("Error: in writeVtkFile() in particles.cpp.");
  }
  ofs << "# vtk DataFile Version 2.0" << std::endl;
  ofs << title << std::endl;
  ofs << "ASCII" << std::endl;
  ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
  ofs << std::endl;
  ofs << "POINTS " << size << " double" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << position(0, i) << " " << position(1, i) << " " << position(2, i) << std::endl;
  }
  ofs << std::endl;
  ofs << "CELLS " << size << " " << size * 2 << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << 1 << " " << i << std::endl;
  }
  ofs << std::endl;
  ofs << "CELL_TYPES " << size << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << 1 << std::endl;
  }
  ofs << std::endl;
  ofs << "POINT_DATA " << size << std::endl;
  ofs << "SCALARS Pressure double" << std::endl;
  ofs << "LOOKUP_TABLE Pressure" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << pressure(i) << std::endl;
  }
  ofs << std::endl;
  ofs << "VECTORS Velocity double" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << velocity(0, i) << " " << velocity(1, i) << " " << velocity(2, i) << std::endl;
  }
  ofs << std::endl;
  ofs << "SCALARS Type int" << std::endl;
  ofs << "LOOKUP_TABLE Type" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << particle_types(i) << std::endl;
  }
  ofs << std::endl;
  ofs << "SCALARS ParticleNumberDensity double" << std::endl;
  ofs << "LOOKUP_TABLE ParticleNumberDensity" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << particle_number_density(i) << std::endl;
  }
  ofs << std::endl;
  ofs << "SCALARS NeighborParticles int" << std::endl;
  ofs << "LOOKUP_TABLE NeighborParticles" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << neighbor_particles(i) << std::endl;
  }
  ofs << std::endl;
  ofs << "SCALARS BoundaryCondition int" << std::endl;
  ofs << "LOOKUP_TABLE BoundaryCondition" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << boundary_types(i) << std::endl;
  }
  ofs << std::endl;
  ofs << "VECTORS CorrectionVelocity double" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << correction_velocity(0, i) << " " << correction_velocity(1, i) << " " << correction_velocity(2, i) << std::endl;
  }

  // extended
  ofs << std::endl;
  ofs << "SCALARS BubbleRadius int" << std::endl;
  ofs << "LOOKUP_TABLE BubbleRadius" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << bubble_radius(i) << std::endl;
  }
  std::cout << "Succeed in writing vtk file: " << path << std::endl;
}

void BubbleParticles::extendStorage(int extra_size) {
  int size = getSize();
  Particles::extendStorage(extra_size);
  bubble_radius.conservativeResize(size + extra_size);
  bubble_radius.segment(size, extra_size) = Eigen::VectorXd::Zero(extra_size);
}

void BubbleParticles::setGhostParticle(int index) {
  Particles::setGhostParticle(index);
  bubble_radius(index) = 0.0;
}

}