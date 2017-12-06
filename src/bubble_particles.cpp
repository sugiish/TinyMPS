// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include "bubble_particles.h"
#define _USE_MATH_DEFINES
#include <cmath>

namespace my_mps{

BubbleParticles::BubbleParticles(const std::string& path, const tiny_mps::Condition& condition)
    : Particles(path, condition) {
  init_bubble_radius = cbrt((3 * condition.initial_void_fraction) / (4 * M_PI * condition.bubble_density * (1 - condition.initial_void_fraction)));
  bubble_radius = Eigen::VectorXd::Constant(getSize(), init_bubble_radius);
  void_fraction = Eigen::VectorXd::Constant(getSize(), condition.initial_void_fraction);
  free_surface_type = Eigen::VectorXi::Zero(getSize());

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
  for(int i = 0; i < size; ++i) {void calculateBubbles();
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
  ofs << std::endl;
  ofs << "SCALARS SourceTerm double" << std::endl;
  ofs << "LOOKUP_TABLE SourceTerm" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << source_term(i) << std::endl;
  }
  ofs << std::endl;
  ofs << "SCALARS VoxelsRatio double" << std::endl;
  ofs << "LOOKUP_TABLE VoxelsRatio" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << voxel_ratio(i) << std::endl;
  }

  // extended
  ofs << std::endl;
  ofs << "SCALARS BubbleRadius double" << std::endl;
  ofs << "LOOKUP_TABLE BubbleRadius" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << bubble_radius(i) << std::endl;
  }
  std::cout << "Succeed in writing vtk file: " << path << std::endl;
  ofs << std::endl;
  ofs << "SCALARS VoidFraction double" << std::endl;
  ofs << "LOOKUP_TABLE VoidFraction" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << void_fraction(i) << std::endl;
  }
  ofs << std::endl;
  ofs << "SCALARS FreeSurfaceType int" << std::endl;
  ofs << "LOOKUP_TABLE FreeSurfaceType" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << free_surface_type(i) << std::endl;
  }
  
  std::cout << "Succeed in writing vtk file: " << path << std::endl;
}

void BubbleParticles::extendStorage(int extra_size) {
  int size = getSize();
  Particles::extendStorage(extra_size);
  bubble_radius.conservativeResize(size + extra_size);
  bubble_radius.segment(size, extra_size).setZero();
  void_fraction.conservativeResize(size + extra_size);
  void_fraction.segment(size, extra_size).setZero();
  free_surface_type.conservativeResize(size + extra_size);
  free_surface_type.segment(size, extra_size).setZero();
}

void BubbleParticles::setGhostParticle(int index) {
  Particles::setGhostParticle(index);
  bubble_radius(index) = 0.0;
  void_fraction(index) = condition_.initial_void_fraction;
}

void BubbleParticles::checkSurface(double shift){
  using namespace tiny_mps;
  for(int i = 0; i < getSize(); ++i) {
    if (particle_types(i) == ParticleType::NORMAL || particle_types(i) == ParticleType::WALL
        || particle_types(i) == ParticleType::INFLOW) {
      if (particle_number_density(i) < condition_.surface_threshold_pnd * initial_particle_number_density
              && neighbor_particles(i) < condition_.surface_threshold_number * initial_neighbor_particles
              && temporary_position(1, i) < shift) {
        boundary_types(i) = BoundaryType::SURFACE;
      } else {
        boundary_types(i) = BoundaryType::INNER;
      }
    } else {
      boundary_types(i) = BoundaryType::OTHERS;
    }
  }
}

void BubbleParticles::calculateBubbles() {
  for (int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if (particle_types(i_particle) == tiny_mps::ParticleType::NORMAL) {
      double del_p = (condition_.vapor_pressure - condition_.head_pressure) - pressure(i_particle);
      if (del_p > 0) bubble_radius(i_particle) += sqrt(2 * abs(del_p) / (3 * condition_.mass_density));
      else bubble_radius(i_particle) -= sqrt(2 * abs(del_p) / (3 * condition_.mass_density));
      if (bubble_radius(i_particle) > condition_.average_distance) bubble_radius(i_particle) = condition_.average_distance;
      if (bubble_radius(i_particle) < 0) bubble_radius(i_particle) = 0;

      double bubble_vol = 4 * M_PI * condition_.bubble_density * bubble_radius(i_particle) * bubble_radius(i_particle) * bubble_radius(i_particle) / 3;
      void_fraction(i_particle) = bubble_vol / (1 + bubble_vol);
      if (void_fraction(i_particle) < condition_.min_void_fraction) void_fraction(i_particle) = condition_.min_void_fraction;
      if (void_fraction(i_particle) > 0.5) void_fraction(i_particle) = 0.5;
    }
  }
}


void BubbleParticles::solvePressurePoisson(const tiny_mps::Timer& timer) {
  using namespace tiny_mps;
  Grid grid(condition_.laplacian_pressure_weight_radius, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
  using T = Eigen::Triplet<double>;
  double lap_r = grid.getGridWidth();
  int n_size = (int)(std::pow(lap_r * 2, dimension));
  double delta_time = timer.getCurrentDeltaTime();
  Eigen::SparseMatrix<double> p_mat(size, size);
  source_term.setZero();
  std::vector<T> coeffs(size * n_size);
  std::vector<int> neighbors(n_size * 2);
  for (int i_particle = 0; i_particle < size; ++i_particle) {
    if (boundary_types(i_particle) == BoundaryType::OTHERS) {
      coeffs.push_back(T(i_particle, i_particle, 1.0));
      continue;
    } else if (boundary_types(i_particle) == BoundaryType::SURFACE) {
      coeffs.push_back(T(i_particle, i_particle, 1.0));
      continue;
    }
    grid.getNeighbors(i_particle, neighbors);
    double sum = 0.0;
    double div_vel = 0.0;
    for (int j_particle : neighbors) {
      if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
      Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
      double mat_ij = weightForLaplacianPressure(r_ij) * 2 * dimension
              / (laplacian_lambda_pressure * initial_particle_number_density);
      sum -= mat_ij;
      div_vel += (temporary_velocity.col(j_particle) - temporary_velocity.col(i_particle)).dot(r_ij)
              * weightForLaplacianPressure(r_ij) * condition_.dimension / (r_ij.squaredNorm() * initial_particle_number_density);
      if (boundary_types(j_particle) == BoundaryType::INNER) {
        coeffs.push_back(T(i_particle, j_particle, mat_ij));
      }
    }
    sum -= condition_.weak_compressibility * condition_.mass_density / (delta_time * delta_time);
    coeffs.push_back(T(i_particle, i_particle, sum));
    double initial_pnd_i = initial_particle_number_density * (1 - void_fraction(i_particle));
    source_term(i_particle) = div_vel * condition_.mass_density * condition_.relaxation_coefficient_vel_div / delta_time
                - (particle_number_density(i_particle) - initial_pnd_i)
                * condition_.relaxation_coefficient_pnd * condition_.mass_density / (delta_time * delta_time * initial_particle_number_density);
  }
  p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
  solveConjugateGradient(p_mat);
}

}