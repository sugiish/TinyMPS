// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include "bubble_particles.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/LU>

namespace my_mps{

BubbleParticles::BubbleParticles(const std::string& path, const tiny_mps::Condition& condition)
    : Particles(path, condition) {
  init_bubble_radius = cbrt((3 * condition.initial_void_fraction) / (4 * M_PI * condition.bubble_density * (1 - condition.initial_void_fraction)));
  normal_vector = Eigen::Matrix3Xd::Zero(3, getSize());
  modified_pnd = Eigen::VectorXd::Zero(getSize());
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
  ofs << "VECTORS NormalVector double" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << normal_vector(0, i) << " " << normal_vector(1, i) << " " << normal_vector(2, i) << std::endl;
  }
  ofs << std::endl;
  ofs << "SCALARS BubbleRadius double" << std::endl;
  ofs << "LOOKUP_TABLE BubbleRadius" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << bubble_radius(i) << std::endl;
  }
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
  ofs << std::endl;
  ofs << "SCALARS ModifiedParticleNumberDensity double" << std::endl;
  ofs << "LOOKUP_TABLE ModifiedParticleNumberDensity" << std::endl;
  for(int i = 0; i < size; ++i) {
    ofs << modified_pnd(i) << std::endl;
  }

  std::cout << "Succeed in writing vtk file: " << path << std::endl;
}

void BubbleParticles::extendStorage(int extra_size) {
  int size = getSize();
  Particles::extendStorage(extra_size);
  normal_vector.conservativeResize(3, size + extra_size);
  normal_vector.block(0, size, 3, extra_size).setZero();
  modified_pnd.conservativeResize(size + extra_size);
  modified_pnd.segment(size, extra_size).setZero();
  bubble_radius.conservativeResize(size + extra_size);
  bubble_radius.segment(size, extra_size).setZero();
  void_fraction.conservativeResize(size + extra_size);
  void_fraction.segment(size, extra_size).setZero();
  free_surface_type.conservativeResize(size + extra_size);
  free_surface_type.segment(size, extra_size).setZero();
}

void BubbleParticles::setGhostParticle(int index) {
  Particles::setGhostParticle(index);
  normal_vector.col(index).setZero();
  modified_pnd(index) = 0.0;
  bubble_radius(index) = 0.0;
  void_fraction(index) = condition_.initial_void_fraction;
  free_surface_type(index) = SurfaceLayer::OTHERS;
}

void BubbleParticles::checkSurface(){
  // First step.
  using namespace tiny_mps;
  for(int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if (particle_types(i_particle) == ParticleType::NORMAL || particle_types(i_particle) == ParticleType::WALL
        || particle_types(i_particle) == ParticleType::INFLOW) {
      if (particle_number_density(i_particle) < condition_.surface_threshold_pnd * initial_particle_number_density
              && neighbor_particles(i_particle) < condition_.surface_threshold_number * initial_neighbor_particles) {
              // && temporary_position(1, i) < shift) {
        boundary_types(i_particle) = BoundaryType::SURFACE;
        free_surface_type(i_particle) = SurfaceLayer::OUTER_SURFACE;
      } else {
        boundary_types(i_particle) = BoundaryType::INNER;
        free_surface_type(i_particle) = SurfaceLayer::INNER;
      }
    } else {
      boundary_types(i_particle) = BoundaryType::OTHERS;
      free_surface_type(i_particle) = SurfaceLayer::OTHERS;
    }
  }
  // Second step.
  Grid grid(condition_.pnd_weight_radius, temporary_position, particle_types.array() != ParticleType::GHOST, dimension);
  const double root2 = std::sqrt(2);
  normal_vector.setZero();
  for(int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if(boundary_types(i_particle) == BoundaryType::SURFACE) {
      Grid::Neighbors neighbors;
      grid.getNeighbors(i_particle, neighbors);
      if (neighbors.empty()) continue;
      for (int j_particle : neighbors) {
        Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
        normal_vector.col(i_particle) += r_ij.normalized() * weightForParticleNumberDensity(r_ij);
      }
      normal_vector.col(i_particle) /= particle_number_density(i_particle);
      for (int j_particle : neighbors) {
        Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
        if (r_ij.norm() >= root2 * condition_.average_distance
            && (temporary_position.col(i_particle) + condition_.average_distance * normal_vector.col(i_particle).normalized() - temporary_position.col(j_particle)).norm() < condition_.average_distance) {
          boundary_types(i_particle) = BoundaryType::INNER;
          free_surface_type(i_particle) = SurfaceLayer::INNER;
          break;
        }
        if (r_ij.norm() < root2 * condition_.average_distance
            && r_ij.normalized().dot(normal_vector.col(i_particle).normalized()) > 1.0 / root2) {
          boundary_types(i_particle) = BoundaryType::INNER;
          free_surface_type(i_particle) = SurfaceLayer::INNER;
          break;
        }
      }
    }
  }
  // Third step.
  Grid judge_inner_surface(condition_.average_distance * condition_.secondary_surface_eta, temporary_position, particle_types.array() != ParticleType::GHOST, dimension);
  for (int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if (boundary_types(i_particle) == BoundaryType::SURFACE) {
      Grid::Neighbors neighbors;
      judge_inner_surface.getNeighbors(i_particle, neighbors);
      if (neighbors.empty()) continue;
      for (int j_particle : neighbors) {
        if (boundary_types(j_particle) == BoundaryType::INNER && free_surface_type(j_particle) == SurfaceLayer::INNER)
          free_surface_type(j_particle) = SurfaceLayer::INNER_SURFACE;
      }
    }
  }
}

void BubbleParticles::checkSurface2(){
  // First step.
  using namespace tiny_mps;
  for(int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if (particle_types(i_particle) == ParticleType::NORMAL || particle_types(i_particle) == ParticleType::WALL
        || particle_types(i_particle) == ParticleType::INFLOW) {
      if (particle_number_density(i_particle) < condition_.surface_threshold_pnd * initial_particle_number_density
              && neighbor_particles(i_particle) < condition_.surface_threshold_number * initial_neighbor_particles
              && temporary_position(1, i_particle) < condition_.pnd_weight_radius) {
        boundary_types(i_particle) = BoundaryType::SURFACE;
        free_surface_type(i_particle) = SurfaceLayer::OUTER_SURFACE;
      } else {
        boundary_types(i_particle) = BoundaryType::INNER;
        free_surface_type(i_particle) = SurfaceLayer::INNER;
      }
    } else {
      boundary_types(i_particle) = BoundaryType::OTHERS;
      free_surface_type(i_particle) = SurfaceLayer::OTHERS;
    }
  }
  // Second step.
  Grid grid(condition_.pnd_weight_radius, temporary_position, particle_types.array() != ParticleType::GHOST, dimension);
  normal_vector.setZero();
  for(int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if(boundary_types(i_particle) == BoundaryType::SURFACE) {
      Grid::Neighbors neighbors;
      grid.getNeighbors(i_particle, neighbors);
      if (neighbors.empty()) continue;
      for (int j_particle : neighbors) {
        Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
        normal_vector.col(i_particle) += r_ij.normalized() * weightForParticleNumberDensity(r_ij);
      }
      normal_vector.col(i_particle) /= particle_number_density(i_particle);
    }
  }
  // Third step.
  Grid judge_inner_surface(condition_.average_distance * condition_.secondary_surface_eta, temporary_position, particle_types.array() != ParticleType::GHOST, dimension);
  for (int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if (boundary_types(i_particle) == BoundaryType::SURFACE) {
      Grid::Neighbors neighbors;
      judge_inner_surface.getNeighbors(i_particle, neighbors);
      if (neighbors.empty()) continue;
      for (int j_particle : neighbors) {
        if (boundary_types(j_particle) == BoundaryType::INNER && free_surface_type(j_particle) == SurfaceLayer::INNER)
          free_surface_type(j_particle) = SurfaceLayer::INNER_SURFACE;
      }
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

void BubbleParticles::calculateModifiedParticleNumberDensity() {
  using namespace tiny_mps;
  Grid grid(condition_.average_distance, temporary_position, particle_types.array() != ParticleType::GHOST, condition_.dimension);
  Eigen::Vector3d l0_vec(condition_.average_distance, 0.0, 0.0);
  for (int i_particle = 0; i_particle < getSize(); ++i_particle) {
    if (particle_types(i_particle) == ParticleType::GHOST) modified_pnd(i_particle) = 0.0;
    double n_hat = initial_particle_number_density;
    Grid::Neighbors neighbors;
    grid.getNeighbors(i_particle, neighbors);
    if (neighbors.empty()) continue;
    for (int j_particle : neighbors) {
      Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
      n_hat += weightForParticleNumberDensity(r_ij) - weightForParticleNumberDensity(l0_vec);
    }
    modified_pnd(i_particle) = std::max(particle_number_density(i_particle), n_hat);
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

void BubbleParticles::solvePressurePoissonDuan(const tiny_mps::Timer& timer) {
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
    if (free_surface_type(i_particle) == SurfaceLayer::INNER_SURFACE) {
      sum -= (modified_pnd(i_particle) - particle_number_density(i_particle)) * 2 * dimension / (laplacian_lambda_pressure * initial_particle_number_density);
      coeffs.push_back(T(i_particle, i_particle, sum));
      // double initial_pnd_i = initial_particle_number_density * (1 - void_fraction(i_particle));
      source_term(i_particle) = div_vel * condition_.mass_density * condition_.relaxation_coefficient_vel_div / delta_time
                - (modified_pnd(i_particle) - initial_particle_number_density)
                * condition_.relaxation_coefficient_pnd * condition_.mass_density / (delta_time * delta_time * initial_particle_number_density);
    } else {
      coeffs.push_back(T(i_particle, i_particle, sum));
      // double initial_pnd_i = initial_particle_number_density * (1 - void_fraction(i_particle));
      source_term(i_particle) = div_vel * condition_.mass_density * condition_.relaxation_coefficient_vel_div / delta_time
                - (particle_number_density(i_particle) - initial_particle_number_density)
                * condition_.relaxation_coefficient_pnd * condition_.mass_density / (delta_time * delta_time * initial_particle_number_density);
    }
  }
  p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
  solveConjugateGradient(p_mat);
}

void BubbleParticles::correctVelocityDuan(const tiny_mps::Timer& timer) {
  using namespace tiny_mps;
  Grid grid(condition_.gradient_radius, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
  correction_velocity.setZero();
  int tensor_count = 0;
  int not_tensor_count = 0;
  for (int i_particle = 0; i_particle < size; ++i_particle) {
    if (particle_types(i_particle) != ParticleType::NORMAL) continue;
    if (boundary_types(i_particle) == BoundaryType::OTHERS) continue;
    Grid::Neighbors neighbors;
    grid.getNeighbors(i_particle, neighbors);
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
    Eigen::Vector3d tmp_vel(0.0, 0.0, 0.0);
    if (free_surface_type(i_particle) == SurfaceLayer::INNER_SURFACE) {
      for (int j_particle : neighbors) {
        if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
        Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
        tmp_vel += r_ij * (pressure(j_particle) + pressure(i_particle)) * weightForGradientPressure(r_ij) / r_ij.squaredNorm();
      }
      if (dimension == 2) tmp_vel(2) = 0;
      correction_velocity.col(i_particle) -= tmp_vel * dimension * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition_.mass_density);
    } else {
      double p_min = pressure(i_particle);
      double p_max = pressure(i_particle);
      for (int j_particle : neighbors) {
        if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
        p_min = std::min(pressure(j_particle), p_min);
        p_max = std::max(pressure(j_particle), p_max);
      }
      for (int j_particle : neighbors) {
        if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
        Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
        Eigen::Vector3d n_ij = r_ij.normalized();
        Eigen::Matrix3d tmp_tensor = Eigen::Matrix3d::Zero();
        tmp_tensor << n_ij(0) * n_ij(0), n_ij(0) * n_ij(1), n_ij(0) * n_ij(2),
                      n_ij(1) * n_ij(0), n_ij(1) * n_ij(1), n_ij(1) * n_ij(2),
                      n_ij(2) * n_ij(0), n_ij(2) * n_ij(1), n_ij(2) * n_ij(2);
        tensor += tmp_tensor * weightForGradientPressure(r_ij) / initial_particle_number_density;
        double xi = 0.2 + 2 * normal_vector.col(j_particle).norm();
        tmp_vel += r_ij * (pressure(j_particle) - pressure(i_particle) + xi * (p_max - p_min)) * weightForGradientPressure(r_ij) / r_ij.squaredNorm();
      }
      if (dimension == 2) {
        tmp_vel(2) = 0;
        tensor(2, 2) = 1.0;
      }
      if (tensor.determinant() > 1.0e-10) {
        correction_velocity.col(i_particle) -= tensor.inverse() * tmp_vel * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition_.mass_density);
        ++tensor_count;
      } else {
        correction_velocity.col(i_particle) -= tmp_vel * dimension * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition_.mass_density);
        ++not_tensor_count;
      }
    }
  }
  std::cout << "Tensor: " << tensor_count << ", Not Tensor: " << not_tensor_count << std::endl;
  temporary_velocity += correction_velocity;
}

}