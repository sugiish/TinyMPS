// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_CONDITION_H_INCLUDED
#define MPS_CONDITION_H_INCLUDED

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>
#include <Eigen/Core>

namespace tiny_mps {

// Holds analysis conditions.
class Condition {
 public:
  explicit Condition(std::string path);
  // Condition is neither copyable nor movable.
  Condition(const Condition&) = delete;
  Condition& operator=(const Condition&) = delete;
  virtual ~Condition(){}

  double average_distance;
  int dimension;
  Eigen::Vector3d gravity;
  double mass_density;
  double temperature;
  double head_pressure;
  bool viscosity_calculation;
  double kinematic_viscosity;

  double courant_number;
  double diffusion_number;

  double initial_time;
  double finish_time;
  double delta_time;
  double min_delta_time;
  double output_interval;

  double surface_threshold_pnd;
  double surface_threshold_number;

  double weak_compressibility;

  double relaxation_coefficient_pnd;
  double relaxation_coefficient_vel_div;

  int extra_ghost_particles;
  int additional_ghost_particles;
  Eigen::Vector3d inflow_velocity;

  double collision_influence;
  double restitution_coefficent;

  double pnd_influence;
  double gradient_influence;
  double laplacian_pressure_influence;
  double laplacian_viscosity_influence;

  double pnd_weight_radius;
  double gradient_radius;
  double laplacian_pressure_weight_radius;
  double laplacian_viscosity_weight_radius;

  double initial_void_fraction;
  double min_void_fraction;
  double bubble_density;
  double vapor_pressure;
  double secondary_surface_eta;

 private:
  void readDataFile(std::string path);

  int getValue(const std::string& item, int& value) const;
  int getValue(const std::string& item, double& value) const;
  int getValue(const std::string& item, bool& value) const;
  int getValue(const std::string& item, std::string& value) const;

  std::unordered_map<std::string, std::string> data;
};

}

#endif //MPS_CONDITION_H_INCLUDED