// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include "condition.h"

namespace tiny_mps{

Condition::Condition(std::string path) {
  readDataFile(path);
  getValue("average_distance",  average_distance);
  getValue("dimension", dimension);
  if(dimension != 2 && dimension != 3) {
    std::cerr << "Error: " << dimension << "-dimension is not supported." << std::endl;
    throw std::out_of_range("Error: dimension is out of range.");
  }
  double gx, gy, gz;
  getValue("gravity_x", gx);
  getValue("gravity_y", gy);
  getValue("gravity_z", gz);
  gravity(0) = gx;
  gravity(1) = gy;
  gravity(2) = (dimension == 3)? gz : 0;
  double ix, iy, iz;
  getValue("inflow_x", ix);
  getValue("inflow_y", iy);
  getValue("inflow_z", iz);
  inflow_velocity(0) = ix;
  inflow_velocity(1) = iy;
  inflow_velocity(2) = (dimension == 3)? iz : 0;

  getValue("temperature", temperature);
  getValue("head_pressure", head_pressure);
  getValue("mass_density", mass_density);
  getValue("viscosity_calculation", viscosity_calculation);
  getValue("kinematic_viscosity", kinematic_viscosity);

  getValue("courant_number", courant_number);
  getValue("diffusion_number", diffusion_number);

  getValue("pnd_influence", pnd_influence);
  getValue("gradient_influence", gradient_influence);
  getValue("laplacian_pressure_influence", laplacian_pressure_influence);
  getValue("laplacian_viscosity_influence", laplacian_viscosity_influence);
  getValue("additional_ghost_particles", additional_ghost_particles);
  getValue("extra_ghost_particles", extra_ghost_particles);
  getValue("collision_influence", collision_influence);
  getValue("restitution_coefficent", restitution_coefficent);

  getValue("initial_time", initial_time);
  getValue("delta_time", delta_time);
  getValue("finish_time", finish_time);
  getValue("min_delta_time", min_delta_time);
  getValue("output_interval", output_interval);
  getValue("relaxation_coefficient_pnd", relaxation_coefficient_pnd);
  getValue("relaxation_coefficient_vel_div", relaxation_coefficient_vel_div);
  getValue("weak_compressibility", weak_compressibility);
  getValue("surface_threshold_pnd", surface_threshold_pnd);
  getValue("surface_threshold_number", surface_threshold_number);

  getValue("initial_void_fraction", initial_void_fraction);
  getValue("min_void_fraction", min_void_fraction);
  getValue("bubble_density", bubble_density);
  getValue("vapor_pressure", vapor_pressure);

  pnd_weight_radius = pnd_influence * average_distance;
  gradient_radius = gradient_influence * average_distance;
  laplacian_pressure_weight_radius = laplacian_pressure_influence * average_distance;
  laplacian_viscosity_weight_radius = laplacian_viscosity_influence * average_distance;
}

void Condition::readDataFile(std::string path) {
  std::ifstream ifs(path);
  if(ifs.fail()) {
    std::cerr << "Error: in readDataFile() in condition.h" << std::endl;
    std::cerr << "Failed to read files: " << path << std::endl;
    throw std::ios_base::failure("Error: in readDataFile() in condition.h");
  }
  std::cout << "Succeed in reading data file: " << path << std::endl;
  std::string tmp_str;
  std::regex re("\\(.*\\)");          // For removing (**)
  std::regex re2("-+\\w+-+");         // For removing like --**--
  while(getline(ifs, tmp_str)) {
    if(tmp_str.empty()) continue;
    std::stringstream ss;
    ss.str(tmp_str);
    std::string item, value;
    ss >> item;
    {
      char first = item.at(0);
      if(first == '#') continue;      // Lines that begin with '#' are comments
    }
    item = std::regex_replace(item, re, "");
    item = std::regex_replace(item, re2, "");
    ss >> value;
    data[item] = value;
    std::cout << "    " << item << ": " << value << std::endl;
  }
  std::cout << std::endl;
}

int Condition::getValue(const std::string& item, int& value) const {
  if(data.find(item) == data.end()) return 1;
  std::stringstream ss;
  ss << data.at(item);
  ss >> value;
  return 0;
}

int Condition::getValue(const std::string& item, double& value) const {
  if(data.find(item) == data.end()) return 1;
  std::stringstream ss;
  ss << data.at(item);
  ss >> value;
  return 0;
}

int Condition::getValue(const std::string& item, bool& value) const {
  if(data.find(item) == data.end()) return 1;
  std::string tmp = data.at(item);
  std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
  if(tmp == "on" || tmp == "true") {
    value = true;
  } else {
    value = false;
  }
  return 0;
}

int Condition::getValue(const std::string& item, std::string& value) const {
  if(data.find(item) == data.end()) return 1;
  value = data.at(item);
  return 0;
}

} // namespace tiny_mps