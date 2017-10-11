// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include <iostream>
#include <boost/format.hpp>
#include <Eigen/Core>
#include "condition.h"
#include "grid.h"
#include "particles.h"

// Sample code using TinyMPS library.
int main(int argc, char* argv[]) {
  try {
    std::string input_data = "./input/input_tensor_verification.data";
    if (argc >= 2) input_data = argv[2];
    tiny_mps::Condition condition(input_data);
    tiny_mps::Particles regular(10000, condition);
    tiny_mps::Timer timer(condition);
    for (int y = 0; y < 100; ++y) {
      for (int x = 0; x < 100; ++x) {
        regular.position.col(100 * y + x) = Eigen::Vector3d(condition.average_distance * x, condition.average_distance * y, 0);
        regular.temporary_position = regular.position;
        regular.pressure(100 * y + x) = condition.mass_density * condition.gravity.norm() * regular.position(1, 100 * y + x);
      }
    }
    regular.correctVelocity(timer);
    regular.writeVtkFile("./output/regular_standard.vtk", "regular_standard");
    {
      std::ofstream ofs("./output/regular_standard_err.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(regular.correction_velocity(1, i) - 9.8);
        ofs << regular.position(0, i) << ", " << regular.position(1, i) << ", " << err << std::endl;
      }
    }
    regular.correctVelocityWithTensor(timer);
    regular.writeVtkFile("./output/regular_tensor.vtk", "regular_tensor");
    {
      std::ofstream ofs("./output/regular_tensor_err.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(regular.correction_velocity(1, i) - 9.8);
        ofs << regular.position(0, i) << ", " << regular.position(1, i) << ", " << err << std::endl;
      }
    }

    tiny_mps::Particles irregular = regular;
    for (int y = 0; y < 100; ++y) {
      for (int x = 0; x < 100; ++x) {
        Eigen::Vector3d rnd = Eigen::Vector3d::Random() * condition.average_distance * 0.1;
        rnd(2) = 0;
        irregular.position.col(100 * y + x) += rnd;
        irregular.temporary_position = irregular.position;
        irregular.pressure(100 * y + x) = condition.mass_density * condition.gravity.norm() * irregular.position(1, 100 * y + x);
      }
    }
    irregular.correctVelocity(timer);
    irregular.writeVtkFile("./output/irregular_standard.vtk", "irregular_standard");
    {
      std::ofstream ofs("./output/irregular_standard_err.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(irregular.correction_velocity(1, i) - 9.8);
        ofs << irregular.position(0, i) << ", " << irregular.position(1, i) << ", " << err << std::endl;
      }
    }
    irregular.correctVelocityWithTensor(timer);
    irregular.writeVtkFile("./output/irregular_tensor.vtk", "irregular_tensor");
    {
      std::ofstream ofs("./output/irregular_tensor_err.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(irregular.correction_velocity(1, i) - 9.8);
        ofs << irregular.position(0, i) << ", " << irregular.position(1, i) << ", " << err << std::endl;
      }
    }

    for (int y = 0; y < 100; ++y) {
      for (int x = 0; x < 100; ++x) {
        regular.pressure(100 * y + x) = cos(sqrt(regular.position(0, 100 * y + x) * regular.position(0, 100 * y + x) +
                                                 regular.position(1, 100 * y + x) * regular.position(1, 100 * y + x)));
      }
    }
    for (int i = 0; i < 10000; ++i) {
      regular.pressure(i) = cos(sqrt(regular.position(0, i) * regular.position(0, i) +
      regular.position(1, i) * regular.position(1, i)));
    }
    regular.correctVelocity(timer);
    regular.writeVtkFile("./output/regular_standard2.vtk", "regular_standard2");
    {
      std::ofstream ofs("./output/regular_standard_err2.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(regular.correction_velocity(1, i) - 9.8);
        ofs << regular.position(0, i) << ", " << regular.position(1, i) << ", " << err << std::endl;
      }
    }
    regular.correctVelocityWithTensor(timer);
    regular.writeVtkFile("./output/regular_tensor2.vtk", "regular_tensor2");
    {
      std::ofstream ofs("./output/regular_tensor_err2.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(regular.correction_velocity(1, i) + sin(sqrt(regular.position(0, i) * regular.position(0, i) +
        regular.position(1, i) * regular.position(1, i))) * regular.position(1, i) / sqrt(regular.position(0, i) * regular.position(0, i) +
        regular.position(1, i) * regular.position(1, i)));
        ofs << regular.position(0, i) << ", " << regular.position(1, i) << ", " << err << std::endl;
      }
    }

    for (int i = 0; i < 10000; ++i) {
      irregular.pressure(i) = cos(sqrt(irregular.position(0, i) * irregular.position(0, i) +
                                                  irregular.position(1, i) * irregular.position(1, i)));
    }
    irregular.correctVelocity(timer);
    irregular.writeVtkFile("./output/irregular_standard2.vtk", "irregular_standard2");
    {
      std::ofstream ofs("./output/irregular_standard_err2.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(irregular.correction_velocity(1, i) + sin(sqrt(irregular.position(0, i) * irregular.position(0, i) +
        irregular.position(1, i) * irregular.position(1, i))) * irregular.position(1, i) / sqrt(irregular.position(0, i) * irregular.position(0, i) +
        irregular.position(1, i) * irregular.position(1, i)));
        ofs << irregular.position(0, i) << ", " << irregular.position(1, i) << ", " << err << std::endl;
      }
    }
    irregular.correctVelocityWithTensor(timer);
    irregular.writeVtkFile("./output/irregular_tensor2.vtk", "irregular_tensor2");
    {
      std::ofstream ofs("./output/irregular_tensor_err2.csv");
      for (int i = 0; i < 10000; ++i) {
        double err = abs(irregular.correction_velocity(1, i) - 9.8);
        ofs << irregular.position(0, i) << ", " << irregular.position(1, i) << ", " << err << std::endl;
      }
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }
  return 1;
}