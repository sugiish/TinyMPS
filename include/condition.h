// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_CONDITION_H_INCLUDED
#define MPS_CONDITION_H_INCLUDED

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>
#include <Eigen/Dense>

namespace tiny_mps {

// Holds analysis conditions.
class Condition {
public:
    Condition(std::string path) {
        readDataFile(path);
        getValue("average_distance",  average_distance);
        getValue("dimension", dimension);
        if(dimension != 2 && dimension != 3) {
            std::cerr << "Error: " << dimension << "-dimension is not supported." << std::endl;
            exit(EXIT_FAILURE);
        }

        double gx, gy, gz;
        getValue("gravity_x", gx);
        getValue("gravity_y", gy);
        getValue("gravity_z", gz);
        gravity(0) = gx;
        gravity(1) = gy;
        gravity(2) = gz;
        if (dimension == 2) gravity(2) = 0;

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

        getValue("initial_time", initial_time);
        getValue("delta_time", delta_time);
        getValue("finish_time", finish_time);
        getValue("min_delta_time", min_delta_time);
        getValue("output_interval", output_interval);
        getValue("inner_particle_index", inner_particle_index);
        getValue("relaxation_coefficient_lambda" ,relaxation_coefficient_lambda);
        getValue("surface_parameter", surface_parameter);
        getValue("tanaka_masunaga_method", tanaka_masunaga_method);
        getValue("tanaka_masunaga_gamma", tanaka_masunaga_gamma);
        getValue("tanaka_masunaga_c", tanaka_masunaga_c);
        getValue("tanaka_masunaga_beta", tanaka_masunaga_beta);

        pnd_weight_radius = pnd_influence * average_distance;
        gradient_radius = gradient_influence * average_distance;
        laplacian_pressure_weight_radius = laplacian_pressure_influence * average_distance;
        laplacian_viscosity_weight_radius = laplacian_viscosity_influence * average_distance;
    }

    virtual ~Condition(){}

    inline int getValue(const std::string& item, int& value) {
        if(data.find(item) == data.end()) return 1;
        std::stringstream ss;
        ss << data[item];
        ss >> value;
        return 0;
    }

    inline int getValue(const std::string& item, double& value) {
        if(data.find(item) == data.end()) return 1;
        std::stringstream ss;
        ss << data[item];
        ss >> value;
        return 0;
    }

    inline int getValue(const std::string& item, bool& value) {
        if(data.find(item) == data.end()) return 1;
        std::string tmp = data[item];
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
        if(tmp == "on" || tmp == "true") {
            value = true;
        } else {
            value = false;
        }
        return 0;
    }
    
    inline int getValue(const std::string& item, std::string& value) {
        if(data.find(item) == data.end()) return 1;
        value = data[item];
        return 0;
    }

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

    int inner_particle_index;
    double surface_parameter;
    double relaxation_coefficient_lambda;

    double pnd_weight_radius;
    double gradient_radius;
    double laplacian_pressure_weight_radius;
    double laplacian_viscosity_weight_radius;

    bool tanaka_masunaga_method;
    double tanaka_masunaga_gamma;
    double tanaka_masunaga_c;
    double tanaka_masunaga_beta;

private:
    inline int readDataFile(std::string path) {
        std::ifstream ifs(path);
        if(ifs.fail()) {
            std::cerr << "Error: in readDataFile() in condition.h" << std::endl;
            std::cerr << "Failed to read files: " << path << std::endl;
            exit(EXIT_FAILURE);
            return 1;
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
                if(first == '#') continue;  // Lines that begin with '#' are comments
            }
            item = std::regex_replace(item, re, "");
            item = std::regex_replace(item, re2, "");
            ss >> value;
            data[item] = value;
            std::cout << "    " << item << ": " << value << std::endl;
        }
        std::cout << std::endl;
        return 0;
    }
    
    std::unordered_map<std::string, std::string> data;
    double pnd_influence;
    double gradient_influence;
    double laplacian_pressure_influence;
    double laplacian_viscosity_influence;
};

}

#endif //MPS_CONDITION_H_INCLUDED