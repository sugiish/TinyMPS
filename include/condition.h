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
        getValue("inner_particle_index", inner_particle_index);
        getValue("surface_parameter", surface_parameter);
    }

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

    double pnd_influence;
    double gradient_influence;
    double laplacian_pressure_influence;
    double laplacian_viscosity_influence;

    double initial_time;
    double finish_time;
    double delta_time;

    int inner_particle_index;
    double surface_parameter;

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

private:
    std::unordered_map<std::string, std::string> data;

    inline int readDataFile(std::string path) {
        std::ifstream ifs(path);
        if(ifs.fail()) {
            std::cerr << "Error: in Reader()" << std::endl;
            std::cerr << "Failed to read files: " << path << std::endl;
            return 1;
        }

        std::string tmp_str;
        std::regex re("\\(.*\\)"); // For removing (**)
        std::regex re2("-+\\w+-+");// For removing like --**--
        while(getline(ifs, tmp_str)) {
            if(tmp_str.empty()) continue;
            std::stringstream ss;
            ss.str(tmp_str);

            std::string item, value;
            ss >> item;
            {
                // Lines that begin with '#' are comments
                char first = item.at(0);
                if(first == '#') continue;
            }

            item = std::regex_replace(item, re, "");
            item = std::regex_replace(item, re2, "");
            
            ss >> value;
            data[item] = value;
        }
        return 0;
    }
};

}

#endif //MPS_CONDITION_H_INCLUDED