#ifndef MPS_CONDITION_H_INCLUDED
#define MPS_CONDITION_H_INCLUDED

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>

#include <Eigen/Dense>

using namespace Eigen;

namespace tiny_mps
{

class Condition
{
public:
	Condition(std::string path)
	{
		readDataFile(path);
		getValue("average_distance",  average_distance);
		getValue("dimension", dimension);
		if(dimension != 2 && dimension != 3)
		{
			std::cerr << "Error: " << dimension << "-dimension is not supported." << std::endl;
		}

		double gx, gy, gz;
		getValue("gravity_x", gx);
		getValue("gravity_y", gy);
		getValue("gravity_z", gz);
		gravity(0) = gx;
		gravity(1) = gy;
		if(dimension == 3) gravity(2) = gz;
		else gravity(2) = gz;

		getValue("temperature", temperature);
		getValue("head_pressure", head_pressure);

		getValue("viscosity_calculation", viscosity_calculation);
		getValue("kinematic_viscosity", kinematic_viscosity);

		getValue("courant_number", courant_number);
		getValue("diffusion_number", diffusion_number);

		getValue("pnd_influence", pnd_influence);
		getValue("gradient_influence", gradient_influence);
		getValue("laplacian_influence", laplacian_influence);

		getValue("initial_time", initial_time);
		getValue("delta_time", delta_time);
		getValue("finish_time", finish_time);
		getValue("initial_pnd_index", initial_pnd_index);
	}

	virtual ~Condition(){}

	double average_distance;
	int dimension;

	Vector3d gravity;
	double temperature;
	double head_pressure;

	bool viscosity_calculation;
	double kinematic_viscosity;

	double courant_number;
	double diffusion_number;

	double pnd_influence;
	double gradient_influence;
	double laplacian_influence;

	double initial_time;
	double finish_time;
	double delta_time;

	int initial_pnd_index;

	inline int getValue(const std::string& item, int& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		std::stringstream ss;
		ss << data[item];
		ss >> value;
		return 0;
	}

	inline int getValue(const std::string& item, double& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		std::stringstream ss;
		ss << data[item];
		ss >> value;
		return 0;
	}

	inline int getValue(const std::string& item, bool& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		std::string tmp = data[item];
		std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
		if(tmp == "on" || tmp == "true")
		{
			value = true;
		}
		else
		{
			value = false;
		}
		return 0;
	}
	
	inline int getValue(const std::string& item, std::string& value)
	{
		if(data.find(item) == data.end())
		{
			return 1;
		}
		value = data[item];
		return 0;
	}
private:
	std::unordered_map<std::string, std::string> data;

	inline int readDataFile(std::string path)
	{
		std::ifstream ifs(path);
		
		if(ifs.fail())
		{
			std::cerr << "Error: in Reader()" << std::endl;
			std::cerr << "Failed to read files: " << path << std::endl;
			return 1;
		}

		std::string tmp_str;
		std::regex re("\\(.*\\)"); // For removing (**)
		std::regex re2("-+\\w+-+");// For removing like --**--
		while(getline(ifs, tmp_str))
		{
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