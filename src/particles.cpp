#include "particles.h"

#include <fstream>
#include <iostream>
#include <sstream>

namespace tiny_mps {

/**
 * Class for describing particles status.
 */
Particles::Particles(const std::string& path, const Condition& condition) {
	readGridFile(path, condition.dimension);
	VectorXb valid = particle_types.array() != ParticleType::GHOST;
	pnd_weight_radius = condition.pnd_influence * condition.average_distance;
	pnd_weight = std::bind(&Particles::weightFunction, this, std::placeholders::_1, std::placeholders::_2, pnd_weight_radius);
	laplacian_pressure_weight_radius = condition.laplacian_pressure_influence * condition.average_distance;
	Grid pnd_grid(pnd_weight_radius, position, valid, condition.dimension);
	Grid lap_grid(laplacian_pressure_weight_radius, position, valid, condition.dimension);
	updateParticleNumberDensity(pnd_grid);
	calculateInitialParticleNumberDensity(condition.initial_pnd_index);
	calculateLaplacianLambda(condition.initial_pnd_index, lap_grid);
}

Particles::~Particles() {}

void Particles::initialize(int size) {
	this->size = size;
	position = Eigen::MatrixXd::Zero(3, size);
	velocity = Eigen::MatrixXd::Zero(3, size);
	pressure = Eigen::VectorXd::Zero(size);
	particle_number_density = Eigen::VectorXd::Zero(size);
	temporary_position = Eigen::MatrixXd::Zero(3, size);
	temporary_velocity = Eigen::MatrixXd::Zero(3, size);
	particle_types = Eigen::VectorXi::Zero(size);
}

int Particles::readGridFile(const std::string& path, int dimension) {
	std::ifstream ifs(path);
	if (ifs.fail()) {
		std::cerr << "Error: in readGridFile()" << std::endl;
		std::cerr << "Failed to read files: " << path << std::endl;
		return 1;
	}

	std::string tmp_str;
	//Line 0: Start time
	getline(ifs, tmp_str);
	//Line 1: particles_number
	getline(ifs, tmp_str);

	{
		std::stringstream ss;
		int ptcl_num = 0;

		ss.str(tmp_str);
		ss >> ptcl_num;
		initialize(ptcl_num);
	}

	int i_counter = 0;
	while(getline(ifs, tmp_str)) {
		std::stringstream ss;
		ss.str(tmp_str);
		ss >> particle_types(i_counter);

		for(int i_dim = 0; i_dim < 3; i_dim++) {
			double tmp;
			ss >> tmp;
			if (i_dim < dimension) position(i_dim, i_counter) = tmp;
		}
		for(int i_dim = 0; i_dim < 3; i_dim++) {
			double tmp;
			ss >> tmp;
			if (i_dim < dimension) velocity(i_dim, i_counter) = tmp;
		}
		ss >> pressure(i_counter);
		i_counter++;
	}
	return 0;
}

int Particles::writeVtkFile(const std::string& path, const std::string& title) {
	std::ofstream ofs(path);
	if(ofs.fail()) {
		std::cerr << "Error: in writeVtkFile()" << std::endl;
		return 1;
	}

	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << title << std::endl;
	ofs << "ASCII" << std::endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	ofs << std::endl;

	ofs << "POINTS " << size << " double" << std::endl;
	for(int i = 0; i < size; i++) {
		ofs << position(0, i) << " " << position(1, i) << " " << position(2, i) << std::endl;
	}
	ofs << std::endl;

	ofs << "CELL_TYPES " << size << std::endl;
	for(int i = 0; i < size; i++) {
		ofs << 1 << std::endl;
	}
	ofs << std::endl;

	ofs << "POINT_DATA " << size << std::endl;
	ofs << "SCALARS Pressure double" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	for(int i = 0; i < size; i++) {
		ofs << pressure(i) << std::endl;
	}
	ofs << std::endl;

	ofs << "VECTORS Velocity double" << std::endl;
	for(int i = 0; i < size; i++) {
		ofs << velocity(0, i) << " " << velocity(1, i) << " " << velocity(2, i) << std::endl;
	}
	ofs << std::endl;

	ofs << "SCALARS Type int" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	for(int i = 0; i < size; i++) {
		ofs << particle_types(i) << std::endl;
	}
	ofs << std::endl;

	ofs << "SCALARS ParticleNumberDensity double" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	for(int i = 0; i < size; i++) {
		ofs << particle_number_density(i) << std::endl;
	}
	
	return 0;
}

void Particles::updateParticleNumberDensity(Grid& grid) {
	grid.sumAllNeighborScalars(pnd_weight, particle_number_density);
}

void Particles::updateParticleNumberDensity(Grid& grid, std::function<double(int, int)> weight) {
	grid.sumAllNeighborScalars(weight, particle_number_density);
}

void Particles::calculateInitialParticleNumberDensity(int index) {
	initial_particle_number_density = particle_number_density(index);
}

void Particles::calculateLaplacianLambda(int index, Grid& grid) {
	
}

void Particles::moveParticlesExplicitly(const Eigen::Vector3d& force, Timer timer) {
	double delta_time = timer.getCurrentDeltaTime();
	temporary_velocity = velocity;
	temporary_velocity.colwise() += delta_time * force;
	temporary_position = position + delta_time * temporary_velocity;
}

double Particles::weightFunction(int i_particle, int j_particle, double influence_radius) {
	Eigen::Vector3d v = position.col(j_particle) - position.col(i_particle);
	double r = v.norm();
	if(r < influence_radius) return (influence_radius / r - 1.0);
	else return 0.0;
}

} // namespace tiny_mps