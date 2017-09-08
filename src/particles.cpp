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
	dimension = condition.dimension;
	VectorXb valid = particle_types.array() != ParticleType::GHOST;
	pnd_weight_radius = condition.pnd_influence * condition.average_distance;
	gradient_radius = condition.gradient_influence * condition.average_distance;
	laplacian_pressure_weight_radius = condition.laplacian_pressure_influence * condition.average_distance;
	laplacian_viscosity_weight_radius = condition.laplacian_viscosity_influence * condition.average_distance;
	pnd_weight = std::bind(&Particles::weightFunction, this, std::placeholders::_1, std::placeholders::_2, pnd_weight_radius);
	Grid pnd_grid(pnd_weight_radius, position, valid, condition.dimension);
	Grid lap_grid(laplacian_pressure_weight_radius, position, valid, condition.dimension);
	updateParticleNumberDensity(pnd_grid);
	setInitialParticleNumberDensity(condition.inner_particle_index);
	calculateLaplacianLambda(condition.inner_particle_index, lap_grid);
	checkSurfaceParticles(condition.surface_parameter);
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
	boundary_types = Eigen::VectorXi::Zero(size);
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
	ofs << std::endl;
	
	ofs << "SCALARS BoundaryCondition int" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	for(int i = 0; i < size; i++) {
		ofs << boundary_types(i) << std::endl;
	}
	return 0;
}

bool Particles::checkNeedlessCalculation() {
	for (int i_particle = 0; i_particle < size; ++i_particle) {
		if (particle_types(i_particle) == ParticleType::NORMAL) return false;
	}
	return true;
}

void Particles::updateParticleNumberDensity(Grid& grid) {
	grid.sumAllNeighborScalars(pnd_weight, particle_number_density);
	/*std::vector<int> neighbors(initial_neighbors_size);
	for (int i_particle = 0; i_particle < size; ++i_particle) {
		if (particle_types(i_particle) == ParticleType::GHOST) {
			particle_number_density(i_particle) = 0.0;
			continue;
		}
		grid.getNeighbors(i_particle, neighbors);
		double pnd = 0.0;
		for (int j_particle : neighbors) {
			if (particle_types(i_particle) == ParticleType::GHOST) continue;
			Eigen::Vector3d r_ji = particle.position.col(j_particle) - particle.position.col(i_particle);

		}
	}*/
}

void Particles::updateParticleNumberDensity(Grid& grid, std::function<double(int, int)> weight) {
	grid.sumAllNeighborScalars(weight, particle_number_density);
}

void Particles::setInitialParticleNumberDensity(int index) {
	initial_particle_number_density = particle_number_density(index);
}

void Particles::calculateLaplacianLambda(int index, Grid& grid) {
	auto lap_weight = std::bind(&Particles::weightFunction, this, std::placeholders::_1, std::placeholders::_2, laplacian_pressure_weight_radius);
	double sum_weight = grid.sumNeighborScalars(index, lap_weight);
	auto lap_weight_norm2 = std::bind(&Particles::laplacianWeightWithNorm2, this, std::placeholders::_1, std::placeholders::_2);
	double sum_weight_norm2 = grid.sumNeighborScalars(index, lap_weight_norm2);
	laplacian_lambda = sum_weight_norm2 / sum_weight;
}

void Particles::calculateTemporaryVelocity(const Eigen::Vector3d& force, Grid& grid, const Timer& timer, const Condition& condition) {
	double delta_time = timer.getCurrentDeltaTime();
	temporary_velocity = velocity;
	temporary_velocity.colwise() += delta_time * force;
	if (condition.viscosity_calculation) {
		Eigen::Matrix3Xd lap_vel;
		auto lap_vel_func = std::bind(&Particles::laplacianViscosity, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		grid.sumAllNeighborVectors(lap_vel_func, lap_vel);
		temporary_velocity += lap_vel * condition.kinematic_viscosity * delta_time;
	}
	for (int i = 0; i < size; i++) {
		if (particle_types(i) != ParticleType::NORMAL) temporary_velocity(i) = 0;
	}
	temporary_position = position + delta_time * temporary_velocity;
}

void Particles::moveExplicitly() {
	velocity = temporary_velocity;
	position = temporary_position;
}

void Particles::solvePressurePoission(Grid& grid, const Timer& timer, const Condition& condition) {
	using T = Eigen::Triplet<double>;
	double lap_r = condition.laplacian_pressure_influence * condition.average_distance;
	int n_size = (int)(std::pow(lap_r * 2, condition.dimension));
	double delta_time = timer.getCurrentDeltaTime();
	Eigen::SparseMatrix<double> p_mat(size, size);
	Eigen::VectorXd source(size);
	Eigen::VectorXd x(size);
	std::vector<T> coeffs(size * n_size);
	std::vector<int> neighbors(n_size * 2);
	for (int i_particle = 0; i_particle < size; i_particle++) {
		if (boundary_types(i_particle) == BoundaryType::OTHERS) {
			continue;
		} else if (boundary_types(i_particle) == BoundaryType::SURFACE) {
			coeffs.push_back(T(i_particle, i_particle, 1.0));
			continue;
		}
		grid.getNeighbors(i_particle, neighbors);
		double sum = 0.0;
		for (int j_particle : neighbors) {
			if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
			Eigen::Vector3d r_ji = position.col(j_particle) - position.col(i_particle);
			// double mat_ij = weightFunc(r_ji, lap_r) * 2 * condition.dimension
					// / (laplacian_lambda * initial_particle_number_density);
			double mat_ij = -weightFunc(r_ji, lap_r);
			sum -= mat_ij;
			if (boundary_types(j_particle) == BoundaryType::INNER) {
				coeffs.push_back(T(i_particle, j_particle, mat_ij));
			}
		}
		coeffs.push_back(T(i_particle, i_particle, sum));
<<<<<<< .merge_file_N55Z8T
		source(i_particle) = - (particle_number_density(i_particle) - initial_particle_number_density) * condition.mass_density
					/ (delta_time * delta_time * initial_particle_number_density);
		// std::cout << boundary_types(i_particle) << source(i_particle) << ", " << particle_number_density(i_particle) << ", " << initial_particle_number_density << std::endl;
=======
		// double tmp_source = - (particle_number_density(i_particle) - initial_particle_number_density) * condition.mass_density
		// / (delta_time * delta_time * initial_particle_number_density);
		double tmp_source = (particle_number_density(i_particle) - initial_particle_number_density) * condition.mass_density * laplacian_lambda
		/ (delta_time * delta_time * 2 * dimension);
		if (particle_number_density(i_particle) - initial_particle_number_density > 1.0e-9) source(i_particle) = tmp_source;
		// std::cout << i_particle << ": " <<  particle_number_density(i_particle) - initial_particle_number_density << ", " << source(i_particle) << std::endl;
>>>>>>> .merge_file_wXslES
	}
	p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
	
	// Solving a problem
<<<<<<< .merge_file_N55Z8T
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
=======
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower> cg;
	cg.setMaxIterations(60);
	cg.setTolerance(1.0e-9);
>>>>>>> .merge_file_wXslES
	cg.compute(p_mat);
	if (cg.info() != Eigen::ComputationInfo::Success) std::cout << "decompostion failed." << std::endl;
	x = cg.solve(source);
	if (cg.info() != Eigen::ComputationInfo::Success) std::cout << "solving failed." << std::endl;
	pressure = x;
	std::cout << "#iterations:     " << cg.iterations() << std::endl;
	std::cout << "estimated error: " << cg.error()      << std::endl;
	// solveConjugateGradient(p_mat, source, pressure, 2000, 10e-9);
	for (int i = 0; i < size; i++) {
		if (pressure(i) < 0) pressure(i) = 0;
	}
}

void Particles::solvePressurePoission2(Grid& grid, const Timer& timer, const Condition& condition) {
	using T = Eigen::Triplet<double>;
	double lap_r = condition.laplacian_pressure_influence * condition.average_distance;
	int mat_size = 0;
	std::vector<int> map(size);
	std::vector<int> toJ(size);
	for (int i_particle = 0; i_particle < size; i_particle++) {
		if (boundary_types(i_particle) == BoundaryType::OTHERS) {
			continue;
		}
		map.push_back(i_particle);
		toJ[i_particle] = mat_size++;
		mat_size++;
	}
	int n_size = (int)(std::pow(lap_r * 2, condition.dimension));
	double delta_time = timer.getCurrentDeltaTime();
	Eigen::SparseMatrix<double> p_mat(mat_size, mat_size);
	Eigen::VectorXd source(mat_size);
	std::vector<T> coeffs(size * n_size);
	std::vector<int> neighbors(n_size * 2);
	int count = 0;
	for (int i_particle = 0; i_particle < size; i_particle++) {
		if (boundary_types(i_particle) == BoundaryType::OTHERS) {
			coeffs.push_back(T(i_particle, i_particle, 1.0));
			continue;
		} else if (boundary_types(i_particle) == BoundaryType::SURFACE) {
			coeffs.push_back(T(i_particle, i_particle, 1.0));
			count++;
			continue;
		}
		grid.getNeighbors(i_particle, neighbors);
		double sum = 0.0;
		for (int j_particle : neighbors) {
			if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
			Eigen::Vector3d r_ji = position.col(j_particle) - position.col(i_particle);
			double mat_ij = weightFunc(r_ji, lap_r) * 2 * condition.dimension
				/ (laplacian_lambda * initial_particle_number_density);
			sum -= mat_ij;
			if (boundary_types(j_particle) == BoundaryType::INNER) {
				coeffs.push_back(T(toJ[i_particle], toJ[j_particle], mat_ij));
			}
		}
		coeffs.push_back(T(toJ[i_particle], toJ[i_particle], sum));
		source(toJ[i_particle]) = - (particle_number_density(i_particle) - initial_particle_number_density) * condition.mass_density
					/ (delta_time * delta_time * initial_particle_number_density);
		// std::cout << boundary_types(i_particle) << source(i_particle) << ", " << particle_number_density(i_particle) << ", " << initial_particle_number_density << std::endl;
		count++;
	}
	p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
	
	// Solving a problem
	// Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
	// cg.compute(p_mat);
	// pressure = cg.solve(source);
	// std::cout << "#iterations:     " << cg.iterations() << std::endl;
	// std::cout << "estimated error: " << cg.error()      << std::endl;
	Eigen::VectorXd x;
	solveConjugateGradient(p_mat, source, x, 2000, 10e-9);
	for (int i = 0; i < mat_size; i++) { 
		pressure(map[i]) = x(i);
	}
	for (int i = 0; i < size; i++) {
		if (pressure(i) < 0) pressure(i) = 0;
	}
}

/// solve Ax = b
void Particles::solveConjugateGradient (const Eigen::SparseMatrix<double> & A, const Eigen::VectorXd& b, Eigen::VectorXd& x, int itr, double eps) {
	x = Eigen::VectorXd::Zero(b.size());
	Eigen::VectorXd r;
	Eigen::VectorXd p;	
	r = b - A * x;
	p = r;
	int k;
	for (k = 0; k < itr; ++k) {
		Eigen::VectorXd y = A * p;
		double alpha = r.squaredNorm() / p.dot(y);
		x += alpha * p;
		Eigen::VectorXd r_next = r - alpha * y;
		if (r.norm() < eps) {	
			break;
		}
		double beta = r_next.squaredNorm() / r.squaredNorm();
		p *= beta;
		p += r_next;
	}
	std::cout << "it: " << k << std::endl;
	std::cout << "err: " << r.norm() << std::endl;
}

void Particles::advectVelocity(Grid& grid, const Timer& timer, const Condition& condition) {
	// int n_size = (int)(std::pow(condition.gradient_influence * 2, condition.dimension));
	std::vector<int> neighbors(30);
	for (int i_particle = 0; i_particle < size; i_particle++) {
		if (particle_types(i_particle) != ParticleType::NORMAL) continue;
		grid.getNeighbors(i_particle, neighbors);
		double p_min = pressure(i_particle);
		for (int j_particle : neighbors) {
			if (particle_types(j_particle) != ParticleType::NORMAL) continue;
			p_min = std::min(pressure(j_particle), p_min);
		}
		Eigen::Vector3d adv_v;
		for (int j_particle : neighbors) {
			if (particle_types(j_particle) != ParticleType::NORMAL) continue;
			Eigen::Vector3d r_ji = position.col(j_particle) - position.col(i_particle);
			adv_v += r_ji * (pressure(j_particle) - p_min) * weightFunc(r_ji, condition.gradient_influence * condition.average_distance) / r_ji.squaredNorm();
		}
		if (condition.dimension == 2) adv_v(2) = 0;
		temporary_velocity.col(i_particle) += -adv_v * condition.dimension * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition.kinematic_viscosity);
		std::cout << adv_v.norm() << std::endl;
		temporary_position.col(i_particle) += adv_v * timer.getCurrentDeltaTime();
	}
	velocity = temporary_velocity;
	position = temporary_position;
}

void Particles::checkSurfaceParticles(double surface_parameter) {
	for(int i = 0; i < getSize(); i++) {
		if(particle_types(i) == ParticleType::NORMAL || particle_types(i) == ParticleType::WALL) {
			if (particle_number_density(i) < surface_parameter * initial_particle_number_density) boundary_types(i) = BoundaryType::SURFACE;
			else boundary_types(i) = BoundaryType::INNER;
		} else {
			boundary_types(i) = BoundaryType::OTHERS;
		}
	}
}

double Particles::weightFunc(Eigen::Vector3d& vec, double influence_radius) {
	double r = vec.norm();
	if(r < influence_radius) return (influence_radius / r - 1.0);
	else return 0.0;
}

double Particles::weightFunction(int i_particle, int j_particle, double influence_radius) {
	Eigen::Vector3d v = position.col(j_particle) - position.col(i_particle);
	double r = v.norm();
	if(r < influence_radius) return (influence_radius / r - 1.0);
	else return 0.0;
}

double Particles::laplacianWeightWithNorm2(int i_particle, int j_particle) {
	Eigen::Vector3d v = position.col(j_particle) - position.col(i_particle);
	double r = v.norm();
	if(r < laplacian_pressure_weight_radius) return r * r * (laplacian_pressure_weight_radius / r - 1.0);
	else return 0.0;
}

void Particles::laplacianViscosity(int i_particle, int j_particle, Eigen::Vector3d& output) {
	Eigen::Vector3d vel = velocity.col(j_particle) - velocity.col(i_particle);
	double w = weightFunction(i_particle, j_particle, laplacian_pressure_weight_radius);
	output = vel * (w * 2 * dimension / (laplacian_lambda * initial_particle_number_density));
}

} // namespace tiny_mps