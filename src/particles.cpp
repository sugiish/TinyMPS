// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include "particles.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>

namespace tiny_mps {

Particles::Particles(const std::string& path, const Condition& condition) : condition_(condition) {
    readGridFile(path, condition.dimension);
    dimension = condition.dimension;
    VectorXb valid = particle_types.array() != ParticleType::GHOST;
    Grid lap_grid(condition.laplacian_pressure_weight_radius, position, valid, condition.dimension);
    updateParticleNumberDensity(condition);
    setInitialParticleNumberDensity(condition.inner_particle_index);
    calculateLaplacianLambda(condition.inner_particle_index, lap_grid);
    checkSurfaceParticles(condition.surface_parameter);
    std::cout << "Laplacian lambda:" << laplacian_lambda << std::endl;
    std::cout << "initial_pnd" << initial_particle_number_density << std::endl;
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
    correction_velocity = Eigen::MatrixXd::Zero(3, size);
}

int Particles::readGridFile(const std::string& path, int dimension) {
    std::ifstream ifs(path);
    if (ifs.fail()) {
        std::cerr << "Error: in readGridFile() in particles.cpp" << std::endl;
        std::cerr << "Failed to read files: " << path << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    std::string tmp_str;
    getline(ifs, tmp_str);      //Line 0: Start time
    getline(ifs, tmp_str);      //Line 1: particles_number
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
        for(int i_dim = 0; i_dim < 3; ++i_dim) {
            double tmp;
            ss >> tmp;
            if (i_dim < dimension) position(i_dim, i_counter) = tmp;
        }
        for(int i_dim = 0; i_dim < 3; ++i_dim) {
            double tmp;
            ss >> tmp;
            if (i_dim < dimension) velocity(i_dim, i_counter) = tmp;
        }
        ss >> pressure(i_counter);
        ++i_counter;
    }
    return 0;
}

int Particles::writeVtkFile(const std::string& path, const std::string& title) {
    std::ofstream ofs(path);
    if(ofs.fail()) {
        std::cerr << "Error: in writeVtkFile()" << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    ofs << "# vtk DataFile Version 2.0" << std::endl;
    ofs << title << std::endl;
    ofs << "ASCII" << std::endl;
    ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ofs << std::endl;
    ofs << "POINTS " << size << " double" << std::endl;
    for(int i = 0; i < size; ++i) {
        ofs << position(0, i) << " " << position(1, i) << " " << position(2, i) << std::endl;
    }
    ofs << std::endl;
    ofs << "CELL_TYPES " << size << std::endl;
    for(int i = 0; i < size; ++i) {
        ofs << 1 << std::endl;
    }
    ofs << std::endl;
    ofs << "POINT_DATA " << size << std::endl;
    ofs << "SCALARS Pressure double" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
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
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int i = 0; i < size; ++i) {
        ofs << particle_types(i) << std::endl;
    }
    ofs << std::endl;
    ofs << "SCALARS ParticleNumberDensity double" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int i = 0; i < size; ++i) {
        ofs << particle_number_density(i) << std::endl;
    }
    ofs << std::endl;
    ofs << "SCALARS BoundaryCondition int" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int i = 0; i < size; ++i) {
        ofs << boundary_types(i) << std::endl;
    }
    ofs << std::endl;
    ofs << "VECTORS CorrectionVelocity double" << std::endl;
    for(int i = 0; i < size; ++i) {
        ofs << correction_velocity(0, i) << " " << correction_velocity(1, i) << " " << correction_velocity(2, i) << std::endl;
    }
    return 0;
}

void Particles::saveInterval(const std::string& path, const Timer& timer) {
    if (!timer.isOutputTime()) return;
    writeVtkFile((boost::format(path) % timer.getOutputCount()).str(), (boost::format("Time: %s") % timer.getCurrentTime()).str());
}

bool Particles::checkNeedlessCalculation() {
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::NORMAL) return false;
    }
    return true;
}

void Particles::calculateTemporaryParticleNumberDensity(const Condition& condition) {
    Grid grid(condition.pnd_weight_radius, temporary_position, particle_types.array() != ParticleType::GHOST, dimension);
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::GHOST) {
            particle_number_density(i_particle) = 0.0;
            continue;
        }
        Grid::Neighbors neighbors;
        grid.getNeighbors(i_particle, neighbors);
        double pnd = 0.0;
        for (int j_particle : neighbors) {
            if (particle_types(i_particle) == ParticleType::GHOST) continue;
            Eigen::Vector3d r_ji = temporary_position.col(j_particle) - temporary_position.col(i_particle);
            pnd += weightFunction(r_ji, grid.getGridWidth());
        }
        particle_number_density(i_particle) = pnd;
    }
}

void Particles::updateParticleNumberDensity(const Condition& condition) {
    Grid grid(condition.pnd_weight_radius, position, particle_types.array() != ParticleType::GHOST, dimension);
    updateParticleNumberDensity(grid);
}

void Particles::updateParticleNumberDensity(const Grid& grid) {
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::GHOST) {
            particle_number_density(i_particle) = 0.0;
            continue;
        }
        Grid::Neighbors neighbors;
        grid.getNeighbors(i_particle, neighbors);
        double pnd = 0.0;
        for (int j_particle : neighbors) {
            if (particle_types(i_particle) == ParticleType::GHOST) continue;
            Eigen::Vector3d r_ji = position.col(j_particle) - position.col(i_particle);
            pnd += weightFunction(r_ji, grid.getGridWidth());
        }
        particle_number_density(i_particle) = pnd;
    }
}

void Particles::setInitialParticleNumberDensity(int index) {
    initial_particle_number_density = particle_number_density(index);
}

void Particles::calculateLaplacianLambda(int index, const Grid& grid) {
    double numerator = 0.0;
    double denominator = 0.0;
    Grid::Neighbors neighbors;
    grid.getNeighbors(index, neighbors);
    for (int j_particle : neighbors) {
        Eigen::Vector3d r_ji = position.col(j_particle) - position.col(index);
        double w =  weightFunction(r_ji, grid.getGridWidth());
        numerator += r_ji.squaredNorm() * w;
        denominator += w;
    }
    laplacian_lambda = numerator / denominator;
}

void Particles::calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer) {
    Grid grid(condition_.laplacian_viscosity_weight_radius, position, particle_types.array() == ParticleType::NORMAL, condition_.dimension);
    calculateTemporaryVelocity(force, grid, timer, condition_);
}

void Particles::calculateTemporaryVelocity(const Eigen::Vector3d& force, Grid& grid, const Timer& timer, const Condition& condition) {
    double delta_time = timer.getCurrentDeltaTime();
    temporary_velocity = velocity;
    temporary_velocity.colwise() += delta_time * force;
    if (condition.viscosity_calculation) {
        for (int i_particle = 0; i_particle < size; ++i_particle) {
            Grid::Neighbors neighbors;
            if (particle_types(i_particle) != ParticleType::NORMAL) continue;
            grid.getNeighbors(i_particle, neighbors);
            Eigen::Vector3d lap_vec(0.0, 0.0, 0.0);
            for (int j_particle : neighbors) {
                if (particle_types(j_particle) != ParticleType::NORMAL) continue;
                Eigen::Vector3d u_ji = velocity.col(j_particle) - velocity.col(i_particle);
                Eigen::Vector3d r_ji = position.col(j_particle) - position.col(i_particle);
                lap_vec += u_ji * weightFunction(r_ji, grid.getGridWidth()) * 2 * dimension / (laplacian_lambda * initial_particle_number_density);
            }
            temporary_velocity.col(i_particle) += lap_vec * condition.kinematic_viscosity * delta_time;
        }
    }
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) != ParticleType::NORMAL) temporary_velocity.col(i_particle).setZero();
    }
    temporary_position = position + delta_time * temporary_velocity;
}

void Particles::solvePressurePoission(const Timer& timer) {
    Grid grid(condition_.laplacian_pressure_weight_radius, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    solvePressurePoission(grid, timer, condition_);
}

void Particles::solvePressurePoission(const Grid& grid, const Timer& timer, const Condition& condition) {
    using T = Eigen::Triplet<double>;
    double lap_r = grid.getGridWidth();
    int n_size = (int)(std::pow(lap_r * 2, dimension));
    double delta_time = timer.getCurrentDeltaTime();
    Eigen::SparseMatrix<double> p_mat(size, size);
    Eigen::VectorXd source(size);
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
        for (int j_particle : neighbors) {
            if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
            Eigen::Vector3d r_ji = temporary_position.col(j_particle) - temporary_position.col(i_particle);
            double mat_ij = weightFunction(r_ji, lap_r) * 2 * dimension
                    / (laplacian_lambda * initial_particle_number_density);
            sum -= mat_ij;
            if (boundary_types(j_particle) == BoundaryType::INNER) {
                coeffs.push_back(T(i_particle, j_particle, mat_ij));
            }
        }
        coeffs.push_back(T(i_particle, i_particle, sum));
        source(i_particle) = - (particle_number_density(i_particle) - initial_particle_number_density) * condition.mass_density
                    / (delta_time * delta_time * initial_particle_number_density);
    }
    p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
    
    // Solving a problem
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.compute(p_mat);
    if (cg.info() != Eigen::ComputationInfo::Success) std::cerr << "decompostion failed." << std::endl;
    pressure = cg.solve(source);
    if (cg.info() != Eigen::ComputationInfo::Success) std::cerr << "solving failed." << std::endl;
    std::cout << "#iterations:     " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    for (int i = 0; i < size; ++i) {
        if (pressure(i) < 0) pressure(i) = 0;
    }
}

void Particles::correctVelocity(const Timer& timer) {
    Grid grid(condition_.gradient_radius, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    correctVelocity(grid, timer, condition_);
}

void Particles::correctVelocity(const Grid& grid, const Timer& timer, const Condition& condition) {
    correction_velocity.setZero();
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) != ParticleType::NORMAL) continue;
        if (boundary_types(i_particle) == BoundaryType::OTHERS) continue;
        Grid::Neighbors neighbors;
        grid.getNeighbors(i_particle, neighbors);
        double p_min = pressure(i_particle);
        for (int j_particle : neighbors) {
            if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
            p_min = std::min(pressure(j_particle), p_min);
        }
        Eigen::Vector3d tmp(0.0, 0.0, 0.0);
        for (int j_particle : neighbors) {
            if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
            Eigen::Vector3d r_ji = temporary_position.col(j_particle) - temporary_position.col(i_particle);
            tmp += r_ji * (pressure(j_particle) - p_min) * weightFunction(r_ji, grid.getGridWidth()) / r_ji.squaredNorm();
        }
        if (dimension == 2) tmp(2) = 0;
        correction_velocity.col(i_particle) -= tmp * dimension * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition.mass_density);
    }
    velocity = temporary_velocity + correction_velocity;
    position = temporary_position + correction_velocity * timer.getCurrentDeltaTime();
}

void Particles::checkSurfaceParticles() {
    checkSurfaceParticles(condition_.surface_parameter);
}

void Particles::checkSurfaceParticles(double surface_parameter) {
    for(int i = 0; i < getSize(); ++i) {
        if(particle_types(i) == ParticleType::NORMAL || particle_types(i) == ParticleType::WALL) {
            if (particle_number_density(i) < surface_parameter * initial_particle_number_density) boundary_types(i) = BoundaryType::SURFACE;
            else boundary_types(i) = BoundaryType::INNER;
        } else {
            boundary_types(i) = BoundaryType::OTHERS;
        }
    }
}

} // namespace tiny_mps