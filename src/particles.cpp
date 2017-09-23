// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#include "particles.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>

namespace tiny_mps {

Particles::Particles(const std::string& path, const Condition& condition) : condition_(condition) {
    dimension = condition.dimension;
    readGridFile(path, condition);
    inflow_stride = 0.0;
    VectorXb valid = particle_types.array() != ParticleType::GHOST;
    updateParticleNumberDensity();
    setInitialParticleNumberDensity(condition.inner_particle_index);
    calculateLaplacianLambda(condition.inner_particle_index, condition);
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
    correction_velocity = Eigen::MatrixXd::Zero(3, size);
    neighbor_particles = Eigen::VectorXi::Zero(size);
}

int Particles::readGridFile(const std::string& path, const Condition& condition) {
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
    int ptcl_num = 0;
    {
        std::stringstream ss;
        ss.str(tmp_str);
        ss >> ptcl_num;
        initialize(ptcl_num + condition.extra_ghost_particles);
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
    for (int i_particle = ptcl_num; i_particle < size; ++i_particle) {
        particle_types(i_particle) = ParticleType::GHOST;
        ghost_stack.push(i_particle);
    }
    std::cout << "Succeed in reading grid file: " << path << std::endl;
    return 0;
}

int Particles::writeVtkFile(const std::string& path, const std::string& title) const {
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
    std::cout << "Succeed in writing vtk file: " << path << std::endl;
    return 0;
}

bool Particles::saveInterval(const std::string& path, const Timer& timer) const {
    if (!timer.isOutputTime()) return false;
    std::string output_index = (boost::format("%04d") % timer.getOutputCount()).str();
    writeVtkFile((boost::format(path) % output_index).str(), (boost::format("Time: %s") % timer.getCurrentTime()).str());
    return true;
}

void Particles::setInitialParticleNumberDensity(int index) {
    initial_particle_number_density = particle_number_density(index);
    initial_neighbor_particles = neighbor_particles(index);
    std::cout << "Initial particle number density: " << initial_particle_number_density << std::endl;
}

void Particles::calculateLaplacianLambda(int index, const Condition& condition) {
    VectorXb valid = particle_types.array() != ParticleType::GHOST;
    Grid lap_pressure(condition.laplacian_pressure_weight_radius, position, valid, condition.dimension);
    Grid::Neighbors neighbors;
    lap_pressure.getNeighbors(index, neighbors);
    double numerator = 0.0;
    double denominator = 0.0;
    for (int j_particle : neighbors) {
        Eigen::Vector3d r_ij = position.col(j_particle) - position.col(index);
        double w =  weightFunction(r_ij, lap_pressure.getGridWidth());
        numerator += r_ij.squaredNorm() * w;
        denominator += w;
    }
    laplacian_lambda_pressure = numerator / denominator;

    Grid lap_viscosity(condition.laplacian_viscosity_weight_radius, position, valid, condition.dimension);
    lap_pressure.getNeighbors(index, neighbors);
    for (int j_particle : neighbors) {
        Eigen::Vector3d r_ij = position.col(j_particle) - position.col(index);
        double w =  weightFunction(r_ij, lap_pressure.getGridWidth());
        numerator += r_ij.squaredNorm() * w;
        denominator += w;
    }
    laplacian_lambda_viscosity = numerator / denominator;
    std::cout << "Laplacian lambda for Pressure: " << laplacian_lambda_pressure << std::endl;
    std::cout << "Laplacian lambda for Viscosity: " << laplacian_lambda_viscosity << std::endl;
    std::cout << "Relaxation coefficient of lambda: " << condition.relaxation_coefficient_lambda << std::endl;
}

bool Particles::nextLoop(const std::string& path, Timer& timer) {
    std::cout << std::endl;
    timer.limitCurrentDeltaTime(getMaxSpeed(), condition_);
    timer.printCompuationTime();
    timer.printTimeInfo();
    std::cout << boost::format("Max velocity: %f") % getMaxSpeed() << std::endl;
    saveInterval(path, timer);
    if (checkNeedlessCalculation()) {
        std::cerr << "Error: All particles have become ghost." << std::endl;
        writeVtkFile((boost::format(path) % "err").str(), (boost::format("Time: %s") % timer.getCurrentTime()).str());
        exit(EXIT_FAILURE);
    }
    if (timer.isUnderMinDeltaTime()) {
        std::cerr << "Error: Delta time has become so small." << std::endl;
        writeVtkFile((boost::format(path) % "err").str(), (boost::format("Time: %s") % timer.getCurrentTime()).str());
        exit(EXIT_FAILURE);
    }
    if (!timer.hasNextLoop()) {
        std::cout << std::endl << "Total ";
        timer.printCompuationTime();
        std::cout << "Succeed in simulation." << std::endl;
        return false;
    }
    timer.update();
    temporary_velocity = velocity;
    temporary_position = position;
    return true;
}

bool Particles::checkNeedlessCalculation() const {
    if (pressure.hasNaN() || velocity.hasNaN() || position.hasNaN()) {
        std::cerr << "Error: Data contains NaN." << std::endl;
        return true;
    }
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::NORMAL) return false;
    }
    std::cerr << "Error: No normal particles." << std::endl;
    return true;
}

void Particles::setGhostParticle(int index) {
    if (index < 0 || index >= size) {
        std::cerr << "Error: Index is out of range." << std::endl;
        std::cerr << "Size: " << size << ", Index: " << std::endl;
    }
    particle_types(index) = ParticleType::GHOST;
    boundary_types(index) = BoundaryType::OTHERS;
    position.col(index).setZero();
    velocity.col(index).setZero();
    pressure(index) = 0.0;
    particle_number_density(index) = 0.0;
    neighbor_particles(index) = 0.0;
    temporary_position.col(index).setZero();
    temporary_velocity.col(index).setZero();
    correction_velocity.col(index).setZero();
    ghost_stack.push(index);
    std::cout << "Changed ghost particle: " << index << std::endl;
}

void Particles::removeOutsideParticles(const Eigen::Vector3d& minpos, const Eigen::Vector3d& maxpos) {
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::GHOST) continue;
        for (int i_dim = 0; i_dim < dimension; ++i_dim) {
            if (position.col(i_particle)(i_dim) < minpos(i_dim) || position.col(i_particle)(i_dim) > maxpos(i_dim)) {
                setGhostParticle(i_particle);
            }
        }
    }
}

void Particles::removeFastParticles(double max_speed) {
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::GHOST) continue;
        if (velocity.col(i_particle).norm() > max_speed) setGhostParticle(i_particle);
    }
}

void Particles::calculateTemporaryParticleNumberDensity() {
    Grid grid(condition_.pnd_weight_radius, temporary_position, particle_types.array() != ParticleType::GHOST, dimension);
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::GHOST) {
            particle_number_density(i_particle) = 0.0;
            neighbor_particles(i_particle) = 0;
            continue;
        }
        Grid::Neighbors neighbors;
        grid.getNeighbors(i_particle, neighbors);
        double pnd = 0.0;
        int count = 0;
        for (int j_particle : neighbors) {
            if (particle_types(i_particle) == ParticleType::GHOST) continue;
            Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
            pnd += weightFunction(r_ij, grid.getGridWidth());
            ++count;
        }
        particle_number_density(i_particle) = pnd;
        neighbor_particles(i_particle) = count;
    }
}

void Particles::updateParticleNumberDensity() {
    Grid grid(condition_.pnd_weight_radius, position, particle_types.array() != ParticleType::GHOST, dimension);
    updateParticleNumberDensity(grid);
}

void Particles::updateParticleNumberDensity(const Grid& grid) {
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::GHOST) {
            particle_number_density(i_particle) = 0.0;
            neighbor_particles(i_particle) = 0;
            continue;
        }
        Grid::Neighbors neighbors;
        grid.getNeighbors(i_particle, neighbors);
        double pnd = 0.0;
        int count = 0;
        for (int j_particle : neighbors) {
            if (particle_types(i_particle) == ParticleType::GHOST) continue;
            Eigen::Vector3d r_ij = position.col(j_particle) - position.col(i_particle);
            pnd += weightFunction(r_ij, grid.getGridWidth());
            ++count;
        }
        particle_number_density(i_particle) = pnd;
        neighbor_particles(i_particle) = count;
    }
}

void Particles::moveInflowParticles(const Timer& timer) {
    inflow_stride += condition_.inflow_velocity.norm() * timer.getCurrentDeltaTime();
    Eigen::Vector3d inflow_normalized = condition_.inflow_velocity.normalized();
    if (inflow_stride >= condition_.average_distance) {
        for (int i_particle = 0; i_particle < size; ++i_particle) {
            if (particle_types(i_particle) == ParticleType::INFLOW) {
                if (ghost_stack.empty()) {
                    std::cerr << "Error: Can't make new particles." << std::endl
                        << "Extra ghost particles has run out." << std::endl;
                    exit(EXIT_FAILURE);
                } else {
                    int new_index = ghost_stack.top();
                    ghost_stack.pop();
                    particle_types(new_index) = ParticleType::NORMAL;
                    boundary_types(new_index) = boundary_types(i_particle);
                    position.col(new_index) = position.col(i_particle);
                    velocity.col(new_index) = velocity.col(i_particle);
                    pressure(new_index) = pressure(i_particle);
                    particle_number_density(new_index) = particle_number_density(i_particle);
                    neighbor_particles(new_index) = neighbor_particles(i_particle);
                    temporary_position.col(new_index) = temporary_position.col(i_particle);
                    temporary_velocity.col(new_index) = temporary_velocity.col(i_particle);
                    correction_velocity.col(new_index) = correction_velocity.col(i_particle);
                    temporary_velocity.col(i_particle) = condition_.inflow_velocity;
                    velocity.col(i_particle) = condition_.inflow_velocity;
                    temporary_position.col(i_particle) -= inflow_normalized  * condition_.average_distance;
                    position.col(i_particle) -= inflow_normalized  * condition_.average_distance;
                }
            }
            if (particle_types(i_particle) == ParticleType::DUMMY_INFLOW) {
                temporary_velocity.col(i_particle) = condition_.inflow_velocity;
                velocity.col(i_particle) = condition_.inflow_velocity;
                temporary_position.col(i_particle) -= inflow_normalized * condition_.average_distance;
                position.col(i_particle) -= inflow_normalized * condition_.average_distance;
            }
        }
        inflow_stride -= condition_.average_distance;
    } else {
        for (int i_particle = 0; i_particle < size; ++i_particle) {
            if (particle_types(i_particle) == ParticleType::INFLOW || particle_types(i_particle) == ParticleType::DUMMY_INFLOW ) {
                temporary_velocity.col(i_particle) = condition_.inflow_velocity;
                velocity.col(i_particle) = condition_.inflow_velocity;
            }
        }
    }
}

void Particles::calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer) {
    Grid grid(condition_.laplacian_viscosity_weight_radius, position,
        particle_types.array() == ParticleType::NORMAL || particle_types.array() == ParticleType::INFLOW,
        condition_.dimension);
    calculateTemporaryVelocity(force, timer, grid);
}

void Particles::calculateTemporaryVelocity(const Eigen::Vector3d& force, const Timer& timer, Grid& grid) {
    double delta_time = timer.getCurrentDeltaTime();
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) == ParticleType::NORMAL) {
            temporary_velocity.col(i_particle) += delta_time * force;
            if (condition_.viscosity_calculation) {
                Grid::Neighbors neighbors;
                grid.getNeighbors(i_particle, neighbors);
                Eigen::Vector3d lap_vec(0.0, 0.0, 0.0);
                for (int j_particle : neighbors) {
                    Eigen::Vector3d u_ij = velocity.col(j_particle) - velocity.col(i_particle);
                    Eigen::Vector3d r_ij = position.col(j_particle) - position.col(i_particle);
                    lap_vec += u_ij * weightFunction(r_ij, grid.getGridWidth()) * 2 * dimension / (laplacian_lambda_viscosity * initial_particle_number_density);
                }
                temporary_velocity.col(i_particle) += lap_vec * condition_.kinematic_viscosity * delta_time;
            }
        }
    }
}

void Particles::updateTemporaryPosition(const Timer& timer) {
    temporary_position = position + timer.getCurrentDeltaTime() * temporary_velocity;
}

void Particles::solvePressurePoission(const Timer& timer) {
    Grid grid(condition_.laplacian_pressure_weight_radius, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    solvePressurePoission(timer, grid);
}

void Particles::solvePressurePoission(const Timer& timer, const Grid& grid) {
    using T = Eigen::Triplet<double>;
    double lap_r = grid.getGridWidth();
    int n_size = (int)(std::pow(lap_r * 2, dimension));
    double delta_time = timer.getCurrentDeltaTime();
    Eigen::SparseMatrix<double> p_mat(size, size);
    Eigen::VectorXd source(size);
    std::vector<T> coeffs(size * n_size);
    std::vector<int> neighbors(n_size * 2);
    source.setZero();
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
            Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
            double mat_ij = weightFunction(r_ij, lap_r) * 2 * dimension
                    / (laplacian_lambda_pressure * initial_particle_number_density);
            sum -= mat_ij;
            if (boundary_types(j_particle) == BoundaryType::INNER) {
                coeffs.push_back(T(i_particle, j_particle, mat_ij));
            }
        }
        sum -= condition_.weak_compressibility * condition_.mass_density / (delta_time * delta_time);
        coeffs.push_back(T(i_particle, i_particle, sum));
        source(i_particle) = - (particle_number_density(i_particle) - initial_particle_number_density) 
                    * condition_.relaxation_coefficient_lambda * condition_.mass_density
                    / (delta_time * delta_time * initial_particle_number_density);
    }
    p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
    
    // Solving a problem
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.compute(p_mat);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed decompostion." << std::endl;
    }
    pressure = cg.solve(source);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed solving." << std::endl;
    }
    std::cout << "Solver - iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    for (int i = 0; i < size; ++i) {
        if (pressure(i) < 0) pressure(i) = 0;
    }
}

void Particles::solveTanakaMasunagaPressurePoission(const Timer& timer) {
    Grid grid(condition_.laplacian_pressure_weight_radius, position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    using T = Eigen::Triplet<double>;
    double lap_r = grid.getGridWidth();
    int n_size = (int)(std::pow(lap_r * 2, dimension));
    double delta_time = timer.getCurrentDeltaTime();
    Eigen::SparseMatrix<double> p_mat(size, size);
    Eigen::VectorXd source(size);
    std::vector<T> coeffs(size * n_size);
    std::vector<int> neighbors(n_size * 2);
    source.setZero();
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
            Eigen::Vector3d r_ij = position.col(j_particle) - position.col(i_particle);
            double mat_ij = weightFunction(r_ij, lap_r) * 2 * dimension
                    / (laplacian_lambda_pressure * initial_particle_number_density);
            sum -= mat_ij * condition_.tanaka_masunaga_c;
            div_vel += (temporary_velocity.col(j_particle) - temporary_velocity.col(i_particle)).dot(r_ij) 
                    * weightFunction(r_ij, lap_r) * condition_.dimension / (r_ij.squaredNorm() * initial_particle_number_density);
            if (boundary_types(j_particle) == BoundaryType::INNER) {
                coeffs.push_back(T(i_particle, j_particle, mat_ij));
            }
        }
        coeffs.push_back(T(i_particle, i_particle, sum));
        source(i_particle) = div_vel * condition_.mass_density / delta_time 
                    - (particle_number_density(i_particle) - initial_particle_number_density) 
                    * condition_.tanaka_masunaga_gamma * condition_.mass_density / (delta_time * delta_time * initial_particle_number_density);
    }
    p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
    
    // Solving a problem
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.compute(p_mat);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed decompostion." << std::endl;
    }
    pressure = cg.solve(source);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed solving." << std::endl;
    }
    std::cout << "Solver - iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    for (int i = 0; i < size; ++i) {
        if (pressure(i) < 0) pressure(i) = 0;
    }
}

void Particles::solveConjugateGradient(Eigen::SparseMatrix<double> p_mat, Eigen::VectorXd source) {
    Eigen::VectorXd solutions;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.compute(p_mat);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed decompostion." << std::endl;
    }
    solutions = cg.solve(source);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed solving." << std::endl;
    }
    std::cout << "Solver - iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    for (int i = 0; i < size; ++i) {
        if (pressure(i) < 0) pressure(i) = 0;
    }
}

void Particles::solvePressurePoissionWithTanakaMasunaga(const Timer& timer) {
    Grid grid(condition_.laplacian_pressure_weight_radius, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    using T = Eigen::Triplet<double>;
    double lap_r = grid.getGridWidth();
    int n_size = (int)(std::pow(lap_r * 2, dimension));
    double delta_time = timer.getCurrentDeltaTime();
    Eigen::SparseMatrix<double> p_mat(size, size);
    Eigen::VectorXd source(size);
    std::vector<T> coeffs(size * n_size);
    std::vector<int> neighbors(n_size * 2);
    source.setZero();
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
            double mat_ij = weightFunction(r_ij, lap_r) * 2 * dimension
                    / (laplacian_lambda_pressure * initial_particle_number_density);
            sum -= mat_ij * condition_.tanaka_masunaga_c;
            div_vel += (temporary_velocity.col(j_particle) - temporary_velocity.col(i_particle)).dot(r_ij) 
                    * weightFunction(r_ij, lap_r) * condition_.dimension / (r_ij.squaredNorm() * initial_particle_number_density);
            if (boundary_types(j_particle) == BoundaryType::INNER) {
                coeffs.push_back(T(i_particle, j_particle, mat_ij));
            }
        }
        coeffs.push_back(T(i_particle, i_particle, sum));
        source(i_particle) = div_vel * condition_.mass_density / delta_time 
                    - (particle_number_density(i_particle) - initial_particle_number_density) 
                    * condition_.tanaka_masunaga_gamma * condition_.mass_density / (delta_time * delta_time * initial_particle_number_density);
    }
    p_mat.setFromTriplets(coeffs.begin(), coeffs.end()); // Finished setup matrix
    
    // Solving a problem
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.compute(p_mat);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed decompostion." << std::endl;
    }
    pressure = cg.solve(source);
    if (cg.info() != Eigen::ComputationInfo::Success) {
        std::cerr << "Error: Failed solving." << std::endl;
    }
    std::cout << "Solver - iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    for (int i = 0; i < size; ++i) {
        if (pressure(i) < 0) pressure(i) = 0;
    }
}

void Particles::correctVelocity(const Timer& timer) {
    Grid grid(condition_.gradient_radius, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    correctVelocity(timer, grid);
}

void Particles::correctVelocity(const Timer& timer, const Grid& grid) {
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
            Eigen::Vector3d r_ij = temporary_position.col(j_particle) - temporary_position.col(i_particle);
            tmp += r_ij * (pressure(j_particle) - p_min) * weightFunction(r_ij, grid.getGridWidth()) / r_ij.squaredNorm();
        }
        if (dimension == 2) tmp(2) = 0;
        correction_velocity.col(i_particle) -= tmp * dimension * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition_.mass_density);
    }
    temporary_velocity += correction_velocity;
}

void Particles::correctVelocityExplicitly(const Timer& timer) {
    Grid grid(condition_.gradient_radius, position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
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
            Eigen::Vector3d r_ij = position.col(j_particle) - position.col(i_particle);
            tmp += r_ij * (pressure(j_particle) - p_min) * weightFunction(r_ij, grid.getGridWidth()) / r_ij.squaredNorm();
        }
        if (dimension == 2) tmp(2) = 0;
        correction_velocity.col(i_particle) -= tmp * dimension * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition_.mass_density);
    }
    temporary_velocity += correction_velocity;
    temporary_position += correction_velocity * timer.getCurrentDeltaTime();
}

void Particles::correctTanakaMasunagaVelocity(const Timer& timer) {
    Grid grid(condition_.gradient_radius, position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    correction_velocity.setZero();
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) != ParticleType::NORMAL) continue;
        if (boundary_types(i_particle) == BoundaryType::OTHERS) continue;
        Grid::Neighbors neighbors;
        grid.getNeighbors(i_particle, neighbors);
        Eigen::Vector3d tmp(0.0, 0.0, 0.0);
        for (int j_particle : neighbors) {
            if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
            Eigen::Vector3d r_ij = position.col(j_particle) - position.col(i_particle);
            tmp += r_ij * (pressure(j_particle) + pressure(i_particle)) * weightFunction(r_ij, grid.getGridWidth()) / r_ij.squaredNorm();
        }
        if (dimension == 2) tmp(2) = 0;
        correction_velocity.col(i_particle) -= tmp * dimension * timer.getCurrentDeltaTime() / (initial_particle_number_density * condition_.mass_density);
    }
    temporary_velocity += correction_velocity;
}

void Particles::updateVelocityAndPosition() {
    velocity = temporary_velocity;
    position = temporary_position;
}

void Particles::checkSurfaceParticles() {
    checkSurfaceParticles(condition_.surface_parameter);
}

void Particles::checkSurfaceParticles(double surface_parameter) {
    for(int i = 0; i < getSize(); ++i) {
        if (particle_types(i) == ParticleType::NORMAL || particle_types(i) == ParticleType::WALL || particle_types(i) == ParticleType::INFLOW) {
            if (particle_number_density(i) < surface_parameter * initial_particle_number_density) boundary_types(i) = BoundaryType::SURFACE;
            else boundary_types(i) = BoundaryType::INNER;
        } else {
            boundary_types(i) = BoundaryType::OTHERS;
        }
    }
}

void Particles::checkSurfaceParticlesWithTanakaMasunaga() {
    for(int i = 0; i < getSize(); ++i) {
        if (particle_types(i) == ParticleType::NORMAL || particle_types(i) == ParticleType::WALL
            || particle_types(i) == ParticleType::INFLOW) {
            if (particle_number_density(i) < condition_.surface_parameter * initial_particle_number_density 
                    && neighbor_particles(i) < condition_.tanaka_masunaga_beta * initial_neighbor_particles) {
                boundary_types(i) = BoundaryType::SURFACE;
            } else {
                boundary_types(i) = BoundaryType::INNER;
            }
        } else {
            boundary_types(i) = BoundaryType::OTHERS;
        }
    }
}

void Particles::checkTanakaMasunagaSurfaceParticles(double surface_parameter) {
    for(int i = 0; i < getSize(); ++i) {
        if (particle_types(i) == ParticleType::NORMAL || particle_types(i) == ParticleType::WALL || particle_types(i) == ParticleType::INFLOW) {
            if (neighbor_particles(i) < surface_parameter * initial_neighbor_particles) boundary_types(i) = BoundaryType::SURFACE;
            else boundary_types(i) = BoundaryType::INNER;
        } else {
            boundary_types(i) = BoundaryType::OTHERS;
        }
    }
}

void Particles::giveCollisionRepulsionForce() {
    giveCollisionRepulsionForce(condition_.collision_influence, condition_.restitution_coefficent);
}

void Particles::giveCollisionRepulsionForce(double influence_ratio, double restitution_coefficient) {
    Grid grid(influence_ratio * condition_.average_distance, temporary_position, boundary_types.array() != BoundaryType::OTHERS, condition_.dimension);
    Eigen::Matrix3Xd impulse_vel = Eigen::MatrixXd::Zero(3, size);
    for (int i_particle = 0; i_particle < size; ++i_particle) {
        if (particle_types(i_particle) != ParticleType::NORMAL) continue;
        if (boundary_types(i_particle) == BoundaryType::OTHERS) continue;
        Grid::Neighbors neighbors;
        grid.getNeighbors(i_particle, neighbors);
        for (int j_particle : neighbors) {
            if (boundary_types(j_particle) == BoundaryType::OTHERS) continue;
            Eigen::Vector3d n_ij = (temporary_position.col(j_particle) - temporary_position.col(i_particle)).normalized();
            Eigen::Vector3d u_ij = temporary_velocity.col(j_particle) - temporary_velocity.col(i_particle);
            impulse_vel.col(i_particle) += n_ij * u_ij.dot(n_ij) * (restitution_coefficient + 1) / 2;
        }
    }
    temporary_velocity += impulse_vel;
}

void Particles::showParticlesInfo() {
    
}

} // namespace tiny_mps