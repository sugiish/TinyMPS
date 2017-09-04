#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <string>

#include <Eigen/Dense>

#include "condition.h"
#include "timer.h"
#include "grid.h"

namespace tiny_mps {

enum ParticleType {
	NORMAL = 0,
	WALL = 2,
	DUMMY_WALL = 3,
	GHOST = -1
};

enum BoundaryType {
	INNER = 0,
	SURFACE = 1,
	OTHERS = -1
};

class Particles {
public:
	Particles(const std::string& path, const Condition& condition);
	virtual ~Particles();
	void updateParticleNumberDensity(Grid& grid);
	void updateParticleNumberDensity(Grid& grid, std::function<double(int, int)> weight);
	void calculateTemporaryVelocity(const Eigen::Vector3d& force, Grid& grid, const Timer& timer, const Condition& condition);
	void moveExplicitly();
	void solvePressurePoission();
	void checkSurfaceParticles(double surface_parameter);
	int writeVtkFile(const std::string& path, const std::string& title);
	inline double getMaxSpeed() {
		return velocity.colwise().norm().maxCoeff();
	}
	inline int getSize() const { return size; }

	Eigen::Matrix3Xd position;
	Eigen::Matrix3Xd velocity;
	Eigen::VectorXd pressure;
	Eigen::VectorXd particle_number_density;
	Eigen::Matrix3Xd temporary_position;
	Eigen::Matrix3Xd temporary_velocity;
	Eigen::VectorXi particle_types;
	Eigen::VectorXi boundary_types;

	double pnd_weight_radius;
	double laplacian_pressure_weight_radius;

private:
	using VectorXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;

	void initialize(int particles_number);
	int readGridFile(const std::string& path, int dimension);
	void calculateInitialParticleNumberDensity(int index);
	void calculateLaplacianLambda(int index, Grid& grid);

	double weightFunction(int i_particle, int j_particle, double influence_radius);
	double laplacianWeightWithNorm2(int i_particle, int j_particle);
	void laplacianViscosity(int i_particle, int j_particle, Eigen::Vector3d& output);

	std::function<double(int, int)> pnd_weight;
	int size;
	int dimension;
	double initial_particle_number_density;
	double laplacian_lambda;
};

} // namespace tiny_mps
#endif // MPS_PARTICLES_H_INCLUDED