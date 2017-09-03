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

class Particles {
public:
	Particles(const std::string& path, const Condition& condition);
	virtual ~Particles();
	void updateParticleNumberDensity(Grid& grid);
	void updateParticleNumberDensity(Grid& grid, std::function<double(int, int)> weight);
	void calculateInitialParticleNumberDensity(int index);
	void calculateLaplacianLambda(int index, Grid& grid);
	void moveParticlesExplicitly(const Eigen::Vector3d& force, Timer timer);
	int writeVtkFile(const std::string& path, const std::string& title);	
	inline int getSize() const { return size; }

	Eigen::Matrix3Xd position;
	Eigen::Matrix3Xd velocity;
	Eigen::VectorXd pressure;
	Eigen::VectorXd particle_number_density;
	Eigen::Matrix3Xd temporary_position;
	Eigen::Matrix3Xd temporary_velocity;
	Eigen::VectorXi particle_types;

private:
	using VectorXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;

	void initialize(int particles_number);
	int readGridFile(const std::string& path, int dimension);
	double weightFunction(int i_particle, int j_particle);

	std::function<double(int, int)> weight_function;
	int size;
	double weight_function_radius;
	double initial_particle_number_density;
	double laplacian_lambda;
};

} // namespace tiny_mps
#endif // MPS_PARTICLES_H_INCLUDED