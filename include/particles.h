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
	void updateParticleNumberDensity(Grid & grid, std::function<double(int, int)> weight);
	void setInitialParticleNumberDensity(int index);
	void moveParticlesExplicitly(const Eigen::Vector3d& force, Timer timer);
	int writeVtkFile(const std::string& path, const std::string& title);	
	inline int getSize() const { return size; }

	Eigen::Matrix<double, 3, Eigen::Dynamic> position;
	Eigen::Matrix<double, 3, Eigen::Dynamic> velocity;
	Eigen::VectorXd pressure;
	Eigen::VectorXd particle_number_density;
	Eigen::Matrix<double, 3, Eigen::Dynamic> temporary_position;
	Eigen::Matrix<double, 3, Eigen::Dynamic> temporary_velocity;
	Eigen::VectorXi particle_types;

private:
	void initialize(int particles_number);
	int readGridFile(const std::string& path, int dimension);

	int size;
	double initial_particle_number_density;
};

} // namespace tiny_mps
#endif // MPS_PARTICLES_H_INCLUDED