#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <string>

#include <Eigen/Dense>

#include "condition.h"
#include "timer.h"
#include "grid.h"

namespace tiny_mps
{

enum ParticleType
{
	NORMAL = 0,
	WALL = 2,
	DUMMY_WALL = 3,
	GHOST = -1
};

class Particles
{
public:
	Eigen::Matrix<double, 3, Eigen::Dynamic> position;
	Eigen::Matrix<double, 3, Eigen::Dynamic> velocity;
	Eigen::VectorXd pressure;
	Eigen::VectorXd particle_number_density;

	Eigen::Matrix<double, 3, Eigen::Dynamic> temporary_position;
	Eigen::Matrix<double, 3, Eigen::Dynamic> temporary_velocity;

	Eigen::VectorXi particle_types;

	Particles(const std::string& path, Condition& condition);
	virtual ~Particles();

	void updateParticleNumberDensity(Grid & grid, std::function<double(int, int)> weight);
	void setInitialParticleNumberDensity(int index);
	
	void moveParticlesExplicitly(const Eigen::Vector3d& force, Timer timer);

	inline int getSize() const
	{
		return size;
	}

	inline void getMaxPosition(Eigen::Vector3d& answer) const
	{
		answer = position.rowwise().maxCoeff();
	}

	inline void getMinPosition(Eigen::Vector3d& answer) const
	{
		answer = position.rowwise().minCoeff();
	}

	int writeVtkFile(const std::string& path, const std::string& title);

private:
	Condition& condition;
	int size;
	double initial_particle_number_density;

	void initialize(int particles_number);
	int readGridFile(const std::string& path, int dimension);

	void laplacianVelocity(int i_particle, int j_particle, Eigen::Vector3d& output);

	inline double weightFunction(double r, double r_e)
	{
		return (0 <= r && r < r_e) ? 
			r / r_e - 1
			: 0;
	}

};

}
#endif // MPS_PARTICLES_H_INCLUDED