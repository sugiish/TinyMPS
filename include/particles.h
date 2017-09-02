#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <string>

#include <Eigen/Dense>

#include "condition.h"
#include "timer.h"
#include "grid.h"

using namespace Eigen;

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
	Matrix<double, 3, Dynamic> position;
	Matrix<double, 3, Dynamic> velocity;
	VectorXd pressure;
	VectorXd particle_number_density;

	Matrix<double, 3, Dynamic> temporary_position;
	Matrix<double, 3, Dynamic> temporary_velocity;

	VectorXi particle_types;

	Particles(const std::string& path, Condition& condition);
	virtual ~Particles();

	void updateParticleNumberDensity(Grid & grid, std::function<double(int, int)> weight);
	void setInitialParticleNumberDensity(int index);
	
	void moveParticlesExplicitly(const Vector3d& force, Timer timer);

	inline int getSize() const
	{
		return size;
	}

	inline void getMaxPosition(Vector3d& answer) const
	{
		answer = position.rowwise().maxCoeff();
	}

	inline void getMinPosition(Vector3d& answer) const
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

	void laplacianVelocity(int i_particle, int j_particle, Vector3d& output);

	inline double weightFunction(double r, double r_e)
	{
		return (0 <= r && r < r_e) ? 
			r / r_e - 1
			: 0;
	}

};

#endif // MPS_PARTICLES_H_INCLUDED