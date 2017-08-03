#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <string>

#include <Eigen/Dense>

#include "condition.h"
#include "timer.h"

using namespace Eigen;
using namespace std;

enum ParticleType
{
	
};

class Particles
{
public:
	Matrix<double, 3, Dynamic> position;
	Matrix<double, 3, Dynamic> velocity;
	VectorXd pressure;

	Matrix<double, 3, Dynamic> temporary_position;
	Matrix<double, 3, Dynamic> temporary_velocity;

	VectorXi particles_type;
	VectorXi particles_valid;

	Timer timer;

	Particles(const string& path, Condition& condition);
	virtual ~Particles();

	void moveParticlesExplicitly(double delta_time, const VectorXd& force);

	inline int getParticlesNumber() const
	{
		return particles_number;
	}

	inline void getMaxPosition(Vector3d& answer) const
	{
		answer = position.rowwise().maxCoeff();
	}

	inline void getMinPosition(Vector3d& answer) const
	{
		answer = position.rowwise().minCoeff();
	}

private:
	Condition& condition;
	int particles_number;

	void initialize(int particles_number);
	int readGridFile(const string& path, int dimension);

};

#endif // MPS_PARTICLES_H_INCLUDED
