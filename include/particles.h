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
	MatrixXd position;
	MatrixXd velocity;
	VectorXd pressure;

	MatrixXd temporary_position;
	MatrixXd temporary_velocity;

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

	inline void getMaxPosition(VectorXd& answer) const
	{
		answer = position.rowwise().maxCoeff();
	}

	inline void getMinPosition(VectorXd& answer) const
	{
		answer = position.rowwise().minCoeff();
	}

private:
	Condition& condition;
	int particles_number;

	void initialize(int particles_number, int dimension);
	int readGridFile(const string& path, int dimension);

};

#endif // MPS_PARTICLES_H_INCLUDED
