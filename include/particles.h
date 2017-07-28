#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <string>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;
using namespace std;

class Particles
{
public:
	MatrixXd position;
	MatrixXd velocity;
	VectorXd pressure;

	MatrixXd temporal_position;
	MatrixXd temporal_velocity;

	VectorXi particles_type;
	VectorXi particles_valid;

	Timer timer;

	Particles(int particles_number, int dimension);
	Particles(const string& path, int dimension);
	virtual ~Particles();

	void moveParticlesExplicitly(double delta_time, const VectorXd& force);

	inline int getParticlesNumber() const
	{
		return particles_number;
	}

	inline int getDimension() const
	{
		return dimension;
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
	int particles_number;
	int dimension;

	void initialize(int particles_number, int dimension);
	int readGridFile(const string& path, int dimension);

};

#endif // MPS_PARTICLES_H_INCLUDED
