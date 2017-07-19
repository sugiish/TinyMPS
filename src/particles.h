#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include <string>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;
using namespace std;
static const int DimX = 0;
static const int DimY = 1;
static const int DimZ = 2;

class Particles
{
public:
	int particles_number;
	MatrixXd position;
	MatrixXd velocity;
	VectorXd pressure;

	MatrixXd temporal_position;
	MatrixXd temporal_velocity;

	VectorXi particles_type;
	VectorXi ghost_particles;

	Timer timer;

	Particles(int particles_number, int dimension);
	Particles(string path, int dimension);
	virtual ~Particles();

private:
	int dimension;

	void initialize(int particles_number, int dimension);
	int readGridFile(string path, int dimension);

};

#endif // MPS_PARTICLES_H_INCLUDED
