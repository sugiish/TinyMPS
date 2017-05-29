#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include<Eigen/Dense>

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

	Particles(int particles_number);
	~Partricles();

};

#endif // MPS_PARTICLES_H_INCLUDED
