#ifndef MPS_PARTICLES_H_INCLUDED
#define MPS_PARTICLES_H_INCLUDED

#include<Eigen/Dense>

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

	VectorXi type;

	Particles(int particles_number, int dimension);
	virtual ~Particles();

private:
	int dimension;

};

#endif // MPS_PARTICLES_H_INCLUDED
