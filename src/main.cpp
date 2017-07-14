#include "particles.h"

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

/**
 * This is a sample code using mps library.
 */
int main()
{
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;

    // temp

    Particles* pt;
    pt = new Particles(10, 3);

}