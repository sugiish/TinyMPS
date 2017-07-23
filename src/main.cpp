#include "particles.h"
#include "grid.h"

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
    pt = new Particles(10, 2);

    Particles* pt2;
    pt2 = new Particles("./input/input.grid", 3);
    pt2->timer.setInitialDeltaTime(0.0002);
    //Grid grid(*pt2, 0.0002);

    Particles pt3("./input/input.grid", 2);
    pt3.timer.setInitialDeltaTime(0.00001);

    Grid g(pt3, 0.0002);
}