#include <iostream>

#include <Eigen/Dense>

#include "condition.h"
#include "grid.h"
#include "particles.h"
#include "reader.h"

using namespace Eigen;
using namespace std;

/**
 * This is a sample code using mps library.
 */
int main()
{
    Reader reader("./input/input.data");
    Condition condition(reader);
    Particles particles("./input/dambreak.grid", condition);
    
    while(particles.timer.hasNextLoop())
    {
        break;

    }
}