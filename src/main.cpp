#include "particles.h"
#include "grid.h"
#include "reader.h"

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

    Particles pt3("./input/dambreak.grid", 2);
    pt3.timer.setInitialDeltaTime(0.00001);

    Grid g(pt3, 0.0002);
    int a,b;
    g.getNeighbor(3582, a, b);
    std::cout << a << "," << b << std::endl;

    /*Reader reader("./input/input.data");
    reader.printValues();
    bool abst = false;
    reader.getValue("outputPressureFile", abst);
    std::cout << abst << std::endl;
    reader.getValue("autoSettingOfCapacityOfNeighborTable", abst);
    std::cout << abst << std::endl;
    std::cout << true << std::endl;

    double tmp_d;
    reader.getValue("gassConstant", tmp_d);
    std::cout << tmp_d << std::endl;
    reader.getValue("compressibilityOfType0", tmp_d);
    std::cout << tmp_d << std::endl;
*/
    Reader reader2("./input/new_input.data");
    reader2.printValues();

    double average, compressibility;
    reader2.getValue("finish_time", average);
    reader2.getValue("compressibility", compressibility);
    std::cout << average << std::endl;
    std::cout << compressibility << std::endl;

}