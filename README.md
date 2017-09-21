# TinyMPS
This is a simple library for Moving Particle Semi-implicit (MPS) method.

## Requirements
- Eigen (http://eigen.tuxfamily.org/)
- Boost C++ Libraries (http://www.boost.org)

## Usage
Build the source code.
'''bash
make
'''

Run example codes for dam break problem.
'''bash
mkdir output
./bin/standard_mps
'''
This example code writes vtk files as output.
You can visualize them with Paraview (https://www.paraview.org).

## License
Copyright (c) 2017 Shota SUGIHARA

Distributed under the [MIT License](LICENSE).