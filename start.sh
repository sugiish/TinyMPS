#!/bin/sh

seq 1 8 | parallel mkdir -p ./output/nozzle_{}
cat ./input/inflow | parallel --colsep ' ' --results ./nozzle_out --joblog parallel_nozzle.log "./bin/cavitation_analysis ./output/nozzle_{1}/ -{2}"
#seq 1 8 | parallel mkdir -p ./output/spray_{}
#cat ./input/inflow | parallel --colsep ' ' --results ./spray_out --joblog parallel_spray.log "./bin/cavitation_spray ./output/spray_{1}/ -{2}"
