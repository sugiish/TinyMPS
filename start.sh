#!/bin/sh

seq 1 8 | parallel mkdir -p ./output{}
cat ./input/inflow | parallel --colsep ' ' --results ./parallel_out --joblog parallel.log "./bin/cavitation_analysis ./output{1}/ -{2}"
