#!/bin/bash

# Yalimbah Test
cd yalimbah_test
cp ../../q2d .
./q2d < inputdata2.modin
R CMD BATCH r_compare.R /dev/tty
cd ..

#
cd tambrioni_test
cp ../../q2d .
./q2d < inputdata2.modin
R CMD BATCH r_compare.R /dev/tty
cd ..
