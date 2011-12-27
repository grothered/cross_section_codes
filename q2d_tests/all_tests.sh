#!/bin/bash

# Yalimbah Test
cd yalimbah_test
cp ../../q2d .
echo 'Running Yalimbah test ...'
./q2d < inputdata2.modin > outfile.log
R CMD BATCH --slave r_compare.R /dev/tty
cd ..

# Tambrioni Test
cd tambrioni_test
cp ../../q2d .
echo 'Running Tambrioni test ...'
./q2d < inputdata2.modin > outfile.log
R CMD BATCH --slave r_compare.R /dev/tty
cd ..

# Simple channel test
cd simple_channel
cp ../../q2d .
echo 'Running Simple Channel test ...'
./q2d < inputdata2.modin > outfile.log
R CMD BATCH --slave r_compare.R /dev/tty
cd ..

