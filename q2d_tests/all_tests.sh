#!/bin/bash

# Yalimbah Test
cd yalimbah_test
cp ../../q2d .
echo '#########################'
echo 'Running Yalimbah test ...'
echo 'This just compares the output with a previous model run, 
      to check for significant changes'
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave r_compare.R /dev/tty
cd ..

# Tambrioni Test
cd tambrioni_test
cp ../../q2d .
echo '#########################'
echo 'Running Tambrioni test ...'
echo 'This just compares the output with a previous model run, 
      to check for significant changes'
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave r_compare.R /dev/tty
cd ..

# Simple channel test
cd simple_channel
cp ../../q2d .
echo '#########################'
echo 'Running Simple Channel test ...'
echo 'This just compares the output with a previous model run, 
      to check for significant changes'
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave r_compare.R /dev/tty
cd ..

# Steady uniform test
cd steady_uniform 
cp ../../q2d .
echo '#########################'
echo 'Running Steady uniform test ...'
echo ' This compares with an analytical solution'
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave analytical.R /dev/tty
cd ..

# Dam break (wet downstream region)
cd dam_break/wet
cp ../../../q2d .
echo '#########################'
echo 'Running Wet dam break ...'
echo ' This compares with an analytical solution'
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave analytical_compare2.R 
evince Dam_break_wet.eps &
cd ../../

# FIXME: Dam break (dry downstream region)

# FIXME: Non-uniform hydraulic jump

# 1D suspended sediment tests
cd sus1D/river_inflow
cp ../../../q2d .
echo '#########################'
echo 'Running river_input suspended sediment  ...'
echo ' This compares the with an analytical solution
        of the advection diffusion equation'
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave Rcompare.R 
evince compare.pdf &
cd ../../../

