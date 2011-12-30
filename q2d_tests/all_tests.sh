#!/bin/bash

# Yalimbah Test
cd yalimbah_test
cp ../../q2d .
echo '#########################'
echo 'Running Yalimbah test ...'
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
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave analytical_compare2.R 
evince Dam_break_wet.eps
cd ../../

# FIXME: Dam break (dry downstream region)

# FIXME: 1D suspended sediment tests
cd sus1D
cp ../../q2d .
echo '#########################'
echo 'Running 1D sus  ...'
echo '#########################'
./q2d < inputdata2.modin > outfile.log
tail -n1 outfile.log
R CMD BATCH --slave Rcompare.R 
evince compare.pdf
cd ../../

