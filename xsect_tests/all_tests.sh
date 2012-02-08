## Automatically run tests for cross-sectional code

# Analytical flow solution
cd skm_analytical
cp ../../xsect .
./xsect < skm_analytical.f90
R CMD BATCH R_analytical.R
eog Analytical_compare.png
cd ..


