In this source tree, the cross-sectional morphodynamics program consists of the
following files:

driver_xsect.f90 -- contains the main program to run everything

global_defs.f90 -- contains dp, pi (global parameters)

hydro_xsect.f90 -- contains routines for single cross-section hydrodynamics,
    e.g. shear, calc_friction, rough_mult, 

Pizzutotry_correction.f90 -- countains an alternative hydrodynamic routine.
    Presently (3/7/2011) this code is outdated and not used much or tested.

suspended_xsect.f90 -- contains routines for single cross-section suspended
    sediment distribution

bed_xsect.f90 -- contains routines to compute rates of sediment transport and
    morphological evolution

util_various.f90 -- contains various utility routines that did not sit neatly
    elsewhere

xsect_data.modin -- a fortran namelist which contains input constants. 

makefile -- to stick it all together

We also rely on a lapack library to support various matrix
computations. 

A slatec library is also used, presently (3/7/2011) only for its
cubic spline interpolation, which is used in remeshing (often I run the code
with a static mesh, and so slatec is not used at all). For gfortran, the source
for this can be downloaded from (11/2011):
http://joachimwuttke.de/slatec4gfortran/


There are also utility codes for plotting / analysis / various checks.

cros.R -- contains routines to read and plot outputs in R
wset.R -- contains various useful calculations. 


And there are sub-directories which contain alternative 'inputdata.f90'
namelists, corresponding to different cases of interest. 

There are also various tests sitting around. 
