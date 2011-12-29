This source tree contains routines for 2 programs:

1) A quasi-2D hydrodynamic/morphodynamic solver for a single thread channel.

   The main program is driver_q2d.f90
   It depends on bed_xsect.f90, global_defs.f90, hydro_xsect.f90,
   Pizzutotry_correction.f90, st_venant_solver.f90, susconc6.f90,
   util_various.f90

   It produces output files with extensions of the form .0DO, .1DO, .2DO.

2) A single cross-section morphodynamic/hydrodynamic model.

   The main program is driver_xsect.f90
   It depends on bed_xsect.f90, global_defs.f90, hydro_xsect.f90,
   Pizzutotry_correction.f90, susconc_xsect.f90,
   util_various.f90

   It produces output files with extensions of the form .cst


Summaries of these routines, and other important files, are as follows:

driver_q2d.f90 -- contains the main program to run everything for the quasi-2D
    case.

st_venant_solver.f90 -- contains a solver for the St Venant Equations for the
    quasi-2D case.

susconc6.f90 -- contains solvers for the concentration of suspended sediment in
    the quasi-2D case.

driver_xsect.f90 -- contains the main program to run everything for the single
    cross-section case

global_defs.f90 -- contains dp, pi etc (global parameters)

hydro_xsect.f90 -- contains routines for single cross-section hydrodynamics,
    e.g. shear, calc_friction, rough_mult, 

Pizzutotry_correction.f90 -- countains an alternative hydrodynamic routine to
    the one used in hydro_xsect.f90. Presently (3/7/2011) this code is outdated and
    not used much or tested.

suspended_xsect.f90 -- contains routines for single cross-section suspended
    sediment distribution

bed_xsect.f90 -- contains routines to compute rates of sediment transport and
    morphological evolution

util_various.f90 -- contains various utility routines that did not sit neatly
    elsewhere

xsect_data.modin -- a fortran namelist which contains input constants for the
    single cross-section case.

q2d_data.modin -- a fortran namelist which contains input constants for the
    quasi 2D case.

makefile -- file to compile both programs

We also rely on a lapack library to support various matrix computations. 

A slatec library is also used, presently (3/7/2011) only for its cubic spline
interpolation, which is used in remeshing (often I run the code with a static
mesh, and so slatec is not used at all). The source for this can be downloaded
from (11/2011):
http://joachimwuttke.de/slatec4gfortran/
, and this compiles straightforwardly with gfortran.


There are also utility codes for plotting / analysis / various checks.

cros.R -- contains routines to read and plot single cross-section outputs in R
wset.R -- contains various useful calculations. 
gettide1.R -- contains various utility routines for plotting the quasi2d outputs in R

And there are sub-directories which contain alternative cases of interest.

There are also various tests sitting around. 


HOW TO COMPILE AND RUN THE CODE.

You need a fortran compiler, and a 'make' program . You need to have the lapack
and slatec libraries compiled into shared libraries (with a '.a' extension).
These libraries and the compiler name are referred to in the 'makefile', so you
may need to edit the latter so that the filenames/paths/compiler name are
correct. If all this is set up correctly, then you should be able to type:
 make
to compile the files, and then:
 ./q2d < q2d_data.modin
to run a quasi-2d case, and:
 ./xsect < xsect_data.modin
to run a single cross-section case.

Beyond this, you'll need to read the source to understand what's being done.
Start with looking at the driver routines and the .modin namelists, and move on
from there. 
I typically use R to plot/investigate the outputs, using the .R routines in the
source tree, but you could use any software that can read the ascii output
files.



OUTSTANDING ISSUES
1) Set up boundary conditions to accommodate sub and super critical flow
    - Then add in the Zoppou Roberts hydraulic jump test

