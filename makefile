# Makefile for the xsect codes.

# Options for compiler optimization or testing
compiler_options =  -O3 
#compiler_options = -fbounds-check -g -pg

# Alternative to check for NaN's / overflow / underflow
# This one is a problem with the 1D suspended sediment routine, because the
# latter deliberately uses NaNs to detect overflow/underflow -- and I don't
# know a better method to use in the algorithm.
#compiler_options = -ffpe-trap=invalid,zero,overflow -fbounds-check -g -pg

# Clean up at the end of the make
clean: xsect q2d
	rm *.o *.mod

# This is to make the quasi2d model
q2d: global_defs.o util_various.o Pizzutotry_correction.o hydro_xsect.o bed_xsect.o susconc6.o st_venant_solver.o driver_q2d.f90  
	gfortran $(compiler_options) -o $@ $^ /usr/lib/liblapack.a libslatec.a

# This is to make the single cross-sectional code
xsect: global_defs.o util_various.o bed_xsect.o suspended_xsect.o hydro_xsect.o Pizzutotry_correction.o driver_xsect.f90
	gfortran  $(compiler_options) -o $@ $^ /usr/lib/liblapack.a libslatec.a

# Quasi-2D modules
st_venant_solver.o: st_venant_solver.f90 global_defs.o
	gfortran   -c $(compiler_options)   st_venant_solver.f90 

susconc6.o: susconc6.f90 global_defs.o
	gfortran   -c $(compiler_options)   susconc6.f90 

# Single cross-section modules
hydro_xsect.o: hydro_xsect.f90 Pizzutotry_correction.o global_defs.o
	gfortran -c  $(compiler_options)  hydro_xsect.f90

Pizzutotry_correction.o: Pizzutotry_correction.f90 global_defs.o
	gfortran -c  $(compiler_options)   Pizzutotry_correction.f90

suspended_xsect.o: suspended_xsect.f90 global_defs.o
	gfortran -c  $(compiler_options)   suspended_xsect.f90

bed_xsect.o: bed_xsect.f90 global_defs.o
	gfortran -c  $(compiler_options)   bed_xsect.f90

util_various.o: util_various.f90 global_defs.o
	gfortran -c  $(compiler_options)   util_various.f90

global_defs.o: global_defs.f90
	gfortran -c  $(compiler_options)  $^


