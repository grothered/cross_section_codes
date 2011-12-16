compiler_options = -ffpe-trap=invalid,zero,overflow -fbounds-check 
#I'll try a makefile

clean: xsect
	rm *.o *.mod

xsect: global_defs.o util_various.o bed_xsect.o suspended_xsect.o hydro_xsect.o Pizzutotry_correction.o driver_xsect.f90
	gfortran  $(compiler_options) -o $@ $^ /usr/lib/liblapack.a libslatec.a

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


