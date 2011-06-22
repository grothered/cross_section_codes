#I'll try a makefile

model: global_defs.o util_various.o bed_xsect.o suspended_xsect.o hydro_xsect.o Pizzutotry_correction.o driver_xsect.f90 
	ifort -check bounds  -o $@ $^ lapack_if_LINUX.a blas_if_LINUX.a libslatec_if.a

hydro_xsect.o: hydro_xsect.f90 Pizzutotry_correction.o global_defs.o
	ifort -c -check bounds  hydro_xsect.f90

Pizzutotry_correction.o: Pizzutotry_correction.f90 global_defs.o
	ifort -c -check bounds  Pizzutotry_correction.f90

suspended_xsect.o: suspended_xsect.f90 global_defs.o
	ifort -c -check bounds  suspended_xsect.f90

bed_xsect.o: bed_xsect.f90 global_defs.o
	ifort -c -check bounds  bed_xsect.f90

util_various.o: util_various.f90 global_defs.o
	ifort -c -check bounds  util_various.f90

global_defs.o: global_defs.f90
	ifort -c -check bounds  $^



