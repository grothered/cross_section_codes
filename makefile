#I'll try a makefile

model: global_defs.o util_various.o bed_xsect.o suspended_xsect.o hydro_xsect.o Pizzutotry_correction.o driver_xsect.f90
	gfortran -O3  -o $@ $^ lapack_gf_LINUX.a blas_gf_LINUX.a libslatec_gf.a

hydro_xsect.o: hydro_xsect.f90 Pizzutotry_correction.o global_defs.o
	gfortran -c -O3  hydro_xsect.f90

Pizzutotry_correction.o: Pizzutotry_correction.f90 global_defs.o
	gfortran -c -O3  Pizzutotry_correction.f90

suspended_xsect.o: suspended_xsect.f90 global_defs.o
	gfortran -c -O3  suspended_xsect.f90

bed_xsect.o: bed_xsect.f90 global_defs.o
	gfortran -c -O3  bed_xsect.f90

util_various.o: util_various.f90 global_defs.o
	gfortran -c -O3  util_various.f90

global_defs.o: global_defs.f90
	gfortran -c -O3  $^


