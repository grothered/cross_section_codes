#I'll try a makefile

model: global_defs.o util_various.o bed_xsect.o suspended_xsect.o hydro_xsect.o Pizzutotry_correction.o driver_xsect.f90
	gfortran  -ffpe-trap=invalid,zero,overflow -o $@ $^ /usr/lib/liblapack.a libslatec_gf.a

hydro_xsect.o: hydro_xsect.f90 Pizzutotry_correction.o global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow    hydro_xsect.f90

Pizzutotry_correction.o: Pizzutotry_correction.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow    Pizzutotry_correction.f90

suspended_xsect.o: suspended_xsect.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow    suspended_xsect.f90

bed_xsect.o: bed_xsect.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow    bed_xsect.f90

util_various.o: util_various.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow    util_various.f90

global_defs.o: global_defs.f90
	gfortran -c  -ffpe-trap=invalid,zero,overflow    $^


