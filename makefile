#I'll try a makefile

xsect: global_defs.o util_various.o bed_xsect.o suspended_xsect.o hydro_xsect.o Pizzutotry_correction.o driver_xsect.f90
	gfortran  -ffpe-trap=invalid,zero,overflow -fbounds-check -o $@ $^ /usr/lib/liblapack.a libslatec.a

hydro_xsect.o: hydro_xsect.f90 Pizzutotry_correction.o global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow -fbounds-check    hydro_xsect.f90

Pizzutotry_correction.o: Pizzutotry_correction.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow -fbounds-check    Pizzutotry_correction.f90

suspended_xsect.o: suspended_xsect.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow -fbounds-check    suspended_xsect.f90

bed_xsect.o: bed_xsect.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow -fbounds-check    bed_xsect.f90

util_various.o: util_various.f90 global_defs.o
	gfortran -c  -ffpe-trap=invalid,zero,overflow -fbounds-check    util_various.f90

global_defs.o: global_defs.f90
	gfortran -c  -ffpe-trap=invalid,zero,overflow -fbounds-check    $^


clean: *.o *.mod
	rm *.o *.mod
