#I'll try a makefile

model: Pizzutotry_correction.o shionocros.o  shionoprog2.f90 
	gfortran -O3 -o $@ $^ lapack_gf_LINUX.a blas_gf_LINUX.a libslatec_gf.a
#	gfortran -O3 -o model Pizzutotry_correction.o shionocros.o shionoprog2.f90 lapack_gf_LINUX.a blas_gf_LINUX.a libslatec_gf.a
	
#Pizzutotry_correction.mod: Pizzutotry_correction.o Pizzutotry_correction.f90
#	gfortran    -c  -O3 Pizzutotry_correction.f90 

#shionomodmovs24.mod: Pizzutotry_correction.mod shionocros.o shionocros.f90
#	gfortran    -c  -O3 shionocros.f90


Pizzutotry_correction.o: Pizzutotry_correction.f90
	gfortran -c -O3 Pizzutotry_correction.f90
#	gfortran -c -O3 $^

shionocros.o: Pizzutotry_correction.o shionocros.f90
	gfortran -c -O3 shionocros.f90
#	gfortran -c -O3 $^

#driver2.o: 

#shionomodmovS2.mod:  matrix_solvers.o mrgrnk.o shionomodmovS2.o
#	gfortran    -c  


