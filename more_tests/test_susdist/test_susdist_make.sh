#!/bin/bash
gfortran -o test_susdist -pg -fbounds-check test_susdist.f90 lapack_gf_LINUX.a blas_gf_LINUX.a libslatec_gf.a
