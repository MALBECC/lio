#!/bin/bash

# To install the gsl execute gsl-icc-mic.sh in the xeon-phi folder
CXXFLAGS="-mmic "
CFLAGS="-mmic "
make -j MKLROOT=/opt/intel/mkl/ INTEL=/opt/intel/lib/mic cpu=1 time=1 openmp=1 xeon_phi=1
