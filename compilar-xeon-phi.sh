#!/bin/bash

CXXFLAGS="-mmic "
CFLAGS="-mmic "
make -j MKLROOT=/opt/intel/mkl/ INTEL=/opt/intel/lib/mic cpu=1 time=1 openmp=1 xeon_phi=1 cpu_recompute=0
