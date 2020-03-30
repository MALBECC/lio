#!/bin/bash


file=pyra
touch PYR.*; rm PYR.*; touch INPUT*; rm INPUT*
export LIOHOME=/home/gonzalo/progs/LIOs/PR_devTSH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOHOME/g2g:$LIOHOME/lioamber
export GFORTRAN_UNBUFFERED_ALL=1
export OMP_NUM_THREADS=3

HYB=/home/gonzalo/progs/hybrid_TSH_to_uploaded/bin
$HYB/hybrid < pyra.fdf > pyra.out

