#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
export OMP_NUM_THREADS=11
export MKL_NUM_THREADS=1
export OMP_NESTED=true
export MKL_DYNAMIC=true
export KMP_AFFINITY=granularity=core,scatter
export KMP_BLOCKTIME=0
$LIOBIN -i hemo.in -b DZVP  -c hem.xyz -v | tee salida
