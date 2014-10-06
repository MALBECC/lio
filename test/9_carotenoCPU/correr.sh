#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=1
export KMP_AFFINITY=granularity=core,scatter
export LIO_OUTER_THREADS=12
export LIO_INNER_THREADS=1
$LIOBIN -i carotenox.in -b DZVP  -c caroteno.xyz -v | tee salida

