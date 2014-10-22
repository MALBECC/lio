#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=1
export OMP_NESTED=true
export MKL_DYNAMIC=true
export KMP_AFFINITY=granularity=fine,scatter
export LIO_OUTER_THREADS=12
export LIO_INNER_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -f ../../g2g):$(readlink -f ../../lioamber)
$LIOBIN -i hemo.in -b DZVP  -c hem.xyz -v | tee salida
