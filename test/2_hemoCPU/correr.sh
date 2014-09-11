#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1
$LIOBIN -i hemo.in -b DZVP  -c hem.xyz -v | tee salida
