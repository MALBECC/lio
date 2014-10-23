#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../lio/liosolo
fi


export OMP_NUM_THREADS=240
export MKL_NUM_THREADS=240
export KMP_AFFINITY=explicit,granularity=fine,proclist=[1-239:1,0]
export OMP_NESTED=false
export MKL_DYNAMIC=false
export LIO_OUTER_THREADS=1
export LIO_INNER_THREADS=240
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/xeonphifs/lio
ulimit -s unlimited

$LIOBIN -i hemo.in -b DZVP  -c hem.xyz -v | tee salida
