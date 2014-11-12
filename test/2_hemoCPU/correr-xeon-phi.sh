#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../lio/liosolo
fi

ulimit -s unlimited
export OMP_NUM_THREADS=60
export MKL_NUM_THREADS=60
export KMP_AFFINITY=explicit,granularity=fine,proclist=[1-239:1,0]
export OMP_NESTED=false
export MKL_DYNAMIC=false
export LIO_MINCOST_OFFSET=100000
export LIO_SPLIT_THRESHOLD=20
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/xeonphifs/lio

$LIOBIN -i hemo.in -b DZVP  -c hem.xyz -v | tee salida
