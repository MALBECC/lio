#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

export OMP_NUM_THREADS=$1
export OMP_NESTED=false
export LIO_OPTIONS_FILE=gpu_options.small
export KMP_AFFINITY=granularity=fine,scatter
export LIO_SPLIT_THRESHOLD=80
export LIO_MINCOST_OFFSET=250000
export LIO_OUTER_THREADS=$1
export LIO_INNER_THREADS=$1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -f ../lio-g2gs/latest)

$LIOBIN -i hemo.in -b DZVP  -c hem.xyz -v | tee salida.txt
