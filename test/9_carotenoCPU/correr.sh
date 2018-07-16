#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

SALIDA_LIBXC=salida.libxc
if [ -n "$1" ]
  then
    SALIDA_LIBXC=$1
fi

export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=1
export KMP_AFFINITY=granularity=core,scatter
export LIO_OUTER_THREADS=12
export LIO_INNER_THREADS=1

$LIOBIN -i carotenox.in -b DZVP -c caroteno.xyz -v > $SALIDA
$LIOBIN -i carotenox.libxc.in -b DZVP -c caroteno.xyz -v > $SALIDA_LIBXC
