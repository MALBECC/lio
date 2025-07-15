#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

export OMP_NUM_THREADS=1
export GFORTRAN_UNBUFFERED_ALL=1

source ../../../liohome.sh
$LIOBIN -i Oxy-mol.in -c Oxy-mol.xyz -v > $SALIDA

