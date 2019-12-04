#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

source ../../../liohome.sh
export GFORTRAN_UNBUFFERED_ALL=1
$LIOBIN -i agua.in -b basis -c agua.xyz -v > $SALIDA
