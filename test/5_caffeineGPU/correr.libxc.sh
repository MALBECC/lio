#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida.libxc
if [ -n "$1" ]
  then
    SALIDA=$1
fi

$LIOBIN -i caff.libxc.in -b DZVP  -c caff.xyz -v > $SALIDA

