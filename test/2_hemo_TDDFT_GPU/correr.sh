#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

$LIOBIN -i hemo.in -b DZVP  -c hem.xyz -v > $SALIDA

