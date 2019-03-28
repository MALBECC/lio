#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

$LIOBIN -i fe3h2o.in -b DZVP  -c fe3h2o.xyz -v > $SALIDA

