#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

$LIOBIN -i Oxy-mol.in -c Oxy-mol.xyz -v > $SALIDA

