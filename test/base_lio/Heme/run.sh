#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

$LIOBIN -i heme.in -b DZVP -c heme.xyz -v > $SALIDA

