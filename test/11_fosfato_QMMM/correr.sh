#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=/home/nano/lio/liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

$LIOBIN -i fos.in  -c  fosqmmm.xyz -b basis  -v > $SALIDA

