#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida.libxc
if [ -n "$1" ]
  then
    SALIDA=$1
fi

echo $LIOBIN -i agua.libxc.in -b basis -c agua.xyz -v

$LIOBIN -i agua.libxc.in -b basis -c agua.xyz -v > $SALIDA



