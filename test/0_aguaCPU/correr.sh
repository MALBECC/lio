#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi
SALIDALIO=salida.lio
if [ -n "$1" ]
    then
    SALIDALIO=$1
fi

$LIOBIN -i agua.in -b basis -c agua.xyz -v > $SALIDA
$LIOBIN -i agua.libxc.in basis -c agua.xyz -v > $SALIDALIO