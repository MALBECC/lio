#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi
SALIDALIBXC=salida.libxc
if [ -n "$1" ]
    then
    SALIDALIBXC=$1
fi

$LIOBIN -i agua.in -b basis -c agua.xyz -v > $SALIDA
$LIOBIN -i agua.libxc.in -b basis -c agua.xyz -v > $SALIDALIBXC
