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

echo testing using lio
$LIOBIN -i caff.in -b DZVP  -c caff.xyz -v > $SALIDA
$LIOBIN -i caff.libxc.in -b DZVP -c caff.xyz -v > $SALIDALIBXC

