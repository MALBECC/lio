#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

#$LIOBIN -i cl2.in -b DZVP  -c cl2.xyz -v > $SALIDA

$LIOBIN -i cl2.in -ib -bs 'SBKJC' -fs 'DZVP Coulomb Fitting' -b basis -c cl2.xyz -v > $SALIDA


