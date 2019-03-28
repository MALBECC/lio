#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

export LIOHOME=/opt/progs/LIOs/NewTest
$LIOBIN -i cl2.in -ib -bs 'SBKJC' -fs 'DZVP Coulomb Fitting' -c cl2.xyz -v > $SALIDA


