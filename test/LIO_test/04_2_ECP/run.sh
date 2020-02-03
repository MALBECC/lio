#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

source ../../../liohome.sh
$LIOBIN -i fe3h2o.in -c fe3h2o.xyz -v > $SALIDA

