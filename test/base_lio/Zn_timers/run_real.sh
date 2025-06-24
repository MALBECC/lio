#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

$LIOBIN -i multiZn.in -b basis -c multiZn.xyz -v > ${SALIDA}

