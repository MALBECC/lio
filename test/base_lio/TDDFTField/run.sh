#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../build/liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

#source ../../../liohome.sh
export LIOHOME=${HOME}/progs/lio-GPU_deprecated
export LD_LIBRARY_PATH=$LIOHOME/lib:${LD_LIBRARY_PATH}
$LIOBIN -i chloride.in -c chloride.xyz -v > $SALIDA


