#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

export GFORTRAN_UNBUFFERED_ALL=1
source ../../../liohome.sh
${LIOBIN} -i indol.in -c indol.xyz -v > $SALIDA

#export LIOHOME=/home/gonzalo/progs/LIOs/develop-orig
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOHOME/g2g:$LIOHOME/lioamber
#${LIOHOME}/liosolo/liosolo -i indol.in -c indol.xyz -v > $SALIDA


