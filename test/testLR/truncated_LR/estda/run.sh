#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

export GFORTRAN_UNBUFFERED_ALL=1

export LIOHOME=/home/gonzalo/progs/LIOs/codes/LRguess
#export LIOHOME=/home/gonzalo/progs/LIOs/codes/LRguess_antes
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOHOME/g2g:$LIOHOME/lioamber
${LIOHOME}/liosolo/liosolo -i co.in -c co.xyz -v > $SALIDA

