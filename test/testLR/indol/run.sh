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
export OMP_NUM_THREADS=4
#export LIOHOME=/home/gonzalo/progs/LIOs/code-to-develop/method_old
export LIOHOME=/home/gonzalo/progs/LIOs/code-to-develop/develop-excited2
#export LIO_SPLIT_POINTS=5000
$LIOHOME/liosolo/liosolo -i indol.in -c indol.xyz -v > $SALIDA
