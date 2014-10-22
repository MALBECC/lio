#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=/home/lab8/nano/lio/liosolo/liosolo
fi

$LIOBIN -i fullereno.in -b DZVP  -c fullereno.xyz -v > salida

