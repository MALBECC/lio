#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

$LIOBIN -i fullereno.in -b DZVP  -c fullereno.xyz -v > salida

