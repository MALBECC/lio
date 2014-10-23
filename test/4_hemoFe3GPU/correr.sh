#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
$LIOBIN -i hemo.in -b DZVP  -c hemo.xyz -v > salida

