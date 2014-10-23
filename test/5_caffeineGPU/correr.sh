#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
$LIOBIN -i caff.in -b DZVP  -c caff.xyz -v > salida

