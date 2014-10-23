#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
$LIOBIN -i fe3h2o.in -b DZVP  -c fe3h2o.xyz -v > salida

