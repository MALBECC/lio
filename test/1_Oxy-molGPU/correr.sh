#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
$LIOBIN -i Oxy-mol.in -b DZVP  -c Oxy-mol.xyz -v > salida

