#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

unbuffer $LIOBIN -i chloride.in -c chloride.xyz -v > $SALIDA

if [ -z "$TDAN" ] ; then
  TDAN=../../lioamber/tdanalize/tdanalyze
fi

$TDAN dipole_moment_td 100 900 0 3
