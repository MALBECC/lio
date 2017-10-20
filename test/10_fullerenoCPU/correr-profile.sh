#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

$LIOBIN -i fullereno.profile.in -b DZVP  -c fullereno.xyz -v 
