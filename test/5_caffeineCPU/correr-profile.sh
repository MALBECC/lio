#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

$LIOBIN -i caff.profile.in -b DZVP  -c caff.xyz -v 
