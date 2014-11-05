#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

$LIOBIN -i hemo.profile.in -b DZVP  -c hem.xyz -v 
