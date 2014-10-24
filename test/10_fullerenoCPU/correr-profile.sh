#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(readlink -f ../../g2g):$(readlink -f ../../lioamber)
$LIOBIN -i fullereno.profile.in -b DZVP  -c fullereno.xyz -v 
