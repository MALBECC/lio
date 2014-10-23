#!/bin/bash
export SINK_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
rm -rf /xeonphifs/lio
mkdir -p /xeonphifs/lio
micnativeloadex liosolo/liosolo -l | python deps-xeon.py | while read dep
do
    cp $dep /xeonphifs/lio
done
cp liosolo/liosolo /xeonphifs/lio

rm -rf /xeonphifs/hemo
mkdir -p /xeonphifs/hemo
cp -r test/2_hemoCPU/* /xeonphifs/hemo
