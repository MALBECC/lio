#!/bin/bash
export SINK_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
ssh root@mic0 "rm -rf /tmp/lio"
ssh root@mic0 "mkdir -p /tmp/lio"
micnativeloadex liosolo/liosolo -l | python deps-xeon.py | while read dep
do
    scp $dep root@mic0:/tmp/lio
done
scp liosolo/liosolo root@mic0:/tmp/lio
