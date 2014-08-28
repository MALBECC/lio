#!/bin/bash
ssh root@mic0 "mkdir -p /tmp/lio"
micnativeloadex liosolo/liosolo -l | python depts.py | while read dep
do
    scp $dep root@mic0:/tmp/lio
done
scp liosolo/liosolo root@mic0:/tmp/lio
