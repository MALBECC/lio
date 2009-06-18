#!/bin/bash

primera_grilla=1

for i in $*; do
  cd ${i}
  
  echo "==> $i"
  /usr/bin/time -p ../../garcha/garcha-cpu < t &> cpu.log

  cd ..
done
