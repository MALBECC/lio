#!/bin/bash
cd g2g && make clean && cd ..
rm -rf liosolo/liosolo
make -j cpu=1 time=1 openmp=1 cpu_recompute=0 > >(tee stdout.log) 2> >(tee stderr.log >&2)
