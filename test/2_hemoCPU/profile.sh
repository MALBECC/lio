#!/bin/bash
analysis=${1:-hotspots}
runid=$(date +%Y-%m-%d_%H-%M-%S)
dir=$(readlink -f ../performance-measures/hemo-$runid-$analysis)
pushd ../../g2g
make clean && make -j cpu=1 time=1 openmp=1
cd ..
rm liosolo/liosolo
make cpu=1 time=1
popd
mkdir -p $dir
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
amplxe-cl -collect $analysis -result-dir $dir -- ../../liosolo/liosolo -i hemo.profile.in -b DZVP -c hem.xyz > $dir/output.txt
