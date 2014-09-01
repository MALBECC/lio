#!/bin/bash
runid=$(date +%Y-%m-%d_%H-%M-%S)
dir=$(readlink -f ../performance-measures/hemo-$runid)
pushd ../../g2g
make clean && make -j cpu=1 time=1
cd ..
make cpu=1 time=1
popd
mkdir -p $dir
amplxe-cl -collect advanced-hotspots -result-dir $dir -- ../../liosolo/liosolo -i hemo.profile.in -b DZVP -c hem.xyz > $dir/output.txt
func=$(amplxe-cl -report hotspots -result-dir $dir 2>&1 | tail -n +5 | head -n 1 | cut -d' ' -f1)
amplxe-cl -report hotspots -source-object function="$func" -result-dir $dir > $dir/rundown.txt
amplxe-cl -report summary -result-dir $dir > $dir/summary.txt
