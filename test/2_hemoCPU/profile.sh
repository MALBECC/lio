#!/bin/bash
runid=$(date +%Y-%m-%d_%H-%M-%S)
dir=$(readlink -f ../performance-measures/hemo-$runid)
amplxe-cl -collect advanced-hotspots -result-dir $dir -- ../../liosolo/liosolo -i hemo.profile.in -b DZVP -c hem.xyz
func=$(amplxe-cl -report hotspots -result-dir $dir 2>&1 | tail -n +5 | head -n 1 | cut -d' ' -f1)
echo "Top consuming functions is $func"
echo "Run down for top consuming function: "
echo " "
amplxe-cl -report hotspots -source-object function="$func" -result-dir $dir > $dir/rundown.txt
