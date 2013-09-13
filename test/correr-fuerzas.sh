#!/bin/bash
set -e
rm -f forces.out opt.xyz
sed -ri 's|max_function_exponent.*$|max_function_exponent 8|' gpu_options
cd ../../g2g; ln -sf libg2g-cpu-f.so libg2g.so; cd -
../../fortran/garcha-g2g < tf | tee -i cpu-single.out; mv forces.out forces-cpu-single-8.out; mv opt.xyz cpu-single-8.xyz
cd ../../g2g; ln -sf libg2g-gpu-f.so libg2g.so; cd -
../../fortran/garcha-g2g < tf | tee -i gpu-single.out; mv forces.out forces-gpu-single-8.out; mv opt.xyz gpu-single-8.xyz
cd ../../g2g; ln -sf libg2g-gpu-d.so libg2g.so; cd -
../../fortran/garcha-g2g < tf | tee -i gpu-double.out; mv forces.out forces-gpu-double-8.out; mv opt.xyz gpu-double-8.xyz
cd ../../g2g; ln -sf libg2g-cpu-d.so libg2g.so; cd -
../../fortran/garcha-g2g < tf | tee -i cpu-double.out; mv forces.out forces-cpu-double-8.out; mv opt.xyz cpu-double-8.xyz

sed -ri 's|max_function_exponent.*$|max_function_exponent 15|' gpu_options
cd ../../g2g; ln -sf libg2g-gpu-d.so libg2g.so; cd -
../../fortran/garcha-g2g < tf | tee -i gpu-double.out; mv forces.out forces-gpu-double-15.out; mv opt.xyz gpu-double-15.xyz
cd ../../g2g; ln -sf libg2g-cpu-d.so libg2g.so; cd -
../../fortran/garcha-g2g < tf | tee -i cpu-double.out; mv forces.out forces-cpu-double-15.out; mv opt.xyz cpu-double-15.xyz

