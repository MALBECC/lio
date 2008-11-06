#!/bin/bash

primera_grilla=1

for i in agu12; do
  cd ${i}
  
  echo "${i} (con grilla ${primera_grilla})"
  sed -ri "s|IGRID2\s*=\s*.|IGRID2 = ${primera_grilla}|" ${i}.i
 
  echo 'Single-point GPU'
  echo ${i}.i > t
  /usr/bin/time -p ../../garcha/garcha-gpu < t &> ${i}-igrid${primera_grilla}.g

  echo 'Single-point CPU'
  echo ${i}.i > t
  /usr/bin/time -p ../../garcha/garcha-cpu < t &> ${i}-igrid${primera_grilla}.c

  #echo 'Optimizacion GPU Normal'
  #echo ${i}f.i > t
  #/usr/bin/time -p ../../../garcha/garcha-gpu < t &> ${i}-f.g
  #mv opt.xyz ${i}-gpu.xyz

  #echo 'Optimizacion CPU Normal'
  #echo ${i}f.i > t
  #/usr/bin/time -p ../../../garcha/garcha-cpu < t &> ${i}-f.c
  #mv opt.xyz ${i}-cpu.xyz
  
  #echo "${i} (con grilla ${primera_grilla})"
  #sed -ri "s|IGRID2\s*=\s*.|IGRID2 = ${primera_grilla}|" ${i}f.i  

  #echo "${i} Optimizacion GPU con maximo de 10 geometrias"
  #echo ${i}f.i > t
  #/usr/bin/time -p ../../../garcha/garcha-gpu-10 < t &> ${i}-f10.g
  #mv opt.xyz ${i}-10-gpu.xyz

  #echo "${i} Optimizacion CPU con maximo de 10 geometrias"
  #echo ${i}f.i > t
  #/usr/bin/time -p ../../../garcha/garcha-cpu-10 < t &> ${i}-f10.c
  #mv opt.xyz ${i}-10-cpu.xyz

  cd ..
done
