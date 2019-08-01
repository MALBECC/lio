#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l gpu=1 #n gpus
#$ -pe mpich 1 #n cpu
#$ -q 1080.q
#$ -N TestOPEN #nombre proceso
rm estado

i=16

export LD_LIBRARY_PATH=/share/apps/libs/:/share/apps/GCC-4.9.4/lib/:/share/apps/GCC-4.9.4/lib64
export GFORTRAN_UNBUFFERED_ALL=1
export LIOHOME=/home/nick/LIO-HYB/Test-OPEN/$i/lio/


#Intel
source  /share/apps/intel/composer_xe_2013.0.079/bin/compilervars.sh intel64
export MKL_HOME=/share/apps/intel/composer_xe_2013.0.079/mkl

#CUDA
export CUDA_HOME=/share/apps/cuda-8.0
export PATH=$PATH:$CUDA_HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib:$CUDA_HOME/lib64/:$CUDA_HOME/lib:

#################
### CUDA LIBRE
export CUDA_VISIBLE_DEVICES=`/share/apps/freegpu.sh`
###################



export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/intel/composer_xe_2013.0.079/compiler/lib/intel64:/share/apps/intel/composer_xe_2013.0.079/mpirt/lib/intel64:/share/apps/intel/composer_xe_2013.0.079/ipp/../compiler/lib/intel64:/share/apps/intel/composer_xe_2013.0.079/ipp/lib/intel64:/share/apps/intel/composer_xe_2013.0.079/compiler/lib/intel64:/share/apps/intel/composer_xe_2013.0.079/mkl/lib/intel64:/share/apps/intel/composer_xe_2013.0.079/tbb/lib/intel64/gcc4.4:/share/apps/cuda-8.0/lib:/share/apps/cuda-8.0/lib64/:$LIOHOME/g2g:$LIOHOME/lioamber:$LIOHOME


##################################################################
echo empiezo $i >> estado

echo "LDlib" $LD_LIBRARY_PATH > vars
echo " " >> vars
echo "liohome" $LIOHOME >> vars
echo " " >> vars
#ldd $LIOHOME/hybrid/bin/hybrid >> vars


ulimit -s unlimited
#ulimit -i 127733
#ulimit -u unlimited 
ulimit -a > utest

$LIOHOME/liosolo/liosolo  -i Oxy-mol.in -b DZVP  -c Oxy-mol.xyz -v > salida

cd ..
echo termine $i >> estado

##################################################################


