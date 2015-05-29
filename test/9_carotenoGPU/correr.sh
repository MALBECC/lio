! /bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l gpu=2

export CUDA_VISIBLE_DEVICES=0,1
export LIO_MINCOST_OFFSET=500
export LIO_SPLIT_THRESHOLD=40
export LIO_SPLIT_POINTS=250

export LD_LIBRARY_PATH=/opt/gridengine/lib/linux-x64:/opt/openmpi/lib:/state/partition/apps/composer_xe_2013.0.079/compiler/lib/intel64:/opt/intel/mic/coi/host-linux-release/lib:/opt/intel/mic/myo/lib:/state/partition/apps/composer_xe_2013.0.079/mpirt/lib/intel64:/state/partition/apps/composer_xe_2013.0.079/ipp/../compiler/lib/intel64:/state/partition/apps/composer_xe_2013.0.079/ipp/lib/intel64:/state/partition/apps/composer_xe_2013.0.079/compiler/lib/intel64:/state/partition/apps/composer_xe_2013.0.079/mkl/lib/intel64:/state/partition/apps/composer_xe_2013.0.079/tbb/lib/intel64:/usr/local/cuda/lib64:/usr/local/cuda/lib:/share/apps/Amber14/lib:/home/nano/lio/g2g:/home/nano/lio/lioamber:/share/apps/lib:/opt/python/lib

if [ -z "$LIOBIN" ] ; then
  LIOBIN=/home/nano/lio
fi
$LIOBIN/liosolo/liosolo -i carotenox.in -b DZVP  -c caroteno.xyz -v > salida

