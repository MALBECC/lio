# 0.002394810

LD_LIBRARY_PATH=:../../lioamber/:../../g2g/:/opt/intel/composer_xe_2013.0.079/compiler/lib/intel64:/opt/intel/mic/coi/host-linux-release/lib:/opt/intel/mic/myo/lib:/opt/intel/composer_xe_2013.0.079/mpirt/lib/intel64:/opt/intel/composer_xe_2013.0.079/ipp/../compiler/lib/intel64:/opt/intel/composer_xe_2013.0.079/ipp/lib/intel64:/opt/intel/composer_xe_2013.0.079/compiler/lib/intel64:/opt/intel/composer_xe_2013.0.079/mkl/lib/intel64:/opt/intel/composer_xe_2013.0.079/tbb/lib/intel64:

runprog=../../liosolo/liosolo
 echo 'estas son las bibliotecas que voy a usar'
 echo 'CHEQUEARLAS!!!!!!!!!!'
ldd ${runprog}
${runprog} -i carotenox.in -b DZVP  -c caroteno.xyz -v > salida

