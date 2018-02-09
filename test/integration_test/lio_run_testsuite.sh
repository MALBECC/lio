#!/usr/bin/env bash

if [ -n "$SKIP_CHECK" ]; then
  echo "Skipping checks"
  exit 0
fi

if [ -z "$destdir" ]; then
  destdir="${HOME}/workspace/tesis/experimentos/lio"
fi

# Write the results according to the running mode
# i.e.: lio_cpu_libx_cpu
#	lio_gpu_libxc_cpu
#	lio_gpu_libxc_gpu
# Destination folders for the results
if [ -z "$liocpulibxccpu" ]; then
  liocpulibxccpufolder="/lio_cpu_libxc_cpu"
fi
if [ -z "$liogpulibxccpufolder" ]; then
  liogpulibxccpufolder="/lio_gpu_libxc_cpu"
fi
if [ -z "$liogpulibxcgpufolder" ]; then
  liogpulibxcgpufolder="/lio_gpu_libxc_gpu"
fi

#TODO: determine the running mode
if [ -z "$runningmode" ]; then
  runningmode="/lio_gpu_libxc_cpu"
fi

if [ -z "$builddir" ]; then
 builddir="./"
fi

# Executable file
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

# Current test folder
if [ -z "$currtestdir" ]; then
  currtestdir="../"
fi

# We only print using colors if the output is going to a terminal
if [ -t 1 ] ; then
  NC='\033[0m'
  RED='\033[1;31m'
  GREEN='\033[1;32m'
  YELLOW='\033[1;33m'
  PINK='\033[1;35m'
  CYAN='\033[1;36m'
fi

#The directory to save the results
workdir=${destdir}

fail=0
echo -e "${YELLOW}Test suite runner for lio ${NC}"
echo -e "Using ${PINK}${workdir}${NC} as the results directory"
echo -e "========================================================================="
#echo -e "   Functional                    System      NSpin   Quantity    Result"
echo -e "========================================================================="

# Timestamp for the test folder name
now=$(date +"%Y%m%d%T")

# Run the test one by one (it may take a while)

## AGUA
## LIO
currenttest="agua"
inputfile="${currenttest}"
outputfile=salida.lio
outputdir="${destdir}${runningmode}/${now}/${currenttest}/"

echo -e "Creating dir ${outputdir} if not exists"
mkdir -p ${outputdir}

echo -e "${CYAN}Test name - ${currenttest}${NC}"
#Show the command
#echo $LIOBIN -i agua.libxc.in -b basis -c agua.xyz -v
echo -e "${GREEN}Running LIO ${NC} $LIOBIN -i ${inputfile}.in -b ${currentest}.basis -c ${currenttest}.xyz -v  > ${outputdir}${outputfile}"
#Run the test
$LIOBIN -i ${inputfile}.in -b ${inputfile}.basis -c ${inputfile}.xyz -v > ${outputdir}${outputfile};

## LIBXC
inputfile="${currenttest}.libxc"
outputfile=salida.libxc

#Show the command
echo -e "${GREEN}Running LIBXC ${NC} $LIOBIN -i ${inputfile}.in -b ${currenttest}.basis -c ${currenttest}.xyz -v  > ${outputdir}${outputfile}"
#Run the test
$LIOBIN -i ${inputfile}.in -b ${currenttest}.basis -c ${currenttest}.xyz -v > ${outputdir}${outputfile};

echo -e "========================================================================="
echo -e "========================================================================="

##TODO AGREGAR: caroteno, hemo que son los que puedo comparar timers.
## CAFFEINE
## LIO
currenttest="caffeine"
inputfile="${currenttest}"
outputfile=salida.lio
outputdir="${destdir}${runningmode}/${now}/${currenttest}/"

echo -e "Creating dir ${outputdir} if not exists"
mkdir -p ${outputdir}

echo -e "${CYAN}Test name - ${currenttest}${NC}"
#Show the command
echo -e "${GREEN}Running LIO ${NC} $LIOBIN -i ${inputfile}.in -b ${currentest}.basis -c ${currenttest}.xyz -v  > ${outputdir}${outputfile}"
#Run the test
#$LIOBIN -i ${inputfile}.in -b ${inputfile}.basis -c ${inputfile}.xyz -v > ${outputdir}${outputfile};

## LIBXC
inputfile="${currenttest}.libxc"
outputfile=salida.libxc

#Show the command
echo -e "${GREEN}Running LIBXC ${NC} $LIOBIN -i ${inputfile}.in -b ${currenttest}.basis -c ${currenttest}.xyz -v  > ${outputdir}${outputfile}"
#Run the test
#$LIOBIN -i ${inputfile}.in -b ${currenttest}.basis -c ${currenttest}.xyz -v > ${outputdir}${outputfile};

echo -e "========================================================================="
echo -e "========================================================================="

##TODO AGREGAR: caroteno, hemo que son los que puedo comparar timers.
#$LIOBIN -i fullereno.in -b DZVP  -c fullereno.xyz -v > $SALIDA


echo -e "${NC}"
exit $fail




# Esto de abajo se borra
#for dir in $(ls -d $srcdir/regression/*); do
#  for i in $(ls $dir/*.bz2 | LANG=C sort); do
#    refname=$(basename ${i} .bz2)
#    func=${refname%%.*}
#    testcase=${refname#*.}
#    system=${testcase%%.*}
#    order=${testcase##*.}
#    if [ "x${testcase#*.}" = "xunpol.$order" ]; then
#      nspin=1
#    else
#      nspin=2
#    fi
#    case $order in
#    0)
#      label=exc
#      tol=5e-8
#      ;;
#    1)
#      label=vxc
#      tol=5e-5
#      ;;
#    2)
#      label=fxc
#      tol=5e-4
#      ;;
#    esac
#    name=$(printf '%-30s %-12s %-8d %-5s' ${func} ${system} ${nspin} ${label})
#    echo -ne "${NC} * ${PINK}${name}"

#    bunzip2 -c $dir/$refname.bz2 > $workdir/${refname}_ref

#    $EXEC $builddir/xc-regression $func $nspin $order $srcdir/input/$system $workdir/$refname > /dev/null

#    res=$($builddir/xc-error $workdir/$refname $workdir/${refname}_ref $tol)
#    if [ "x$res" = "x1" ]; then
#      echo -e "      ${GREEN}OK"
#    else
#      echo -e "     ${RED}FAIL"
#      let fail++
#    fi
#  done
#done
echo -e "${NC}"
exit $fail
