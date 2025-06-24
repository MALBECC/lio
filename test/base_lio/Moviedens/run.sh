#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=output
if [ -n "$1" ]
  then
    SALIDA=$1
fi

echo "This test will take 30 min to run in CPU."

FIXINP ()
{
    fxyz_inp=${1}
    fxyz_out=${2}
    tail --lines=+3 ${fxyz_inp} > ${fxyz_out}
    sed -i 's/H / 1/g' ${fxyz_out}
    sed -i 's/C / 6/g' ${fxyz_out}
    sed -i 's/N / 7/g' ${fxyz_out}
    sed -i 's/O / 8/g' ${fxyz_out}
    sed -i 's/S /16/g' ${fxyz_out}
    sed -i 's/Fe/26/g' ${fxyz_out}
    sed -i 's/Cu/29/g' ${fxyz_out}
}

source ../../../liohome.sh

FIXINP chain_test.xyz temp_coord.in
cp restart_grnd.in restart.in
$LIOBIN -i chain_grnd.in -c temp_coord.in -v >  chain_grnd.out
cp liomovie_el0000.out dens_neg.in

FIXINP chain_test.xyz temp_coord.in
cp restart_edyn.in restart.in
$LIOBIN -i chain_edyn.in -c temp_coord.in -v > chain_edyn.out

RUN_1PIC ()
{
    IDN=${1}
    FIXINP liomovie_nu${IDN}.out temp_coord.in
    cp liomovie_el${IDN}.out dens_pos.in
    ../../../tools/densdif.py
    mv dens_dif.out restart.in
    $LIOBIN -i chain_1pic.in -c temp_coord.in -v > chain_1pic.out
    mv chain_1pic.cube frame${IDN}.cube
    rm dens_pos.in
}

RUN_1PIC 0000
RUN_1PIC 0001
RUN_1PIC 0002
RUN_1PIC 0003
RUN_1PIC 0004
RUN_1PIC 0005
RUN_1PIC 0006
RUN_1PIC 0007
RUN_1PIC 0008
RUN_1PIC 0009
RUN_1PIC 0010

