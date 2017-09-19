#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNAMBER ()
{
   sysname="$1"

   prmtop="water_forall.prmtop"
   inpcrd="water_forall.rst7"
   amberi="${sysname}.in"

   stdout="${sysname}.stdo"
   ambero="${sysname}.out"
   rstrto="${sysname}.rst7"
   mdcrdo="${sysname}.mdcrd"
   mdvelo="${sysname}.mdvel"
   mdfrzo="${sysname}.mdfrz"

   sander="${AMBERHOME}/bin/sander -O"
   cmdinp="-p ${prmtop} -c ${inpcrd} -i ${amberi}"
   cmdout="-o ${ambero} -r ${rstrto} -x ${mdcrdo}"
   cmdout=${cmdout}"  -frc ${mdfrzo} -v ${mdvelo}"

   ${sander} ${cmdinp} ${cmdout} > ${stdout}
   rm mdinfo
}
################################################################################
# RUN ORDER
################################################################################
RUNAMBER "water_scratch"
RUNAMBER "water_restart"
################################################################################
