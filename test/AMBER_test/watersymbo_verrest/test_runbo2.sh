#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="watersymbo"

   prmtop="${system}_xx0.prmtop"
   inpcrd="${system}_bo1.rst7"
   amberi="${system}_bo2.in"

   stdout="${system}_bo2.stdo"
   ambero="${system}_bo2.out"
   rstrto="${system}_bo2.rst7"
   mdcrdo="${system}_bo2.mdcrd"
   mdvelo="${system}_bo2.mdvel"
   mdfrzo="${system}_bo2.mdfrz"

   lioinp="liobo2.in"
}
#------------------------------------------------------------------------------#
RUNAMBER ()
{
   sander="${AMBERHOME}/bin/sander -O"
   cmdinp="-p ${prmtop} -c ${inpcrd} -i ${amberi}"
   cmdout="-o ${ambero} -r ${rstrto} -x ${mdcrdo}"
   cmdout=${cmdout}"  -frc ${mdfrzo} -v ${mdvelo}"

   cp ${lioinp} lio.in
   ${sander} ${cmdinp} ${cmdout} > ${stdout}
   rm lio.in
}
################################################################################
# RUN ORDER
################################################################################
RUNSETUP
RUNAMBER
################################################################################
