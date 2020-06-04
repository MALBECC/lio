#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="watersymbo"

   prmtop="${system}_xx0.prmtop"
   inpcrd="${system}_xx0.rst7"
   amberi="${system}_eh0.in"

   stdout="${system}_eh0.stdo"
   ambero="${system}_eh0.out"
   rstrto="${system}_eh0.rst7"
   mdcrdo="${system}_eh0.mdcrd"
   mdvelo="${system}_eh0.mdvel"
   mdfrzo="${system}_eh0.mdfrz"
}
#------------------------------------------------------------------------------#
RUNAMBER ()
{
   sander="${AMBERHOME}/bin/sander -O"
   cmdinp="-p ${prmtop} -c ${inpcrd} -i ${amberi}"
   cmdout="-o ${ambero} -r ${rstrto} -x ${mdcrdo}"
   cmdout=${cmdout}"  -frc ${mdfrzo} -v ${mdvelo}"

   ${sander} ${cmdinp} ${cmdout} > ${stdout}
}
################################################################################
# RUN ORDER
################################################################################
RUNSETUP
RUNAMBER
################################################################################
