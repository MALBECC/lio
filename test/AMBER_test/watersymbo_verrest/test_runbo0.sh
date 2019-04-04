#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="watersymbo"

   prmtop="${system}_xx0.prmtop"
   inpcrd="${system}_xx0.rst7"
   amberi="${system}_bo0.in"

   stdout="${system}_bo0.stdo"
   ambero="${system}_bo0.out"
   rstrto="${system}_bo0.rst7"
   mdcrdo="${system}_bo0.mdcrd"
   mdvelo="${system}_bo0.mdvel"
   mdfrzo="${system}_bo0.mdfrz"

   lioinp="liobo0.in"
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
