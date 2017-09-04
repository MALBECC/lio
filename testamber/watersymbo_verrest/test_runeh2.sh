#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="watersymbo"

   prmtop="${system}_xx0.prmtop"
   inpcrd="${system}_eh1.rst7"
   amberi="${system}_eh2.in"

   stdout="${system}_eh2.stdo"
   ambero="${system}_eh2.out"
   rstrto="${system}_eh2.rst7"
   mdcrdo="${system}_eh2.mdcrd"
   mdvelo="${system}_eh2.mdvel"
   mdfrzo="${system}_eh2.mdfrz"

   lioinp="lioeh2.in"
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
