#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="watersymbo"

   prmtop="${system}_xx0.prmtop"
   inpcrd="${system}_xx0.rst7"
   amberi="${system}_eh1.in"

   stdout="${system}_eh1.stdo"
   ambero="${system}_eh1.out"
   rstrto="${system}_eh1.rst7"
   mdcrdo="${system}_eh1.mdcrd"
   mdvelo="${system}_eh1.mdvel"
   mdfrzo="${system}_eh1.mdfrz"

   lioinp="lioeh1.in"
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
