#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="watersymbo"

   prmtop="${system}_xx.prmtop"
   inpcrd="${system}_xx.rst7"
   amberi="${system}_bo.in"

   stdout="${system}_bo.stdo"
   ambero="${system}_bo.out"
   rstrto="${system}_bo.rst7"
   mdcrdo="${system}_bo.mdcrd"
   mdvelo="${system}_bo.mdvel"
   mdfrzo="${system}_bo.mdfrz"

   basisf="basis_dzvp.in"
   lioinp="liobo.in"
}
#------------------------------------------------------------------------------#
RUNAMBER ()
{
   sander="${AMBERHOME}/bin/sander -O"
   cmdinp="-p ${prmtop} -c ${inpcrd} -i ${amberi}"
   cmdout="-o ${ambero} -r ${rstrto} -x ${mdcrdo}"
   cmdout=${cmdout}"  -frc ${mdfrzo} -v ${mdvelo}"

   cp ${basisf} basis
   cp ${lioinp} lio.in
   ${sander} ${cmdinp} ${cmdout} > ${stdout}
   rm basis lio.in
   rm qm.xyz mdinfo fort.*
}
COMPARES ()
{
   echo "not ready yet"
}
################################################################################
# RUN ORDER
################################################################################
RUNSETUP
RUNAMBER
COMPARES
################################################################################
