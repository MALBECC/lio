#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="watersymbo"

   prmtop="${system}_xx.prmtop"
   inpcrd="${system}_xx.rst7"
   amberi="${system}_eh.in"

   stdout="${system}_eh.stdo"
   ambero="${system}_eh.out"
   rstrto="${system}_eh.rst7"
   mdcrdo="${system}_eh.mdcrd"
   mdvelo="${system}_eh.mdvel"
   mdfrzo="${system}_eh.mdfrz"

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
#------------------------------------------------------------------------------#
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
