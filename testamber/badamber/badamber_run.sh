#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="$1"

   prmtop="waterboth.prmtop"
   inpcrd="waterboth.rst7"
   amberi="${system}.in"

   stdout="${system}.stdo"
   ambero="${system}.out"
   rstrto="${system}.rst7"
   mdcrdo="${system}.mdcrd"
   mdvelo="${system}.mdvel"
   mdfrzo="${system}.mdfrz"

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
################################################################################
# RUN ORDER
################################################################################
RUNSETUP "waterlong"
RUNAMBER
RUNSETUP "watershrt"
RUNAMBER
################################################################################
