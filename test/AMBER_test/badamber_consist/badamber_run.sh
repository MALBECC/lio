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

}
#------------------------------------------------------------------------------#
RUNAMBER ()
{
   sander="${AMBERHOME}/bin/sander -O"
   cmdinp="-p ${prmtop} -c ${inpcrd} -i ${amberi}"
   cmdout="-o ${ambero} -r ${rstrto} -x ${mdcrdo}"
   cmdout=${cmdout}"  -frc ${mdfrzo} -v ${mdvelo}"

   ${sander} ${cmdinp} ${cmdout} > ${stdout}
   rm mdinfo
#   rm qm.xyz mdinfo
}
################################################################################
# RUN ORDER
################################################################################
RUNSETUP "waterlong"
RUNAMBER
RUNSETUP "watershrt"
RUNAMBER
################################################################################
