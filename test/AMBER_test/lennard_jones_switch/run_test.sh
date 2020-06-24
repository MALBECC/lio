#/bin/bash
################################################################################
# CONFIGURABLE CODE
################################################################################
RUNSETUP ()
{
   system="$1"

   prmtop="NH4.prmtop"
   inpcrd="NH4.rst7"
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
   rm mdinfo ${mdcrdo} ${mdvelo} ${rstrto} ${stdout}
}
################################################################################
# RUN ORDER
################################################################################
RUNSETUP "lj_null"
RUNAMBER
RUNSETUP "lj_same"
RUNAMBER
RUNSETUP "lj_same_border"
RUNAMBER
RUNSETUP "lj_same_mid"
RUNAMBER
################################################################################
