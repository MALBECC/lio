!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module faint_cpu
   ! KS Integrals.
   use subm_int1   , only: int1
   use subm_int2   , only: int2
   use subm_int3lu , only: int3lu
   use subm_int3mem, only: int3mem
   use subm_intSol , only: intSol

   ! Gradients.
   use subm_int1G  , only: int1G
   use subm_intSG  , only: intSG
   use subm_intSolG, only: intSolG

   ! Other integrals.
   use subm_intfld , only: intfld
   implicit none
end module faint_cpu
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
