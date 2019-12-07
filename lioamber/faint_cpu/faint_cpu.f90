!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% MODULE FAINT_CPU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This modules contain all the CPU KS analytic integrals (i.e., all non-XC     !
! integrals) and calculations for the density fitting sets, using the          !
! Obara-Saika recursive method.                                                !
! Edit at your own risk, it took two months to properly reformat this. Think   !
! twice about wasting your time here unless it is strictly necessary.          !
!                                                                              !
! Original integrals: Dario Estrin, 1992                                       !
! Original module:    Francisco Ramirez, Oct/2017                              !
! First refactor:     Diego Armiño, May/2018                                   !
! Full refactor:      Federico Pedron, Sep/2018                                !
! Included ECP:       Nicolas Foglia, Nov/2019                                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module faint_cpu
   ! KS Integrals.
   use subm_int1   , only: int1
   use subm_int2   , only: int2
   use subm_int3lu , only: int3lu
   use subm_int3mem, only: int3mem
   use subm_intSol , only: intSol
   use subm_intECP , only: intECP

   ! Gradients.
   use subm_int1G  , only: int1G
   use subm_int2G  , only: int2G
   use subm_int3G  , only: int3G
   use subm_intSG  , only: intSG
   use subm_intSolG, only: intSolG
   use subm_intECPG, only: intECPG, ECP_gradients

   ! Other integrals.
   use subm_intfld , only: intfld
   implicit none
end module faint_cpu
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
