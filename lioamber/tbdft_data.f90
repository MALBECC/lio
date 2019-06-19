#include "complex_type.fh"
module tbdft_data
   implicit none

   logical      :: tbdft_calc =.false.      ! Logical indicator for tbdft calculation
   integer      :: MTB   = 0                ! Size of the two tight-binding subatrices
   integer      :: MTBDFT = 0               ! Size of the DFT-TB matrix
   integer      :: start_tdtb=0             ! Initial time step for evolution of diagonal TB terms
   integer      :: end_tdtb=0               ! Final time step for evolution of diagonal TB terms
   integer      :: end_bTB                  ! Index matrix size
   real(kind=8) :: alfaTB                   ! Fermi Energy
   real(kind=8) :: betaTB                   ! Offdiagonal tight binding param
   real(kind=8) :: gammaTB                  ! DFT-TB terms
   real(kind=8) :: Vbias_TB                 ! Bias potential
   real(kind=8) :: chargeA_TB
   real(kind=8) :: chargeB_TB
   real(kind=8) :: chargeM_TB
   
   integer     , allocatable :: Iend_TB(:,:)        ! Index matrix
   real(kind=8), allocatable :: rhoa_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8), allocatable :: rhob_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8), allocatable :: chimerafock (:,:,:) ! Allocated in the central code
   real(kind=8), allocatable :: gammaW(:)           ! gamma weight
   TDCOMPLEX   , allocatable :: rhold_AOTB(:,:,:)   ! rho in AO to calculate charges
   TDCOMPLEX   , allocatable :: rhonew_AOTB(:,:,:)  ! rho in AO to calculate charges
end module tbdft_data
