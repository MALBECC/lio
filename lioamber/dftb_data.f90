!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_data
   implicit none
   

   logical      :: dftb_calc =.false. !logical indicator for dftb calculation
   integer      :: MTB   = 0 ! Size of the two tight-binding subatrices
   integer      :: MDFTB = 0 ! Size of the DFT-TB matrix

   integer      :: end_basis              ! Index matrix size
   integer, allocatable :: Iend(:,:) ! Index matrix


   real*8       :: alfaTB  ! Diagonal tight binding param (fermi energy)
   real*8       :: betaTB  ! Offdiagonal tight binding param
   real*8       :: gammaTB ! DFT-TB terms
   real*8       :: Vbias   ! Bias potential
      
!#ifdef TD_SIMPLE
   real*8,allocatable   :: chimerafock (:, :) !allocate in the central code
!#else
!   real*16, allocatable :: chimerafock (:, :)
!#endif

end module dftb_data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
