!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_data
   implicit none

   integer      :: MTB    ! Size of the two tight-binding subatrices
   integer      :: MTBQM  ! Size of the DFT-TB matrix

   integer      :: nend              ! Index matrix size
   integer, allocatable :: Iend(:,:) ! Index matrix


   real*8       :: alfaTB  ! Diagonal tight binding param
   real*8       :: betaTB  ! Offdiagonal tight binding param
   real*8       :: gammaTB ! DFT-TB terms
   
      
!#ifdef TD_SIMPLE
   real*8,allocatable   :: chimerafock (:, :) !allocate in the central code
!#else
!   real*16, allocatable :: chimerafock (:, :)
!#endif

end module dftb_data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
