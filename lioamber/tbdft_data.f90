#include "complex_type.fh"
module tbdft_data
   implicit none

!Important input variables:                                                    !
! * tbdft_calc : integer, indicates the different calculation options (0 off - !
!               1 microcanonical dynamic - 2 DLVN rho0 generation -            !
!               3 DLVN dynamic )                                               !
! * MTB        : Integer, total of TB elements.                                !
! * n_biasTB   : Integer, number of electrodes.                                !
! * start_tdtb : Integer, initial TD step for the aplication of the bias.      !
! * end_tdtb   : Integer, final TD step for the aplication of the bias.        !
! * alfaTB     : Double precision real, fermi level for TB part                !
! * betaTB     : Double precision real, coupling parameter of TB elements,     !
!                this also control the band with of TB DOS.                    !
! * gammaTB    : Double precision real, coupling terms between TB and DFT part !
! * driving_rateTB : Double precision real, driving rate for DLVN calculation  !
!                                                                              !
!Input variables readed from the file gamma.in in the next order:              !
! * VbiasTB    : Double precision real array, store the bias of each electrode !
! * n_atperbias: Integer, number of DFT atoms per bias cupled.                 !
! * end_bTB    : Integer, number of coupling basis per atom.                   !
! * linkTB     : Integer array, list of atoms coupled with TB, separated by    !
!                electrode.                                                    !
! * basTB      : Integer array, list of the basis coupled in the LIO orther    !
! * gammaW     : Double precision, 2eight of each interaction between TB and   !
!                the coupled basis                                             !

   integer      :: tbdft_calc = 0
   integer      :: MTB    = 0
   integer      :: n_atTB = 0            ! Number of TB atoms per electrode.
   integer      :: MTBDFT = 0            ! Size of the TBDFT matrix
   integer      :: start_tdtb=0
   integer      :: end_tdtb=0
   integer      :: end_bTB
   integer      :: n_biasTB
   integer      :: n_atperbias
   real(kind=8) :: alfaTB
   real(kind=8) :: betaTB
   real(kind=8) :: gammaTB
   real(kind=8) :: driving_rateTB = 0.00d0
   integer        , allocatable :: Iend_TB(:,:)        ! Index matrix for coupling.
   integer        , allocatable :: linkTB(:,:)
   integer        , allocatable :: basTB(:)
   real(kind=8)   , allocatable :: rhoa_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8)   , allocatable :: rhob_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8)   , allocatable :: chimerafock (:,:,:) ! Allocated in the central code
   real(kind=8)   , allocatable :: gammaW(:)
   real(kind=8)   , allocatable :: VbiasTB(:)
   TDCOMPLEX      , allocatable :: rhold_AOTB(:,:,:)   ! rho in AO to calculate charges
   TDCOMPLEX      , allocatable :: rhonew_AOTB(:,:,:)  ! rho in AO to calculate charges
   TDCOMPLEX      , allocatable :: rhofirst_TB(:,:,:)  ! rhofirst for DLVN calc

end module tbdft_data
