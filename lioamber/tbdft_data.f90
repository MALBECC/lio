#include "complex_type.fh"
module tbdft_data
   implicit none

   integer      :: tbdft_calc = 0        ! 0 off - 1 microcanonical dynamic
                                         ! 2 DLVN rho0 generation - 3 DLVN
                                         ! dynamic
   integer      :: MTB    = 0            ! Total of TB elements.
   integer      :: n_atTB = 0            ! Number of TB atoms per electrode.
   integer      :: MTBDFT = 0            ! Size of the TBDFT matrix
   integer      :: start_tdtb=0          ! Initial time step for evolution of diagonal TB terms
   integer      :: end_tdtb=0            ! Final time step for evolution of diagonal TB terms
   integer      :: end_bTB               ! Number of basis coupled with TB part.
   integer      :: n_biasTB              ! Number of electrodes
   integer      :: n_atperbias           ! Number of coupled DFT atoms per bias
   real(kind=8) :: alfaTB                ! Fermi Energy
   real(kind=8) :: betaTB                ! Offdiagonal tight binding param
   real(kind=8) :: gammaTB               ! Coupling terms of TB-DFT
   real(kind=8) :: driving_rateTB = 0.00d0 ! Driving rate for DLVN
   integer        , allocatable :: Iend_TB(:,:)        ! Index matrix for coupling.
   integer        , allocatable :: linkTB(:,:)         ! Link atoms, separated by bias
   integer        , allocatable :: basTB(:)            ! Coupling basis in the LIO order
   real(kind=8)   , allocatable :: rhoa_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8)   , allocatable :: rhob_TBDFT(:,:)     ! Matrix to store rho TBDFT for TD
   real(kind=8)   , allocatable :: chimerafock (:,:,:) ! Allocated in the central code
   real(kind=8)   , allocatable :: gammaW(:)           ! gamma weight, per atom
   real(kind=8)   , allocatable :: VbiasTB(:)          ! Bias potential for each
                                                       ! electrode
   TDCOMPLEX      , allocatable :: rhold_AOTB(:,:,:)   ! rho in AO to calculate charges
   TDCOMPLEX      , allocatable :: rhonew_AOTB(:,:,:)  ! rho in AO to calculate charges
   TDCOMPLEX      , allocatable :: rhofirst_TB(:,:,:)  ! rhofirst for DLVN calc
end module tbdft_data
