!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_data
   implicit none


   logical      :: dftb_calc =.false.   ! Logical indicator for dftb calculation
   logical      :: TBsave=.false.       ! Logical to save TB rho
   logical      :: TBload=.false.       ! Logical to use TB rho save
   integer      :: MTB   = 0            ! Size of the two tight-binding subatrices
   integer      :: MDFTB = 0            ! Size of the DFT-TB matrix
   integer      :: start_tdtb=0         ! Initial time step for evolution of
                                        !                  diagonal TB terms
   integer      :: end_tdtb=0           ! Final time step for evolution of
                                        ! diagonal TB terms
   integer      :: end_bTB            ! Index matrix size
   integer, allocatable :: Iend_TB(:,:) ! Index matrix


   real*8               :: alfaTB             ! Fermi Energy
   real*8               :: betaTB             ! Offdiagonal tight binding param
   real*8               :: gammaTB            ! DFT-TB terms
   real*8               :: Vbias_TB           ! Bias potential
   real*8               :: chargeA_TB
   real*8               :: chargeB_TB
   real*8               :: chargeM_TB
   real*8, allocatable  :: rho_aDFTB(:,:)     ! Matrix to store rho DFTB for TD
   real*8, allocatable  :: rho_bDFTB(:,:)     ! Matrix to store rho DFTB for TD

   real*8,allocatable   :: chimerafock (:, :) ! Allocated in the central code
#ifdef TD_SIMPLE
   complex*8, allocatable  :: rhold_AOTB(:,:)     ! rho in AO to calculate charges
   complex*8, allocatable  :: rhonew_AOTB(:,:)     ! rho in AO to calculate charges
#else
   complex*16, allocatable  :: rhold_AOTB(:,:)     ! rho in AO to calculate charges
   complex*16, allocatable  :: rhonew_AOTB(:,:)     ! rho in AO to calculate charges
#endif
end module dftb_data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
