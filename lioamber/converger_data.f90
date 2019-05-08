!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module converger_data

   implicit none

!  Covergence methods and criteria, as per input file.
!  Fock damping = 1, DIIS = 2, Hybrid convergence = 3, Biased DIIS = 4
!  Biased DIIS + Hybrid convergence = 5
   integer      :: conver_criter  = 2
   real(kind=8) :: gOld           = 10.0D0
   real(kind=8) :: damping_factor = 10.0D0

   ! DIIS and biased DIIS.
   integer      :: nDIIS          = 15
   logical      :: DIIS           = .true.
   real(kind=8) :: DIIS_bias      = 1.05D0

   ! Hybrid convergence
   logical      :: hybrid_converg = .false.
   real(kind=8) :: good_cut       = 1.0D-3
      
   ! Level shifting
   logical      :: level_shift    = .true.
   real(kind=8) :: lvl_shift_en   = 0.1D0

   ! Rho squared difference cut for each convergence strategy:
   real(kind=8) :: lvl_shift_cut  = 1.0D0
   real(kind=8) :: DIIS_start     = 1.0D0
   real(kind=8) :: bDIIS_start    = 1D-4
   real(kind=8) :: EDIIS_start    = 1D-20

   ! Tolerace for SCF convergence
   integer      :: nMax           = 100
   real(kind=8) :: tolD           = 1.0D-6
   real(kind=8) :: EtolD          = 1.0D-1

   ! Options for linear search. Rho_LS =1 activates
   ! linear search after failed convergence, =2 means
   ! only attempt linear search.
   integer      :: Rho_LS = 0

   ! Internal variables
   real(kind=8), allocatable :: fock_damped(:,:,:)
   real(kind=8)              :: rho_diff = 1.0D0

   ! Internal variables for DIIS (and variants)
   real(kind=8), allocatable :: fockm(:,:,:,:)
   real(kind=8), allocatable :: FP_PFm(:,:,:,:)
   real(kind=8), allocatable :: bcoef(:,:)
   real(kind=8), allocatable :: EMAT2(:,:,:)
   real(kind=8), allocatable :: energy_list(:)

   ! Internal variables for EDIIS
   integer                   :: nediis = 15
   real(kind=8), allocatable :: ediis_fock(:,:,:,:)
   real(kind=8), allocatable :: ediis_dens(:,:,:,:)
   real(kind=8), allocatable :: BMAT(:,:,:)
   real(kind=8), allocatable :: EDIIS_E(:)

   ! Internal variables for Linear Search
   logical                                   :: may_conv      = .true.
   real(kind=8)                              :: Elast         = 0.0D0
   real(kind=8)                              :: Pstepsize     = 0.0D0
   real(kind=8), dimension(:)  , allocatable :: rho_lambda1
   real(kind=8), dimension(:)  , allocatable :: rho_lambda0
   real(kind=8), dimension(:)  , allocatable :: rho_lambda1_alpha
   real(kind=8), dimension(:)  , allocatable :: rho_lambda0_alpha
   real(kind=8), dimension(:)  , allocatable :: rho_lambda1_betha
   real(kind=8), dimension(:)  , allocatable :: rho_lambda0_betha
   real(kind=8), dimension(:,:), allocatable :: P_hist

end module converger_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
