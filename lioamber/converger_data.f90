#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module converger_data

   implicit none

!  Covergence methods and criteria, as per input file.
!  Fock damping = 1, DIIS = 2, Hybrid convergence = 3, Biased DIIS = 4
!  Biased DIIS + Hybrid convergence = 5
   integer      :: conver_method  = 2
   LIODBLE :: gOld           = 10.0D0
   LIODBLE :: damping_factor = 10.0D0

   ! DIIS and biased DIIS.
   integer      :: nDIIS          = 15
   logical      :: DIIS           = .true.
   LIODBLE :: DIIS_bias      = 1.05D0

   ! Hybrid convergence
   logical      :: hybrid_converg = .false.
   LIODBLE :: good_cut       = 1.0D-3

   ! Level shifting
   logical      :: level_shift    = .false.
   LIODBLE :: lvl_shift_en   = 0.25D0
   LIODBLE :: lvl_shift_cut  = 0.005D0

   ! DIIS error cut for each convergence strategy:
   LIODBLE :: EDIIS_start    = 1D-20
   LIODBLE :: DIIS_start     = 0.01D0
   LIODBLE :: bDIIS_start    = 1D-3

   ! Tolerace for SCF convergence
   integer      :: nMax           = 100
   LIODBLE :: tolD           = 1.0D-6
   LIODBLE :: EtolD          = 1.0D-1

   ! Options for linear search. Rho_LS =1 activates
   ! linear search after failed convergence, =2 means
   ! only attempt linear search.
   integer      :: Rho_LS = 0

   ! Internal variables
   LIODBLE, allocatable :: fock_damped(:,:,:)
   LIODBLE              :: rho_diff      = 1.0D0
   LIODBLE              :: DIIS_error    = 100.0D0
   logical                   :: DIIS_started  = .false.
   logical                   :: EDIIS_started = .false.
   logical                   :: bDIIS_started = .false.

   ! Internal variables for DIIS (and variants)
   LIODBLE, allocatable :: fockm(:,:,:,:)
   LIODBLE, allocatable :: FP_PFm(:,:,:,:)
   LIODBLE, allocatable :: bcoef(:)
   LIODBLE, allocatable :: EMAT(:,:)
   LIODBLE, allocatable :: energy_list(:)

   ! Internal variables for EDIIS
   integer                   :: nediis          = 15
   logical                   :: EDIIS_not_ADIIS = .true.
   LIODBLE, allocatable :: ediis_fock(:,:,:,:)
   LIODBLE, allocatable :: ediis_dens(:,:,:,:)
   LIODBLE, allocatable :: BMAT(:,:)
   LIODBLE, allocatable :: EDIIS_E(:)
   LIODBLE, allocatable :: EDIIS_coef(:)

   ! Internal variables for Linear Search
   logical                   :: first_call = .true.
   LIODBLE              :: Elast      = 1000.0D0
   LIODBLE              :: Pstepsize  = 1.0D0
   LIODBLE, allocatable :: rho_lambda1(:)
   LIODBLE, allocatable :: rho_lambda0(:)
   LIODBLE, allocatable :: rhoa_lambda1(:)
   LIODBLE, allocatable :: rhoa_lambda0(:)
   LIODBLE, allocatable :: rhob_lambda1(:)
   LIODBLE, allocatable :: rhob_lambda0(:)

end module converger_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
