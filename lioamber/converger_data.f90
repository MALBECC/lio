!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module converger_data

   implicit none

!  Covergence methods and criteria, as per input file.
!  Fock damping = 1, DIIS = 2, Hybrid convergence = 3, Biased DIIS = 4

   integer      :: conver_criter  = 2
   integer      :: nDIIS          = 15
   logical      :: DIIS           = .true.
   logical      :: hybrid_converg = .false.
   real(kind=8) :: damping_factor = 10.0D0
   real(kind=8) :: DIIS_bias      = 1.05D0
   real(kind=8) :: gOld           = 10.0D0
   real(kind=8) :: good_cut       = 1.0D-3
   real(kind=8) :: tolD           = 1.0D-6
   real(kind=8) :: EtolD          = 1.0D-1


   ! Internal variables
   logical :: hagodiis
   real(kind=8), allocatable :: fock_damped(:,:,:)
   real(kind=8), allocatable :: fockm(:,:,:,:)
   real(kind=8), allocatable :: FP_PFm(:,:,:,:)
   real(kind=8), allocatable :: bcoef(:,:)
   real(kind=8), allocatable :: EMAT2(:,:,:)

end module converger_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
