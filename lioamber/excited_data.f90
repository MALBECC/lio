#include "datatypes/datatypes.fh"
module excited_data
   implicit none

   ! Save Integrals in Memory
   integer :: libint_recalc = 0 ! Recalculated Integrals, 1 = Save in Memory

   ! Save GS density and derivatives in each point of the group 
   integer :: density_point_save = 0 ! Save in memory, 1 = Recalculating in each point

   ! Linear Response
   logical :: lresp = .false.
   integer :: nstates = 4 
   integer :: max_subs = 400 ! max subspace in LR
   logical :: fittExcited = .false. ! poner de input
   LIODBLE :: tolv = 1.0d-7 ! conv. crit. to vectors ! pon inp
   LIODBLE :: tole = 1.0d-7 ! conv. crit. to energy !pon inp
   integer :: root = 0 ! Relaxed Density Matrix of Excited State root
   logical :: print_trdip = .false. ! Print Transition Dipole Moment of all states

   ! Cubegen files in excited states
   logical :: cgPe = .false. !Excited Density
   logical :: cgPd = .false. !Diference Density
   logical :: cgPt = .false. !Transition Density
   logical :: cgMO1= .false. !The principal contribution of MOs
   logical :: cgPg = .false. !Ground Density

   ! Truncated MOs
   integer :: trunc_mos = 0 ! 0 = NO, 1 = FCA, 2 = Reduced MOs
   integer :: nfo = 3 ! occupied in FCA
   integer :: nfv = 3 ! virtual in FCA
   LIODBLE :: thres_occ = 0.6d0 ! threshold occupied in Reduced MOs
   LIODBLE :: thres_vir = 0.4d0 ! threshold virtual in Reduced MOs
   integer, dimension(:), allocatable :: map_occ, map_vir ! map (small) -> big indexes

   ! Energy-Specific TDA
   logical :: estda      = .false.
   LIODBLE :: energy_min = 10.0d0! 3.6 Ha. = 100 eV ( X ray )
   LIODBLE :: d_energy   = 1.0d0 ! 
   integer :: window     = -1    ! Energy window id 

   ! Basis Change
   LIODBLE, dimension(:,:), allocatable :: Coef_trans, Cocc 
   LIODBLE, dimension(:,:), allocatable :: Cocc_trans, Cvir
   LIODBLE, dimension(:,:), allocatable :: Cvir_trans

   ! Using Last step as Initial Guess in LR
   logical :: use_last = .false.
   LIODBLE, dimension(:,:), allocatable :: guessLR

   ! Excited States Forces
   logical :: excited_forces = .false.
   LIODBLE, dimension(:,:), allocatable :: for_exc

   ! Save Relaxed Density matrix in vector form in order to obtain
   ! differents properties
   LIODBLE, dimension(:), allocatable   :: pack_dens_exc

   ! Trajectory Surface Hopping: Only coupling between GS and 1st ES
   logical :: TSH = .false.
   integer :: tsh_Jstate, tsh_Kstate
   LIODBLE :: dE_accum, lambda, tsh_time_dt
   LIODBLE, dimension(:,:), allocatable :: gamma_old
   complex(kind=8) :: B_old = cmplx(0.0d0,0.0d0,8)
   complex(kind=8), dimension(:), allocatable :: tsh_coef

   ! Second Linear Response
   logical :: second_LR = .false.
   integer :: state_LR  = 1
   LIODBLE :: lambda_LR = 0.05d0
   LIODBLE :: Ctol = 0.01d0
   LIODBLE, dimension(:,:), allocatable :: Tdip_save
   
   ! ESA with first LR
   logical :: ESAfosc = .false.

   ! Save Transitions Vectors
   logical :: save_tlr = .false.
   logical :: read_other = .false.

end module excited_data
