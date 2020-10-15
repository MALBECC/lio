!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module garcha_mod

   implicit none
   integer :: natom, ntatom, NCO, NUNP, igrid, igrid2, Iexch, nsol, npas, &
              npasw, watermod, noconverge, converge, nang, propagator, NBCH
   integer :: ex_functional_id, ec_functional_id 
   logical :: use_libxc
   integer :: restart_freq, energy_freq, charge, timers
   character(len=40) :: fcoord, frestart, frestartin
   logical :: MEMO, OPEN, VCINP, writexyz, writeforces, dbug

   ! Grid options.
   LIODBLE :: e_(50,3), wang(50), e_2(116,3), wang2(116), e3(194,3), &
                   wang3(194), Rm2(0:54)
   integer      :: Nr(0:54), Nr2(0:54)

   ! System description.
   LIODBLE, dimension(:,:), allocatable ::  r, v, rqm, d
   LIODBLE, dimension(:)  , allocatable ::  Em, Rm, pc
   integer, dimension (:), allocatable :: Iz

   ! Density, Fock, and similar matrices and vectors.
   LIODBLE, dimension (:), allocatable :: Fmat_vec, Fmat_vec2, Pmat_vec, &
                                               Hmat_vec, Ginv_vec, Gmat_vec,  &
                                               Pmat_en_wgt, rhoalpha, rhobeta
   LIODBLE,allocatable,dimension(:,:) :: Smat, RealRho, sqsm,  &
                                              MO_coef_at, MO_coef_at_b

   ! For ehrenfest?
   logical :: doing_ehrenfest = .false.
   logical :: first_step
   LIODBLE :: total_time
   LIODBLE,allocatable,dimension(:)   :: atom_mass
   LIODBLE,allocatable,dimension(:,:) :: nucpos, nucvel
   LIODBLE,allocatable,dimension(:,:) :: qm_forces_ds, qm_forces_total

   ! Variables for property calculations.
   logical :: hybrid_forces_props = .false.
   logical :: print_coeffs = .false.
   logical :: dipole = .false.

   ! GPU OPTIONS and G2G
   logical      :: assign_all_functions, remove_zero_weights, energy_all_iterations
   LIODBLE :: free_global_memory, sphere_radius, little_cube_size
   integer      :: min_points_per_cube, max_function_exponent, gpu_level = 4

   ! Energy contributions
   LIODBLE :: Enucl
   LIODBLE,dimension(:)  ,allocatable :: Eorbs, Eorbs_b

   ! Variables for distance combination restrains
   integer :: number_restr, number_index
   integer, allocatable, dimension(:,:) :: restr_pairs
   integer, allocatable, dimension(:) ::  restr_index
   LIODBLE, allocatable, dimension(:) :: restr_k, restr_w, restr_r0

end module garcha_mod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
