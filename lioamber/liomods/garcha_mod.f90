!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module garcha_mod

   implicit none
   integer :: natom, ntatom, NCO, NUNP, igrid, igrid2, Iexch, nsol, npas, &
              npasw, watermod, noconverge, converge, nang, propagator, NBCH
   integer :: ex_functional_id, ec_functional_id 
   logical :: use_libxc
   integer :: restart_freq, energy_freq, charge, timers
   character(len=40) :: fcoord, fmulliken, frestart, frestartin
   logical :: MEMO, OPEN, VCINP, writexyz, writeforces, dbug

   ! Grid options.
   real(kind=8) :: e_(50,3), wang(50), e_2(116,3), wang2(116), e3(194,3), &
                   wang3(194), Rm2(0:54)
   integer      :: Nr(0:54), Nr2(0:54)

   ! System description.
   real(kind=8), dimension(:,:), allocatable ::  r, v, rqm, d
   real(kind=8), dimension(:)  , allocatable ::  Em, Rm, pc
   integer, dimension (:), allocatable :: Iz

   ! Density, Fock, and similar matrices and vectors.
   real(kind=8), dimension (:), allocatable :: Fmat_vec, Fmat_vec2, Pmat_vec, &
                                               Hmat_vec, Ginv_vec, Gmat_vec,  &
                                               Pmat_en_wgt, rhoalpha, rhobeta
   real(kind=8),allocatable,dimension(:,:) :: Smat, RealRho, sqsm,  &
                                              MO_coef_at, MO_coef_at_b

   ! For ehrenfest?
   logical :: doing_ehrenfest = .false.
   logical :: first_step
   real(kind=8) :: total_time
   real(kind=8),allocatable,dimension(:)   :: atom_mass
   real(kind=8),allocatable,dimension(:,:) :: nucpos, nucvel
   real(kind=8),allocatable,dimension(:,:) :: qm_forces_ds, qm_forces_total

   ! Variables for property calculations.
   logical :: fukui, dipole, lowdin, mulliken, print_coeffs, calc_propM
   logical :: becke = .false.

   ! GPU OPTIONS and G2G
   logical      :: assign_all_functions, remove_zero_weights, energy_all_iterations
   real(kind=8) :: free_global_memory, sphere_radius, little_cube_size
   integer      :: min_points_per_cube, max_function_exponent, gpu_level = 4

   ! Energy contributions
   real(kind=8) :: Enucl
   real(kind=8),dimension(:)  ,allocatable :: Eorbs, Eorbs_b

   ! Variables for distance combination restrains
   integer :: number_restr, number_index
   integer, allocatable, dimension(:,:) :: restr_pairs
   integer, allocatable, dimension(:) ::  restr_index
   real(kind=8), allocatable, dimension(:) :: restr_k, restr_w, restr_r0

   ! For pbe0 functional
   logical :: PBE0 = .false.

end module garcha_mod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
