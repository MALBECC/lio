!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module garcha_mod

   implicit none
   include 'param.f90'

   integer :: natom, ntatom, NCO, NUNP, igrid, igrid2, Iexch, nsol, npas, &
              npasw, watermod, noconverge, converge, nang, propagator, NBCH
   integer :: ex_functional_id, ec_functional_id 
   logical :: use_libxc
   integer :: restart_freq, energy_freq
   character(len=20) :: fcoord, fmulliken, frestart, frestartin
   logical :: MEMO, OPEN, VCINP, writexyz, writeforces

   ! Cubegen options.
   logical :: cubegen_only, cube_dens, cube_orb, cube_elec, cube_sqrt_orb
   integer :: cube_res, cube_sel
   character(len=20) :: cube_dens_file, cube_orb_file, cube_elec_file

   ! Grid options.
   real(kind=8) :: e_(50,3), wang(50), e_2(116,3), wang2(116), e3(194,3), wang3(194)
   integer      :: Nr(0:54), Nr2(0:54)
   real(kind=8) :: Rm2(0:54)


   real(kind=8), dimension(:,:), allocatable ::  r, v, rqm, d
   real(kind=8), dimension(:)  , allocatable ::  Em, Rm, pc
   integer, dimension (:), allocatable :: Iz

   ! Density, Fock, and similar matrices and vectors.
   real(kind=8), dimension (:), allocatable :: Fmat_vec, Fmat_vec2, Pmat_vec, &
                                               Hmat_vec, Ginv_vec, Gmat_vec,  &
                                               Pmat_en_wgt, rhoalpha, rhobeta
   !real(kind=8), dimension (:,:), allocatable :: X

   real(kind=8),allocatable,dimension(:,:) :: Smat, RealRho

   ! For ehrenfest?
   logical :: doing_ehrenfest=.false.
   logical :: first_step
   real(kind=8) :: total_time
   real(kind=8),allocatable,dimension(:)   :: atom_mass
   real(kind=8),allocatable,dimension(:,:) :: nucpos, nucvel
   real(kind=8),allocatable,dimension(:,:) :: qm_forces_ds, qm_forces_total

   ! Variables for property calculations.
   logical :: fukui, dipole, lowdin, mulliken, print_coeffs
   logical :: becke = .false.

   ! GPU OPTIONS and G2G
   integer :: nng, max_func

   logical :: assign_all_functions, remove_zero_weights, energy_all_iterations
   real(kind=8)  :: free_global_memory, sphere_radius, little_cube_size
   integer :: min_points_per_cube, max_function_exponent

   ! Energy contributions
   real(kind=8) :: Enucl
   real(kind=8),dimension(:)  ,allocatable :: Eorbs, Eorbs_b
   real(kind=8),dimension(:,:),allocatable :: sqsm

   ! Variables for distance combination restrains
   integer :: number_restr, number_index
   integer, allocatable, dimension(:,:) :: restr_pairs
   integer, allocatable, dimension(:) ::  restr_index
   real(kind=8), allocatable, dimension(:) :: restr_k, restr_w, restr_r0

   ! Debug. Activates check of NaN in Fock and Rho
   logical :: Dbug
   integer :: timers

   real(kind=8), dimension (:,:), allocatable :: MO_coef_at, MO_coef_at_b

   ! Geometry optimizations
   logical :: steep !enables steepest decend algorithm
   real(kind=8) :: Force_cut, Energy_cut, minimzation_steep !energy and force convergence crit and initial steep
   integer :: n_points ! number of points scaned for lineal search
   integer :: n_min_steeps !number of optimization steps
   integer :: charge, gpu_level=4
   logical :: lineal_search !enable lineal search

   ! For properties calculation control
   logical :: calc_propM

   ! For pbe0 functional
   logical :: PBE0 = .false.

end module garcha_mod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
