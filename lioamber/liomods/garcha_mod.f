!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module garcha_mod
      implicit none
      INCLUDE 'param.f'
      integer natom,ntatom,NMAX,NCO,NUNP,igrid,igrid2
     >  ,Iexch,nsol,npas,npasw,watermod,noconverge,
     > converge,ndiis,nang,propagator,NBCH
      integer ex_functional_id, ec_functional_id
      logical use_libxc
      integer restart_freq, energy_freq
      real*8 GOLD, TOLD
      character*20 fcoord,fmulliken,frestart,frestartin,solv,solv2
      logical MEMO,predcoef
      logical OPEN,DIRECT,VCINP,DIIS
      logical sol
      logical primera,writexyz
      logical writeforces

      logical cubegen_only,cube_dens,cube_orb,cube_elec, cube_sqrt_orb
      integer cube_res,cube_sel
      character*20 cube_dens_file,cube_orb_file,cube_elec_file

      real*8 e_(50,3),wang(50),e_2(116,3),wang2(116),e3(194,3), ! intg1 e intg2
     > wang3(194)                                               !
      integer Nr(0:54),Nr2(0:54)

      real*8, dimension (:,:), ALLOCATABLE :: r,v,rqm,d
      real*8, dimension (:), ALLOCATABLE ::  Em, Rm, pc
      integer, dimension (:), ALLOCATABLE :: Iz

      real*8 :: Rm2(0:54)
c Everything is dimensioned for 2 basis, normal and density
c ncf, lt,at,ct parameters for atomic basis sets

      real*8, dimension (:), ALLOCATABLE :: Fmat_vec, Fmat_vec2,
     > Pmat_vec, Hmat_vec, Ginv_vec, Gmat_vec, Pmat_en_wgt
       real*8, dimension (:), ALLOCATABLE :: rhoalpha,rhobeta
       real*8, dimension (:,:), ALLOCATABLE :: X

       real*8 :: pi, pi32, rpi, pi5, pi52
       real*8 :: piss, pis32, rpis, pis5, pis52
       parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0, pi5=34.9868366552497108D0,
     >    pi52=34.9868366552497108D0)
       parameter(pis32=5.56832799683170698E0,piss=3.14159265358979312E0,
     >          rpis=1.77245385090551588E0, pis5=34.9868366552497108E0,
     >    pis52=34.9868366552497108E0)


c Angular momenta : up to f functions ( could be easily extended if
c necessary)

! FFR - My global variables
!------------------------------------------------------------------------------!
       real*8,allocatable,dimension(:,:) :: Smat
       real*8,allocatable,dimension(:,:) :: RealRho

       logical                           :: doing_ehrenfest=.false.
       logical                           :: first_step
       real*8,allocatable,dimension(:)   :: atom_mass
       real*8,allocatable,dimension(:,:) :: nucpos, nucvel
       real*8                            :: total_time
       real*8,allocatable,dimension(:,:) :: qm_forces_ds
       real*8,allocatable,dimension(:,:) :: qm_forces_total
!------------------------------------------------------------------------------!

!-Variables for hibrid damping-diis
      logical :: hybrid_converg
      double precision :: good_cut
      double precision :: Etold

!-Variables for property calculations.
      logical :: fukui, dipole, lowdin, mulliken, spinpop, print_coeffs

      integer :: nng, max_func

! GPU OPTIONS
      logical :: assign_all_functions, remove_zero_weights,
     >              energy_all_iterations
      real*8  :: free_global_memory, sphere_radius, little_cube_size
      integer :: min_points_per_cube, max_function_exponent

! Energy contributions
      real*8 :: Enucl
      real*8,dimension(:)  ,allocatable :: Eorbs, Eorbs_b
! need this for lowdin
      real*8,dimension(:,:),allocatable :: sqsm
!-Variables for distance combination restrain
      INTEGER :: number_restr, number_index
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: restr_pairs
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  restr_index
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: restr_k,restr_w,
     > restr_r0
!-Debug. Activates check of NaN in Fock and Rho
      Logical :: Dbug
      integer :: timers

      real*8, dimension (:,:), ALLOCATABLE :: MO_coef_at, MO_coef_at_b
!Geometry optimizations
      logical :: steep !enables steepest decend algorithm
      real*8 :: Force_cut, Energy_cut, minimzation_steep !energy and force convergence crit and initial steep
      integer :: n_points ! number of points scaned for lineal search
      integer :: n_min_steeps !number of optimization steps
      integer :: charge, gpu_level=4
      logical :: lineal_search !enable lineal search

! for properties calculation control
      logical :: calc_propM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      end module
