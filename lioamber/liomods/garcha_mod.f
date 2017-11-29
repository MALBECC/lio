!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module garcha_mod
!------------------------------------------------------------------------------!
       implicit real*8 (a-h,o-z)

      INCLUDE 'param.f'
      integer M,Md,natom,ntatom,NMAX,NCO,NUNP,igrid,igrid2
     >  ,Iexch,nsol,npas,npasw,idip,watermod,noconverge,
     > converge,ndiis,nang,propagator,NBCH
      integer ex_functional_id, ec_functional_id
      logical use_libxc
      integer restart_freq, energy_freq
      real*8 GOLD, TOLD, qmmmcut, dgtrig
!      parameter (nng=100)
      character*65 title
      character*20 basis,whatis,stdbas
      character*40 basis_set, fitting_set
      logical int_basis
      character*20 output,fcoord,fmulliken,frestart,frestartin,solv,
     > solv2
      character*4 ctype
      logical exists,MEMO,predcoef
      logical done(ntq),used,NORM,OPEN,DIRECT,VCINP,SHFT,DIIS
      logical done_fit(ntq)
      logical dens,EXTR,SVD
!write1
      logical Coul,GRAD,BSSE,integ,SVD1,sol,tipe
!Prop,GRAD,BSSE,integ,SVD1,sol,tipe
      logical exter,exter1,resp1,primera,writexyz,intsoldouble
      logical OPEN1
      logical dens1,integ1,sol1,free,free1

      logical writeforces

      logical cubegen_only,cube_dens,cube_orb,cube_elec, cube_sqrt_orb
      integer cube_res,cube_sel
      character*20 cube_dens_file,cube_orb_file,cube_elec_file


      dimension ighost(ntq), ighost1(ntq)
      real*8 e_(50,3),wang(50),e_2(116,3),wang2(116),e3(194,3), ! intg1 e intg2
     > wang3(194)                                               !
      integer Nr(0:54),Nr2(0:54)

      real*8, dimension (:,:), ALLOCATABLE :: r,v,rqm,d
      real*8, dimension (:), ALLOCATABLE ::  Em, Rm, pc
      integer, dimension (:), ALLOCATABLE :: Iz, nnat

      dimension isotop(54)!,Pm(nt)
      dimension Rm2(0:54), STR(880,0:21), FAC(0:16)
      dimension alpha(nss)
c Everything is dimensioned for 2 basis, normal and density
c ncf, lt,at,ct parameters for atomic basis sets
      dimension nshell(0:4)
      dimension Num(0:3),nlb(ng),nld(ngd),nshelld(0:4)
       integer iconst1,idip1,ipop1,ispin1,
     > icharge1,Nsol1,natsol1,Ll(3)

       real*8, dimension (:), ALLOCATABLE :: af
       real*8, dimension (:,:), ALLOCATABLE :: c,a,cx,ax,cd,ad,B
       integer, dimension (:), ALLOCATABLE :: Nuc,ncont,Nucx,ncontx,Nucd
     >  ,ncontd
        integer, dimension (:), ALLOCATABLE :: indexii, indexiid

c

       real*8, dimension (:), ALLOCATABLE :: RMM,RMM1,RMM2,RMM3
       real*8, dimension (:), ALLOCATABLE :: RMM_save ! TODO: delete after use (FFR)
       real*8, dimension (:), ALLOCATABLE :: rhoalpha,rhobeta
       real*8, dimension (:,:), ALLOCATABLE :: X, XX
       real*8, dimension (:), ALLOCATABLE :: old1,old2,old3

c
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0, pi5=34.9868366552497108D0,
     >    pi52=34.9868366552497108D0)
      parameter(pis32=5.56832799683170698E0,piss=3.14159265358979312E0,
     >          rpis=1.77245385090551588E0, pis5=34.9868366552497108E0,
     >    pis52=34.9868366552497108E0)


c Angular momenta : up to f functions ( could be easily extended if
c necessary)
c
c      common /fit/ Nang,dens1,integ1,Iexch1,igridc,igrid2c
c      common /cav/ a01,epsilon1,field1,exter1,Fx1,Fy1,Fz1
c
c      common /index/ ii,iid
c
      Data Num /1,3,6,10/
      dimension jatom(2,100),coef(100),dist(100,3),distt(100)
      integer ndis
      real*8 kjar,xini,xfinal

      integer, dimension(:), ALLOCATABLE :: natomc,nnps,nnpp,nnpd,nns
      integer, dimension(:), ALLOCATABLE :: nnd,nnp
      real*8, dimension (:), ALLOCATABLE :: atmin
      integer, dimension(:,:), ALLOCATABLE :: jatc
      integer kknums,kknumd
      integer, dimension (:), ALLOCATABLE :: kkind,kkinds
      real*8     rmax, rmaxs
      real*8, dimension (:), ALLOCATABLE :: cool
      real*4, dimension (:), ALLOCATABLE :: cools
c      parameter rmintsol=16.0D0
!
!
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

!-Variables for library reading.
      logical :: omit_bas
!-Variables for property calculations.
      logical :: fukui, dipole, lowdin, mulliken, print_coeffs

!     parameter (nng=100)
      integer :: nng, max_func
      integer, allocatable, dimension(:) :: ncf, lt
      real*8, allocatable, dimension(:) :: at, ct

!      dimension ncf(nng),lt(nng)
!      dimension at(nng),ct(nng)

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

!-Temporal variable for the transition of change RMM array int smaller arrays
	logical :: break_rmm
!-Arrays that going to reemplaze RMM
      real*8, dimension (:), ALLOCATABLE :: Fock_Hcore !contains core-Hamiltonian matrix, reemplaze RMM(M11)
      real*8, dimension (:), ALLOCATABLE :: Fock_Overlap ! reemplaze RMM(M5)
      real*8, dimension (:), ALLOCATABLE :: P_density ! reemplaze RMM(M1)
      real*8, dimension (:), ALLOCATABLE :: MO_coef_at, MO_coef_at_b
!Geometry optimizations
      logical :: steep !enables steepest decend algorithm
      real*8 :: Force_cut, Energy_cut, minimzation_steep !energy and force convergence crit and initial steep
      integer :: n_points ! number of points scaned for lineal search
      integer :: n_min_steeps !number of optimization steps
      integer :: charge
      logical :: lineal_search !enable lineal search


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      end module
