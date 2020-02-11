!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! DIRECT VERSION
! Calls all integrals generator subroutines : 1 el integrals,
! 2 el integrals, exchange fitting , so it gets S matrix, F matrix
! and P matrix in lower storage mode (symmetric matrices)
!
! Dario Estrin, 1992
!------------------------------------------------------------------------------!
! Modified to f90
! Nick, 2017
!------------------------------------------------------------------------------!
! Header with new format. Added comments on how to proceed with a cleanup of
! of the subroutines. Other things to do:
! TODO: change to 3 space indentation.
! TODO: break at line 80.
! TODO: change to lowercase (ex: implicit none)
!
! This log can be removed once all improvements have been made.
! FFR, 01/2018
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
   use ehrensubs , only: ehrendyn_init
   use garcha_mod, only : NCO, natom, number_restr, MEMO, &
                          igrid, energy_freq, converge, noconverge, lowdin,    &
                          cubegen_only, VCINP, primera, Nunp, igrid2,    &
                          nsol, r, pc, Iz, &
                          Eorbs, Dbug, doing_ehrenfest, &
                          MO_coef_at, MO_coef_at_b, Smat, &
                          rhoalpha, rhobeta, OPEN, RealRho, d, ntatom,  &
                          Eorbs_b, npas, X, npasw, Fmat_vec, Fmat_vec2,        &
                          Ginv_vec, Gmat_vec, Hmat_vec, Pmat_en_wgt, Pmat_vec, &
                          sqsm, PBE0
   use ECP_mod, only : ecpmode
   use field_data, only: field, fx, fy, fz
   use field_subs, only: field_calc, field_setup_old
   use faint_cpu, only: int1, intsol, int2, int3mem, int3lu
   use tbdft_data, only : tbdft_calc, MTBDFT, MTB,rhoa_tbdft,rhob_tbdft,n_biasTB
   use tbdft_subs, only : getXY_TBDFT, build_chimera_TBDFT, extract_rhoDFT, &
                          construct_rhoTBDFT, tbdft_scf_output,write_rhofirstTB
   use transport_data, only: generate_rho0
   use cubegen       , only: cubegen_matin, cubegen_write
   use mask_ecp      , only: ECP_fock, ECP_energy
   use typedef_sop   , only: sop              ! Testing SOP
   use fockbias_subs , only: fockbias_loads, fockbias_setmat, fockbias_apply
   use SCF_aux       , only: seek_nan, standard_coefs, messup_densmat, fix_densmat
   use liosubs_math  , only: transform
   use converger_data, only: Rho_LS, nMax
   use converger_subs, only: converger_init, converger_fock, converger_setup, &
                             converger_check, rho_ls_init, do_rho_ls,         &
                             rho_ls_switch
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_r
   use trans_Data    , only: gaussian_convert, rho_exc, translation
   use initial_guess_subs, only: get_initial_guess
   use fileio       , only: write_energies, write_energy_convergence, &
                            write_final_convergence, write_ls_convergence, &
                            movieprint
   use fileio_data  , only: verbose
   use basis_data   , only: kkinds, kkind, cools, cool, Nuc, nshell, M, MM, c_raw
   use basis_subs, only: neighbour_list_2e
   use excited_data, only: libint_recalc
   use excitedsubs, only: ExcProp
   use dftd3, only: dftd3_energy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   implicit none
   ! E is the total SCF energy.
   ! The others are fock and rho operators alpha and beta (FOCK/RHO Alpha/Beta
   ! OPerator). In the case of closed shell, rho_aop and fock_aop contain full
   ! Fock and Rho matrices.
   real(kind=8)  , intent(inout)           :: E
   type(operator), intent(inout)           :: rho_aop, fock_aop
   type(operator), intent(inout), optional :: rho_bop, fock_bop


   integer :: niter
   logical :: converged     = .false.
   logical :: changed_to_LS = .false.
   real*8  :: sq2
   integer :: igpu

!  The following two variables are in a part of the code that is never
!  used. Check if these must be taken out...
   real*8, allocatable :: morb_coefon(:,:)

!------------------------------------------------------------------------------!
!  TBDFT: variables to use as input for some subroutines instead of M and NCO
   integer :: M_f
   integer :: NCOa_f
   integer :: NCOb_f

   real*8, allocatable :: rho_a0(:,:), rho_b0(:,:)
   real*8, allocatable :: fock_a0(:,:), fock_b0(:,:)
   real*8, allocatable :: rho_a(:,:), rho_b(:,:)
   real*8, allocatable :: fock_a(:,:), fock_b(:,:)
   real*8, allocatable :: morb_coefat(:,:)
   real*8, allocatable :: X_min(:,:)
   real*8, allocatable :: Y_min(:,:)
   real*8, allocatable :: X_min_trans(:,:)
   real*8, allocatable :: Y_min_trans(:,:)
   real*8, allocatable :: morb_energy(:)
   integer             :: ii, jj, kk, kkk

!------------------------------------------------------------------------------!
! FFR variables
   type(sop)           :: overop
   real*8, allocatable :: tmpmat(:,:)
   real*8  :: HL_gap = 10.0D0

!------------------------------------------------------------------------------!
! Energy contributions and convergence

   real*8 :: E1          ! kinetic + nuclear attraction + e-/MM charge
                         !    interaction + effective core potetial
   real*8 :: E1s = 0.0D0 ! kinetic + nuclear attraction + effective core
                         !    potetial
   real*8 :: E2          ! Coulomb (e- - e-)
   real*8 :: Eecp        ! Efective core potential
   real*8 :: En          ! nuclear-nuclear repulsion
   real*8 :: Ens         ! MM point charge-nuclear interaction
   real*8 :: Es          ! ???
   real*8 :: E_restrain  ! distance restrain
   real*8 :: Exc         ! exchange-correlation
   real*8 :: Etrash      ! auxiliar variable
   real*8 :: Evieja      !


   ! Base change matrices (for ON-AO changes).
   type(cumat_r)       :: Xmat, Ymat

   ! TODO : Variables to eliminate...
   real*8, allocatable :: xnano(:,:)
   integer :: M2

   ! Carlos: Open shell, variables.
   real*8              :: ocupF
   integer             :: NCOa, NCOb

   ! Variables related to VdW
   real(kind=8) :: E_dftd

   ! Variables-PBE0
   real(kind=8) :: Eexact
   real(kind=8), allocatable :: FockEE_a0(:,:), FockEE_b0(:,:)

   call g2g_timer_start('SCF_full')
   call g2g_timer_start('SCF')
   call g2g_timer_sum_start('SCF')
   call g2g_timer_sum_start('Initialize SCF')

   changed_to_LS=.false. ! LINSEARCH
   call rho_ls_init(open, MM)

   E=0.0D0
   E1=0.0D0
   En=0.0D0
   E2=0.0D0
   Es=0.0D0
   Eecp=0.d0
   Ens=0.0D0
   E_restrain=0.d0
   E_dftd=0.0D0
   Eexact=0.D0

   ! Distance Restrain
   IF (number_restr.GT.0) THEN
      call get_restrain_energy(E_restrain)
      WRITE(*,*) "DISTANCE RESTRAIN ADDED TO FORCES"
   END IF

   !carlos: ocupation factor
   !carlos: NCOa works in open shell and close shell
   NCOa   = NCO
   NCOa_f = NCOa
   ocupF  = 2.0d0
   if (OPEN) then
      ! Number of OM down
      NCOb   = NCO + Nunp
      NCOb_f = NCOb
      ocupF = 1.0d0
      allocate(rho_b0(M,M),fock_b0(M,M))
   endif
   allocate(fock_a0(M,M), rho_a0(M,M))

   M_f = M
   if (tbdft_calc /= 0) then
      M_f    = MTBDFT
      NCOa_f = NCOa + MTB / 2
      if (OPEN) NCOb_f = NCOb + MTB / 2
   endif

   allocate(fock_a(M_f,M_f), rho_a(M_f,M_f))
   allocate(morb_energy(M_f), morb_coefat(M_f,M_f))
   if (OPEN) then
      allocate(fock_b(M_f,M_f), rho_b(M_f,M_f))
   end if

!------------------------------------------------------------------------------!
! TODO: damp and gold should no longer be here??
! TODO: Qc should probably be a separated subroutine? Apparently it is only
!       used in dipole calculation so...it only adds noise to have it here.
! TODO: convergence criteria should be set at namelist/keywords setting

      Evieja=0.d0
      niter=0

!------------------------------------------------------------------------------!
! TODO: this whole part which calculates the non-electron depending terms of
!       fock and the overlap matrix should probably be in a separated sub.
!       (diagonalization of overlap, starting guess, should be taken out)
!
! Reformat from here...

! Nano: calculating neighbour list helps to make 2 electrons integral scale
! linearly with natoms/basis
!
      call neighbour_list_2e(natom, ntatom, r, d)

! -Create integration grid for XC here
! -Assign points to groups (spheres/cubes)
! -Assign significant functions to groups
! -Calculate point weights
!
      call g2g_timer_sum_start('Exchange-correlation grid setup')
      call g2g_reload_atom_positions(igrid2, Iz)
      call g2g_timer_sum_stop('Exchange-correlation grid setup')

      call aint_query_gpu_level(igpu)
      if (igpu.gt.1) call aint_new_step()

! Calculate 1e part of F here (kinetic/nuc in int1, MM point charges
! in intsol)
!
      call g2g_timer_sum_start('1-e Fock')
      call g2g_timer_sum_start('Nuclear attraction')
      call int1(En, Fmat_vec, Hmat_vec, Smat, d, r, Iz, natom, &
                ntatom)
      call ECP_fock( MM, Hmat_vec )

! Other terms
!
      call g2g_timer_sum_stop('Nuclear attraction')
      if(nsol.gt.0.or.igpu.ge.4) then
          call g2g_timer_sum_start('QM/MM')
       if (igpu.le.1) then
          call g2g_timer_start('intsol')
          call intsol(Pmat_vec, Hmat_vec, Iz, pc, r, d, natom, ntatom, &
                      E1s, Ens, .true.)
          call g2g_timer_stop('intsol')
        else
          call aint_qmmm_init(nsol,r,pc)
          call g2g_timer_start('aint_qmmm_fock')
          call aint_qmmm_fock(E1s,Ens)
          call g2g_timer_stop('aint_qmmm_fock')
        endif
          call g2g_timer_sum_stop('QM/MM')
      endif


! Initialization of libint
      if ( PBE0 ) then
         call g2g_timer_sum_start('Libint init')
         call g2g_libint_init(c_raw,libint_recalc)
         call g2g_timer_sum_stop('Libint init')
      endif

! test
! TODO: test? remove or sistematize
!
      E1=0.D0
      do kk=1,MM
        E1 = E1 + Pmat_vec(kk) * Hmat_vec(kk)
      enddo
      call g2g_timer_sum_stop('1-e Fock')


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! OVERLAP DIAGONALIZATION
! TODO: Simplify, this has too much stuff going on...
! (maybe trans mats are not even necessary?)
!
        if (allocated(X_min)) deallocate(X_min)
        if (allocated(Y_min)) deallocate(Y_min)
        if (allocated(X_min_trans)) deallocate(X_min_trans)
        if (allocated(Y_min_trans)) deallocate(Y_min_trans)

        allocate(X_min(M,M), Y_min(M,M), X_min_trans(M,M), Y_min_trans(M,M))

        call overop%Sets_smat( Smat )
        if (lowdin) then
!          TODO: inputs insuficient; there is also the symetric orthog using
!                3 instead of 2 or 1. Use integer for onbasis_id
           call overop%Gets_orthog_4m( 2, 0.0d0, X_min, Y_min, X_min_trans, Y_min_trans)
        else
           call overop%Gets_orthog_4m( 1, 0.0d0, X_min, Y_min, X_min_trans, Y_min_trans)
        end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!  Fockbias setup
        if ( allocated(sqsm) ) deallocate(sqsm)
        if ( allocated(tmpmat) ) deallocate(tmpmat)
        allocate( sqsm(M,M), tmpmat(M,M) )
        call overop%Gets_orthog_2m( 2, 0.0d0, tmpmat, sqsm )
        call fockbias_loads( natom, nuc )
        call fockbias_setmat( sqsm )
        deallocate( tmpmat )


!TBDFT: Dimensions of Xmat and Ymat are modified for TBDFT.
!
! TODO: this is nasty, a temporary solution would be to have a Msize variable
!       be assigned M (or, even better, "basis_size") or MTBDFT
!       ("basis_size_dftb") according to the case
! Uses arrays fock_a y rho_a as temporary storage to initialize Xmat and Ymat.

   call getXY_TBDFT(M, X_min, Y_min, fock_a, rho_a)
   call Xmat%init(M_f, fock_a)
   call Ymat%init(M_f, rho_a)

   deallocate(X_min, Y_min, X_min_trans, Y_min_trans)

! Generates starting guess
!
   if ( (.not.VCINP) .and. primera ) then
      call get_initial_guess(M, MM, NCO, NCOb, &
                             Xmat%matrix(MTB+1:MTB+M,MTB+1:MTB+M),        &
                             Hmat_vec, Pmat_vec, rhoalpha, rhobeta, OPEN, &
                             natom, Iz, nshell, Nuc)
      primera = .false.
   endif

!----------------------------------------------------------!
! Precalculate two-index (density basis) "G" matrix used in density fitting
! here (S_ij in Dunlap, et al JCP 71(8) 1979).
! Also, pre-calculate G^-1 if G is not ill-conditioned.
      call g2g_timer_sum_start('Coulomb G matrix')
      call int2(Gmat_vec, Ginv_vec, r, d, ntatom)
      call g2g_timer_sum_stop('Coulomb G matrix')

! Precalculate three-index (two in MO basis, one in density basis) matrix
! used in density fitting / Coulomb F element calculation here
! (t_i in Dunlap)
!
      call aint_query_gpu_level(igpu)
      if (igpu.gt.2) then
        call aint_coulomb_init()
      endif
      if (igpu.eq.5) MEMO = .false.
      !MEMO=.true.
      if (MEMO) then
         call g2g_timer_start('int3mem')
         call g2g_timer_sum_start('Coulomb precalc')
!        Large elements of t_i put into double-precision cool here
!        Size criteria based on size of pre-factor in Gaussian Product Theorem
!        (applied to MO basis indices)
         call int3mem(r, d, natom, ntatom)
         call g2g_timer_stop('int3mem')
         call g2g_timer_sum_stop('Coulomb precalc')
      endif
!
!##########################################################!
! TODO: ...to here
!##########################################################!
!
!
!
!------------------------------------------------------------------------------!
! TODO: the following comment is outdated? Also, hybrid_converg switch should
!       be handled differently.
!
! Now, damping is performed on the density matrix
! The first 4 iterations ( it may be changed, if necessary)
! when the density is evaluated on the grid, the density
! matrix is used ( slower), after that it is calculated
! using the vectors . Since the vectors are not damped,
! only at the end of the SCF, the density matrix and the
! vectors are 'coherent'

      call g2g_timer_sum_stop('Initialize SCF')

!------------------------------------------------------------------------------!
! TODO: Maybe evaluate conditions for loop continuance at the end of loop
!       and condense in a single "keep_iterating" or something like that.
   if (verbose > 1) then
      write(*,*)
      write(*,'(A)') "Starting SCF cycles."
   endif

   converged = .false.
   call converger_init( M_f, OPEN )

   do 999 while ( (.not. converged) .and. (niter <= nMax) )
      call g2g_timer_start('Total iter')
      call g2g_timer_sum_start('Iteration')
      call g2g_timer_sum_start('Fock integrals')

      niter  = niter +1

      ! Test for NaN
      if (Dbug) call SEEK_NaN(Pmat_vec,1,MM,"RHO Start")
      if (Dbug) call SEEK_NaN(Fmat_vec,1,MM,"FOCK Start")

      ! Computes Coulomb part of Fock, and energy on E2
      call g2g_timer_sum_start('Coulomb fit + Fock')
      E2 = 0.0D0
      call int3lu(E2, Pmat_vec, Fmat_vec2, Fmat_vec, Gmat_vec, Ginv_vec, &
                  Hmat_vec, open, MEMO)
      call g2g_timer_sum_pause('Coulomb fit + Fock')

      if (Dbug) then
         call SEEK_NaN(Fmat_vec,1,MM,"FOCK Coulomb")
         if (open) call SEEK_NaN(Fmat_vec2,1,MM,"FOCK B Coulomb")
      endif

      ! XC integration / Fock elements
      call g2g_timer_sum_start('Exchange-correlation Fock')
      Exc = 0.0D0
      call g2g_solve_groups(0,Exc,0)
      call g2g_timer_sum_pause('Exchange-correlation Fock')

      ! Test for NaN
      if (Dbug) then
         call SEEK_NaN(Fmat_vec,1,MM,"FOCK Ex-Corr")
         if (open) call SEEK_NaN(Fmat_vec2,1,MM,"FOCK B Ex-Corr")
      endif


      ! Calculates 1e energy contributions (including solvent)
      E1 = 0.0D0
      if (generate_rho0) then
         ! REACTION FIELD CASE
         if (field) call field_setup_old(1.0D0, 0, fx, fy, fz)
         call field_calc(E1, 0.0D0, Pmat_vec, Fmat_vec2, Fmat_vec, r, d, Iz,&
                         natom, ntatom, open)
         do kk = 1, MM
            E1 = E1 + Pmat_vec(kk) * Hmat_vec(kk)
         enddo
      else
         do kk=1,MM
            E1 = E1 + Pmat_vec(kk) * Hmat_vec(kk)
         enddo
      endif

      ! Calculates total energy
      E = E1 + E2 + En + Exc
      call g2g_timer_sum_pause('Fock integrals')


      if (OPEN) then
         call spunpack_rho('L', M, rhoalpha , rho_a0)
         call spunpack_rho('L', M, rhobeta  , rho_b0)
         call spunpack(    'L', M, Fmat_vec , fock_a0)
         call spunpack(    'L', M, Fmat_vec2, fock_b0)
         call fockbias_apply( 0.0d0, fock_a0)
         call fockbias_apply( 0.0d0, fock_b0)
      else
         call spunpack_rho('L', M, Pmat_vec, rho_a0)
         call spunpack(    'L', M, Fmat_vec, fock_a0)
         call fockbias_apply(0.0d0, fock_a0)
      end if

!     EXACT EXCHANGE - PBE0
      if ( PBE0 ) then
         call g2g_timer_sum_start('Exact Exchange Fock')

         if (allocated(FockEE_a0)) deallocate(FockEE_a0)
         allocate(FockEE_a0(M,M)); FockEE_a0 = 0.0d0

         if ( OPEN ) then
           if (allocated(FockEE_b0)) deallocate(FockEE_b0)
           allocate(FockEE_b0(M,M)); FockEE_b0 = 0.0d0
           call g2g_exact_exchange_open(rho_a0,rho_b0,FockEE_a0,FockEE_b0)
           fock_a0 = fock_a0 - 0.25D0 * FockEE_a0
           fock_b0 = fock_b0 - 0.25D0 * FockEE_b0
         else
           call g2g_exact_exchange(rho_a0,FockEE_a0)
           fock_a0 = fock_a0 - 0.25D0 * FockEE_a0
         endif

         call g2g_timer_sum_pause('Exact Exchange Fock')
      endif

      if (tbdft_calc == 0) then
         fock_a = fock_a0
         rho_a  = rho_a0
         if (OPEN) fock_b = fock_b0
         if (OPEN) rho_b  = rho_b0
      else
         ! TBDFT: We extract rho and fock before convergence acceleration
         ! routines. Then, Fock and Rho for TBDFT are builded.
         call build_chimera_TBDFT (M, fock_a0, fock_a, natom)
         call construct_rhoTBDFT(M, rho_a, rho_a0 ,rhoa_tbdft, niter,OPEN)
         if (OPEN) then
            call build_chimera_TBDFT(M, fock_b0, fock_b, natom)
            call construct_rhoTBDFT(M, rho_b, rho_b0 ,rhob_tbdft,niter, OPEN)
         end if
      endif

      ! Stores matrices in operators, and sets up matrices in convergence
      ! acceleration algorithms (DIIS/EDIIS).
      call rho_aop%Sets_data_AO(rho_a)
      call fock_aop%Sets_data_AO(fock_a)

      if (OPEN) then
         call rho_bop%Sets_data_AO(rho_b)
         call fock_bop%Sets_data_AO(fock_b)
         call converger_setup(niter, M_f, rho_aop, fock_aop, E,  Xmat, Ymat, &
                              rho_bop, fock_bop)
      else
         call converger_setup(niter, M_f, rho_aop, fock_aop, E,  Xmat, Ymat)
      endif

      ! Convergence accelerator processing.
      ! In closed shell, rho_a is the total density matrix; in open shell,
      ! it is the alpha density.
      call g2g_timer_sum_start('SCF acceleration')
      call converger_fock(niter, M_f, fock_aop, 1, NCOa_f, HL_gap, Xmat)
      call g2g_timer_sum_pause('SCF acceleration')

      ! Fock(ON) diagonalization
      if ( allocated(morb_coefon) ) deallocate(morb_coefon)
      allocate( morb_coefon(M_f,M_f) )
      call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')

      call fock_aop%Diagon_datamat( morb_coefon, morb_energy )
      call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')

      ! Base change of coeficients ( (X^-1)*C ) and construction of new
      ! density matrix.
      call g2g_timer_sum_start('SCF - MOC base change (sum)')
      call Xmat%multiply(morb_coefat, morb_coefon)
      call standard_coefs( morb_coefat )
      call g2g_timer_sum_pause('SCF - MOC base change (sum)')

      if ( allocated(morb_coefon) ) deallocate(morb_coefon)
      call rho_aop%Dens_build(M_f, NCOa_f, ocupF, morb_coefat)
      call rho_aop%Gets_data_AO(rho_a)
      call messup_densmat(rho_a)

      Eorbs      = morb_energy
      MO_coef_at = morb_coefat

      if (OPEN) then
         ! In open shell, performs the previous operations for beta operators.
         call g2g_timer_sum_start('SCF acceleration')
         call converger_fock(niter, M_f, fock_bop, 2, NCOb_f, HL_gap, Xmat)
         call g2g_timer_sum_pause('SCF acceleration')

         ! Fock(ON) diagonalization
         if ( allocated(morb_coefon) ) deallocate(morb_coefon)
         allocate( morb_coefon(M_f,M_f) )

         call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
         call fock_bop%Diagon_datamat( morb_coefon, morb_energy )
         call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')

         ! Base change of coeficients ( (X^-1)*C ) and construction of new
         ! density matrix.
         call g2g_timer_sum_start('SCF - MOC base change (sum)')
         call Xmat%multiply(morb_coefat, morb_coefon)
         call standard_coefs( morb_coefat )
         call g2g_timer_sum_pause('SCF - MOC base change (sum)')

         if ( allocated(morb_coefon) ) deallocate(morb_coefon)
         call rho_bop%Dens_build(M_f, NCOb_f, ocupF, morb_coefat)
         call rho_bop%Gets_data_AO(rho_b)
         call messup_densmat( rho_b )

         Eorbs_b      = morb_energy
         MO_coef_at_b = morb_coefat
      endif

      ! Calculates HOMO-LUMO gap.
      HL_gap = abs(Eorbs(NCOa_f+1) - Eorbs(NCOa_f))
      if (OPEN) HL_gap = abs(min(Eorbs(NCOa_f+1),Eorbs_b(NCOb_f+1)) &
                             - max(Eorbs(NCOa_f),Eorbs_b(NCOb_f)))

      ! We are not sure how to translate the sumation over molecular orbitals
      ! and energies when changing from TBDFT system to DFT subsystem. Forces
      ! may be broken due to this. This should not be affecting normal DFT
      ! calculations.

      ! X is DEPRECATED
      do ii=1,M
      do jj=1,M
         X( ii, 2*M+jj ) = morb_coefat( MTB+ii, jj )
      enddo
      enddo

      ! Perfoms TBDFT checks and extracts density matrices. Allocates xnano,
      ! which contains the total (alpha+beta) density matrix.
      allocate ( xnano(M,M) )

      if (tbdft_calc == 0) then
         xnano = rho_a
         if (OPEN) xnano = xnano + rho_b
      else
         rhoa_TBDFT = rho_a
         call extract_rhoDFT(M, rho_a, rho_a0)
         xnano = rho_a0

         if (OPEN) then
            rhob_TBDFT = rho_b
            call extract_rhoDFT(M, rho_b, rho_b0)
            xnano = xnano + rho_b0
         endif
      endif

      if ((rho_LS > 1) .and. (niter > 10)) then
         ! Performs a linear search in Rho if activated. This uses the
         ! vector-form densities as the old densities, and matrix-form
         ! densities as the new ones.
         if (open) then
            call do_rho_ls(En, E1, E2, Exc, xnano, Pmat_vec, Hmat_vec,    &
                           Fmat_vec, Fmat_vec2, Gmat_vec, Ginv_vec, memo, &
                           rho_a, rho_b, rhoalpha, rhobeta)
         else
            call do_rho_ls(En, E1, E2, Exc, xnano, Pmat_vec, Hmat_vec, &
                           Fmat_vec, Fmat_vec2, Gmat_vec, Ginv_vec, memo)
         endif
      endif

      E = E + Eexact
      ! Checks convergence criteria and starts linear search if able.
      call converger_check(Pmat_vec, xnano, Evieja, E, niter, converged, &
                           open, changed_to_LS)

      ! Updates old density matrices with the new ones and updates energy.
      call sprepack('L', M, Pmat_vec, xnano)
      if (OPEN) then
         if (tbdft_calc /= 0) then
            call sprepack('L', M, rhoalpha, rho_a0)
            call sprepack('L', M, rhobeta , rho_b0)
         else
            call sprepack('L', M, rhoalpha, rho_a)
            call sprepack('L', M, rhobeta , rho_b)
         endif
      endif
      deallocate ( xnano )
      Evieja = E

      call g2g_timer_stop('Total iter')
      call g2g_timer_sum_pause('Iteration')
999 continue

   call g2g_timer_sum_start('Finalize SCF')

   ! Checks of convergence
   if (niter >= nMax) then
      call write_final_convergence(.false., nMax, Evieja)
      noconverge = noconverge + 1
      converge   = 0
   else
      call write_final_convergence(.true., niter, Evieja)
      converge   = converge + 1
      noconverge = 0
   endif

   if (changed_to_LS) then
      changed_to_LS = .false.
      nMax          = nMax / 2
      Rho_LS        = 1
   endif

   if (noconverge > 4) then
      write(6,'(A)') "FATAL ERROR - No convergence achieved "&
                    &"4 consecutive times."
      stop
   endif


   if (MOD(npas,energy_freq).eq.0) then
!       Resolve with last density to get XC energy
        call g2g_timer_sum_start('Exchange-correlation energy')
        call g2g_new_grid(igrid)
        call g2g_solve_groups(1, Exc, 0)
        call g2g_timer_sum_stop('Exchange-correlation energy')

!       COmputing the QM/MM contribution to total energy
!       Total SCF energy =
!       E1   - kinetic + nuclear attraction + QM/MM interaction + effective
!              core potential
!       E2   - Coulomb
!       En   - nuclear-nuclear repulsion
!       Ens  - MM point charge - nuclear interaction
!       Exc  - exchange-correlation
!       Eecp - Efective core potential
!       E_restrain - distance restrain

!       NucleusQM-CHarges MM
        Es=Ens

!       One electron Kinetic (with aint >3) or Kinetic + Nuc-elec (aint >=3)
        call int1(En, Fmat_vec, Hmat_vec, Smat, d, r, Iz, natom, &
                  ntatom)

!       Computing the E1-fock without the MM atoms
        if (nsol.gt.0.and.igpu.ge.1) then
          call aint_qmmm_init(0,r,pc)
          call aint_qmmm_fock(E1s,Etrash)
          call aint_qmmm_init(nsol,r,pc)
        endif

!       E1s (here) is the 1e-energy without the MM contribution
        E1s=0.D0
        do kk=1,MM
          E1s = E1s + Pmat_vec(kk) * Hmat_vec(kk)
        enddo


!       Es is the QM/MM energy computated as total 1e - E1s + QMnuc-MMcharges
        Es=Es+E1-E1s

        ! Calculates DTFD3 Grimme's corrections to energy.
        call g2g_timer_sum_start("DFTD3 Energy")
        call dftd3_energy(E_dftd, d, natom, .true.)
        call g2g_timer_sum_pause("DFTD3 Energy")

!       Exact Exchange Energy PBE0
        Eexact = 0.0d0
        if ( PBE0 ) then
           call g2g_timer_sum_start("Exact Exchange Energy")
           do ii=1,M
             Eexact = Eexact + 0.5D0 * rho_a0(ii,ii) * FockEE_a0(ii,ii)
             do jj=1,ii-1
               Eexact = Eexact + 0.5D0 * rho_a0(ii,jj) * FockEE_a0(ii,jj)
               Eexact = Eexact + 0.5D0 * rho_a0(jj,ii) * FockEE_a0(jj,ii)
             enddo
           enddo
           if ( OPEN ) then
              do ii=1,M
                Eexact = Eexact + 0.5D0 * rho_b0(ii,ii) * FockEE_b0(ii,ii)
                do jj=1,ii-1
                  Eexact = Eexact + 0.5D0 * rho_b0(ii,jj) * FockEE_b0(ii,jj)
                  Eexact = Eexact + 0.5D0 * rho_b0(jj,ii) * FockEE_b0(jj,ii)
                enddo
              enddo
           endif
           Eexact = Eexact * (-0.25d0)
           call g2g_timer_sum_pause("Exact Exchange Energy")
        endif

!       Part of the QM/MM contrubution are in E1
        E=E1+E2+En+Ens+Exc+E_restrain+E_dftd+Eexact

!       Write Energy Contributions
        if (npas.eq.1) npasw = 0

        if (npas.gt.npasw) then
           call ECP_energy( MM, Pmat_vec, Eecp, Es )
           call write_energies(E1, E2, En, Ens, Eecp, Exc, ecpmode, E_restrain,&
                               number_restr, nsol, E_dftd, Eexact)
           npasw=npas+10
        end if
      endif ! npas



      ! Calculation of energy weighted density matrix
      call g2g_timer_sum_start('energy-weighted density')
      kkk = 0
      Pmat_en_wgt = 0.0D0
      if (.not. OPEN) then
         ! Closed shell
         do jj = MTB+1, MTB+M
            kkk = kkk +1
            do kk = MTB+1, NCOa_f
               Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 2.0D0 * Eorbs(kk) * &
                                  MO_coef_at(jj,kk) * MO_coef_at(jj,kk)
            enddo

            do ii = MTB+jj+1, M_f
               kkk = kkk +1
               do kk = MTB+1, NCOa_f
                  Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 4.0D0 * Eorbs(kk) * &
                                     MO_coef_at(ii,kk) * MO_coef_at(jj,kk)
               enddo
            enddo
         enddo

      else
         ! Open shell
         do jj = MTB+1, MTB+M
            kkk = kkk +1
            do kk = MTB+1, NCOa_f
               Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - Eorbs(kk) * &
                                  MO_coef_at(jj,kk) * MO_coef_at(jj,kk)
            enddo
            do kk = MTB+1, NCOb_f
               Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - Eorbs_b(kk) * &
                                  MO_coef_at_b(jj,kk) * MO_coef_at_b(jj,kk)
            enddo

            do ii = MTB+jj+1, MTB+M
               kkk = kkk +1
               do kk = MTB+1, NCOa_f
                  Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 2.0D0 * Eorbs(kk) * &
                                     MO_coef_at(ii,kk) * MO_coef_at(jj,kk)
               enddo
               do kk = MTB+1, NCOb_f
                  Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 2.0D0 * Eorbs_b(kk) * &
                                     MO_coef_at_b(ii,kk) * MO_coef_at_b(jj,kk)
               enddo
            enddo
         enddo
      endif

      call g2g_timer_sum_stop('energy-weighted density')

      call cubegen_write(MO_coef_at(MTB+1:MTB+M,1:M))

   if (gaussian_convert) then       ! Density matrix translation from Gaussian09
      allocate(rho_exc(M,M))
      call translation(M,rho_exc)   ! Reorganizes Rho to LIO format.

      do jj=1,M                     ! Stores matrix in vector form.
         Pmat_vec(jj + (2*M-jj)*(jj-1)/2) = rho_exc(jj,jj)
         do kk = jj+1, M
            Pmat_vec( kk + (2*M-jj)*(jj-1)/2) = rho_exc(jj,kk) * 2.0D0
         enddo
      enddo

      deallocate(rho_exc)
   endif                            ! End of translation

!  Excited States routines
   call ExcProp(MO_coef_at,MO_coef_at_b,Eorbs,Eorbs_b,E)

!------------------------------------------------------------------------------!
! TODO: have ehrendyn call SCF and have SCF always save the resulting rho in
!       a module so that ehrendyn can retrieve it afterwards.
!       Remove all of this.
!
      if (doing_ehrenfest) then
         call spunpack('L',M,Pmat_vec,RealRho)
         call fix_densmat(RealRho)
         call ehrendyn_init(natom, M, RealRho)
      endif


!------------------------------------------------------------------------------!
! TODO: Deallocation of variables that should be removed
! TODO: MEMO should be handled differently...

      if (MEMO) then
        deallocate(kkind,kkinds)
        deallocate(cool,cools)
      endif

!------------------------------------------------------------------------------!
! MovieMaker
      call spunpack('L',M,Pmat_vec,RealRho)
      call fix_densmat(RealRho)
      call movieprint( natom, M, npas-1, Iz, r, dcmplx( RealRho ) )


      call Xmat%destroy()
      call Ymat%destroy()

      call g2g_timer_stop('SCF')
      call g2g_timer_sum_stop('Finalize SCF')
      call g2g_timer_sum_stop('SCF')
      call g2g_timer_stop('SCF_full')
      end subroutine SCF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
