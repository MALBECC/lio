!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine SCF(E)
!
!------------------------------------------------------------------------------!
! DIRECT VERSION
! Calls all integrals generator subroutines : 1 el integrals,
! 2 el integrals, exchange fitting , so it gets S matrix, F matrix
! and P matrix in lower storage mode (symmetric matrices)
!
! Dario Estrin, 1992
!
!------------------------------------------------------------------------------!
! Modified to f90
!
! Nick, 2017
!
!------------------------------------------------------------------------------!
! Header with new format. Added comments on how to proceed with a cleanup of
! of the subroutines. Other things to do:
! TODO: change to 3 space indentation.
! TODO: break at line 80.
! TODO: change to lowercase (ex: implicit none)
!
! This log can be removed once all improvements have been made.
!
! FFR, 01/2018
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use ehrensubs , only: ehrendyn_init
   use garcha_mod, only : NCO, natom, Nang, number_restr, hybrid_converg, MEMO,&
                          igrid, energy_freq, converge, noconverge, lowdin,    &
                          cubegen_only, VCINP, primera, Nunp, GOLD, igrid2,    &
                          predcoef, nsol, r, pc, DIIS, told, Etold, Enucl, Iz, &
                          Eorbs, NMAX,Dbug, doing_ehrenfest, first_step,       &
                          total_time, MO_coef_at, MO_coef_at_b, Smat, good_cut,&
                          ndiis, rhoalpha, rhobeta, OPEN, RealRho, d, ntatom,  &
                          Eorbs_b, npas, X, npasw, Fmat_vec, Fmat_vec2,        &
                          Ginv_vec, Gmat_vec, Hmat_vec, Pmat_en_wgt, Pmat_vec
   use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC, &
                       FOCK_ECP_read,FOCK_ECP_write,IzECP
   use field_data, only: field, fx, fy, fz
   use field_subs, only: field_calc, field_setup_old
   use td_data, only: timedep, tdrestart, tdstep
   use transport_data, only : generate_rho0
   use time_dependent, only : TD
   use faint_cpu, only: int1, intsol, int2, int3mem, int3lu
   use tbdft_data, only : tbdft_calc, MTBDFT, MTB, chargeA_TB, chargeB_TB,     &
                         rhoa_tbdft, rhob_tbdft
   use tbdft_subs, only : tbdft_init, getXY_TBDFT, build_chimera_TBDFT,        &
                          extract_rhoDFT, construct_rhoTBDFT, tbdft_scf_output
   use cubegen       , only: cubegen_vecin, cubegen_matin, cubegen_write
   use mask_ecp      , only: ECP_init, ECP_fock, ECP_energy
   use typedef_sop   , only: sop              ! Testing SOP
   use fockbias_subs , only: fockbias_loads, fockbias_setmat, fockbias_apply
   use SCF_aux       , only: neighbour_list_2e, seek_nan, standard_coefs, &
                             messup_densmat, fix_densmat
   use liosubs_math  , only: transform
   use linear_algebra, only: matrix_diagon
   use converger_subs, only: converger_init, conver
   use mask_cublas   , only: cublas_setmat, cublas_release
   use typedef_operator, only: operator !Testing operator
   use trans_Data    , only: gaussian_convert, rho_exc, translation
#  ifdef  CUBLAS
      use cublasmath , only: cumxp_r
#  endif
   use initial_guess_subs, only: get_initial_guess
   use fileio       , only: write_energies, write_energy_convergence, &
                            write_final_convergence
   use fileio_data  , only: verbose
   use basis_data   , only: kkinds, kkind, cools, cool, Nuc, nshell, ncont, a, &
                            c, M, Md
   use lr_data, only: lresp
   use lrtddft, only: linear_response
   use converger_ls , only: Rho_LS, changed_to_LS, P_conver, P_linearsearch_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   implicit none
   integer :: Nel
   integer :: niter
   real*8  :: sq2
   real*8  :: good
   real*8  :: del
   real*8  :: DAMP0
   real*8  :: DAMP
   integer :: igpu

!  The following two variables are in a part of the code that is never
!  used. Check if these must be taken out...
   real*8  :: factor
   integer :: IDAMP

   real*8, allocatable :: rho_test(:,:)
   real*8, allocatable :: fockat(:,:)
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
   real*8, allocatable :: Xmat(:,:)
   real*8, allocatable :: Ymat(:,:)
   real*8, allocatable :: sqsmat(:,:)
   real*8, allocatable :: tmpmat(:,:)

   real*8              :: dipxyz(3)

! FIELD variables (maybe temporary)
   real*8  :: Qc, Qc2, g
   integer :: ng2

!------------------------------------------------------------------------------!
! Energy contributions and convergence
   real*8, intent(inout) :: E ! Total SCF energy

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
   real*8 :: Ex          ! exchange-correlation inside SCF loop
   real*8 :: Exc         ! exchange-correlation
   real*8 :: Etrash      ! auxiliar variable
   real*8 :: Egood       !
   real*8 :: Evieja      !


! CUBLAS
!------------------------------------------------------------------------------!
   integer*8          :: dev_Xmat, dev_Ymat



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! TODO : Variables to eliminate...
   real*8, allocatable :: xnano(:,:)
   integer :: MM, MM2, MMd, Md2
   integer :: M1, M2

   real*8, allocatable :: Y(:,:)
   real*8, allocatable :: Ytrans(:,:)
   real*8, allocatable :: Xtrans(:,:)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!carlos: Operators for matrices with alpha and beta spins.
   type(operator)      :: rho_aop, fock_aop
   type(operator)      :: rho_bop, fock_bop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!carlos: Open shell, variables.
   real*8              :: ocupF
   integer             :: NCOa, NCOb

! LINSEARCH
   integer :: nniter
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   call g2g_timer_start('SCF_full')

   changed_to_LS=.false. ! LINSEARCH

   if (verbose > 1) then
      write(*,*)
      write(*,'(A)') "Starting SCF cycles."
   endif

!------------------------------------------------------------------------------!
!TBDFT: initialisation of variable that could depend of TB.
   allocate (fock_a0(M,M), rho_a0(M,M))

   if (tbdft_calc) then
      call tbdft_init(M, Nuc,natom,OPEN)
      M_f = MTBDFT
   else
      M_f = M
   end if

   allocate(fock_a(M_f,M_f), rho_a(M_f,M_f))
   allocate(morb_energy(M_f), morb_coefat(M_f,M_f))
   if (OPEN) then
      allocate(fock_b(M_f,M_f), rho_b(M_f,M_f))
   end if


!------------------------------------------------------------------------------!
   call ECP_init()

   call g2g_timer_start('SCF')
   call g2g_timer_sum_start('SCF')
   call g2g_timer_sum_start('Initialize SCF')

   E=0.0D0
   E1=0.0D0
   En=0.0D0
   E2=0.0D0
   Es=0.0D0
   Eecp=0.d0
   Ens=0.0D0
   E_restrain=0.d0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%    Distance Restrain     %%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        IF (number_restr.GT.0) THEN
          call get_restrain_energy(E_restrain)
          WRITE(*,*) "DISTANCE RESTRAIN ADDED TO FORCES"
        END IF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


      sq2=sqrt(2.D0)
      MM=M*(M+1)/2
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M

!carlos: ocupation factor
      ocupF = 2.0d0
!carlos:NCOa works in open shell and close shell
      NCOa=NCO
      if (OPEN) then
!Number of OM down
         NCOb=NCO+Nunp
         ocupF=1.0d0
         allocate(rho_b0(M,M),fock_b0(M,M))
      end if

      M1=1 ! first P

!------------------------------------------------------------------------------!
! TODO: I don't like ending timers inside a conditional...
       if (cubegen_only) then
          call cubegen_vecin( M, MO_coef_at )
          call g2g_timer_sum_stop('Initialize SCF')
          call g2g_timer_sum_stop('SCF')
          return
       end if


!------------------------------------------------------------------------------!
! TODO: damp and gold should no longer be here??
! TODO: Qc should probably be a separated subroutine? Apparently it is only
!       used in dipole calculation so...it only adds noise to have it here.
! TODO: convergence criteria should be set at namelist/keywords setting

      Nel=2*NCO+Nunp
!
      good=1.00D0
      Egood=1.00D0
      Evieja=0.d0

      niter=0
      DAMP0=GOLD
      DAMP=DAMP0

      Qc=0.0D0
      do ii=1,natom
         Qc=Qc+Iz(ii)
      enddo
      Qc=Qc-Nel
      Qc2=Qc**2


!------------------------------------------------------------------------------!
! TODO: this whole part which calculates the non-electron depending terms of
!       fock and the overlap matrix should probably be in a separated sub.
!       (diagonalization of overlap, starting guess, the call to TD, should be taken out)
!
! Reformat from here...

! Nano: calculating neighbour list helps to make 2 electrons integral scale
! linearly with natoms/basis
!
      call neighbour_list_2e(natom, ntatom, r, d)

! Goes straight to TD if a restart is used.
      if ((timedep.eq.1).and.(tdrestart)) then
        call g2g_timer_sum_stop('Initialize SCF')
        call g2g_timer_sum_start('TD')
        if(OPEN) then
           call TD(fock_aop, rho_aop, fock_bop, rho_bop)
        else
           call TD(fock_aop, rho_aop)
        endif
        call g2g_timer_sum_stop('TD')
        return
      endif
!
! -Create integration grid for XC here
! -Assign points to groups (spheres/cubes)
! -Assign significant functions to groups
! -Calculate point weights
!
      call g2g_timer_sum_start('Exchange-correlation grid setup')
      call g2g_reload_atom_positions(igrid2)
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

        if ( allocated(Xmat) ) deallocate(Xmat)
        if ( allocated(Ymat) ) deallocate(Ymat)
        allocate(Xmat(M_f,M_f), Ymat(M_f,M_f))

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
        if ( allocated(sqsmat) ) deallocate(sqsmat)
        if ( allocated(tmpmat) ) deallocate(tmpmat)
        allocate( sqsmat(M,M), tmpmat(M,M) )
        call overop%Gets_orthog_2m( 2, 0.0d0, tmpmat, sqsmat )
        call fockbias_loads( natom, nuc )
        call fockbias_setmat( sqsmat )
        deallocate( sqsmat, tmpmat )


!TBDFT: Dimensions of Xmat and Ymat are modified for TBDFT.
!
! TODO: this is nasty, a temporary solution would be to have a Msize variable
!       be assigned M (or, even better, "basis_size") or MTBDFT
!       ("basis_size_dftb") according to the case

       if (tbdft_calc) then
          call getXY_TBDFT(M,X_min,Y_min,Xmat,Ymat)
       else
          do jj = 1, M
          do ii = 1, M
             Xmat(ii,jj) = X_min(ii,jj)
             Ymat(ii,jj) = Y_min(ii,jj)
          enddo
          enddo
      end if

! CUBLAS
   call cublas_setmat( M_f, Xmat, dev_Xmat)
   call cublas_setmat( M_f, Ymat, dev_Ymat)


! Generates starting guess
!
   if ( (.not.VCINP) .and. primera ) then
      call get_initial_guess(M, MM, NCO, NCOb, Xmat(MTB+1:MTB+M,MTB+1:MTB+M),  &
                             Hmat_vec, Pmat_vec, rhoalpha, rhobeta, OPEN, &
                             natom, Iz, nshell, Nuc)
      primera = .false.
   end if

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

      if (hybrid_converg) DIIS=.true. ! cambio para convergencia damping-diis
      call g2g_timer_sum_stop('Initialize SCF')

!------------------------------------------------------------------------------!
! TODO: Maybe evaluate conditions for loop continuance at the end of loop
!       and condense in a single "keep_iterating" or something like that.
      do 999 while ((good.ge.told.or.Egood.ge.Etold).and.niter.le.NMAX)
        call g2g_timer_start('Total iter')
        call g2g_timer_sum_start('Iteration')
        call g2g_timer_sum_start('Fock integrals')

        niter=niter+1
        nniter=niter
        IF (changed_to_LS .and. niter.eq. (NMAX/2 +1)) nniter=1 
!       (first steep of damping after NMAX steeps without convergence)

        E1=0.0D0

!------------------------------------------------------------------------------!
!       Fit density basis to current MO coeff and calculate Coulomb F elements
!
! TODO: Calculation of fock terms should be separated. Also, maybe it would be
!       convenient to add the NaN checks and the timers inside the subroutine
!       calls is a more systematic way.

!       Test for NaN
        if (Dbug) call SEEK_NaN(Pmat_vec,1,MM,"RHO Start")
        if (Dbug) call SEEK_NaN(Fmat_vec,1,MM,"FOCK Start")

!       Computes Coulomb part of Fock, and energy on E2
        call g2g_timer_sum_start('Coulomb fit + Fock')
        call int3lu(E2, Pmat_vec, Fmat_vec2, Fmat_vec, Gmat_vec, Ginv_vec, &
                    Hmat_vec, open, MEMO)
        call g2g_timer_sum_pause('Coulomb fit + Fock')

!       Test for NaN
        if (Dbug) call SEEK_NaN(Pmat_vec,1,MM,"RHO Coulomb")
        if (Dbug) call SEEK_NaN(Fmat_vec,1,MM,"FOCK Coulomb")

!       XC integration / Fock elements
        call g2g_timer_sum_start('Exchange-correlation Fock')
        call g2g_solve_groups(0,Ex,0)
        call g2g_timer_sum_pause('Exchange-correlation Fock')
!       Test for NaN
        if (Dbug) call SEEK_NaN(Pmat_vec,1,MM,"RHO Ex-Corr")
        if (Dbug) call SEEK_NaN(Fmat_vec,1,MM,"FOCK Ex-Corr")


!------------------------------------------------------------------------------!
! REACTION FIELD CASE
!
! TODO: what is reaction field?? in any case, it should be either in its own
!       module and only appear here as a subroutine, or be included in the
!       separated fock calculation subroutine...
!
        if ( generate_rho0 ) then
           if (field) call field_setup_old(1.0D0, 0, fx, fy, fz)
           call field_calc(E1, 0.0D0, Pmat_vec(1:MM), Fmat_vec2, Fmat_vec, &
                           r, d, Iz, natom, ntatom, open)

           do kk=1,MM
               E1 = E1 + Pmat_vec(kk) * Hmat_vec(kk)
           enddo
        else
!          E1 includes solvent 1 electron contributions
           do kk=1,MM
              E1 = E1 + Pmat_vec(kk) * Hmat_vec(kk)
           enddo

        endif
        call g2g_timer_sum_pause('Fock integrals')

!------------------------------------------------------------------------------!
! TBDFT: we extract rho and fock before conver routine
!carlos: extractions for Open Shell and Close Shell.
        if (OPEN) then
           call spunpack_rho('L',M,rhoalpha,rho_a0)
           call spunpack('L', M, Fmat_vec, fock_a0)
           call spunpack_rho('L',M,rhobeta,rho_b0)
           call spunpack('L', M, Fmat_vec2, fock_b0)
           call fockbias_apply( 0.0d0, fock_a0)
           call fockbias_apply( 0.0d0, fock_b0)
        else
           call spunpack_rho('L',M,Pmat_vec,rho_a0)
           call spunpack('L', M, Fmat_vec, fock_a0)
           call fockbias_apply( 0.0d0, fock_a0 )
        end if

!------------------------------------------------------------------------------!
! TBDFT: Fock and Rho for TBDFT are builded.
!
! TODO: this should be wrapped inside a single dftb subroutine. Also, two
!       consecutive tbdft_calc switches? really?
!
      if (tbdft_calc) then
         NCOa_f = NCOa + MTB
         call build_chimera_TBDFT (M, fock_a0, fock_a, natom)
         call construct_rhoTBDFT(M, rho_a, rho_a0 ,rhoa_tbdft, niter,OPEN)
         if (OPEN) then
            NCOb_f = NCOb + MTB
            call build_chimera_TBDFT(M, fock_b0, fock_b, natom)
            call construct_rhoTBDFT(M, rho_b, rho_b0 ,rhob_tbdft,niter, OPEN)
         end if
      else
         NCOa_f = NCOa
         fock_a=fock_a0
         rho_a=rho_a0
         if (OPEN) then
            NCOb_f = NCOb
            fock_b=fock_b0
            rho_b=rho_b0
         end if
      endif
!carlos: storing rho and fock in operator.

   call rho_aop%Sets_data_AO(rho_a)
   call fock_aop%Sets_data_AO(fock_a)

   if (OPEN) then
      call rho_bop%Sets_data_AO(rho_b)
      call fock_bop%Sets_data_AO(fock_b)
   end if
!------------------------------------------------------------------------------!
!  Convergence accelerator processing
        call g2g_timer_sum_start('SCF acceleration')
        if (niter==1) then
           call converger_init( M_f, ndiis, DAMP, DIIS, hybrid_converg, OPEN )
        end if
!carlos: this is repeted twice for open shell

!%%%%%%%%%%%%%%%%%%%%
!CLOSE SHELL OPTION |
!%%%%%%%%%%%%%%%%%%%%
#       ifdef CUBLAS
           call conver(nniter, good, good_cut, M_f, rho_aop, fock_aop,         &
                       dev_Xmat, dev_Ymat, 1)
#       else
           call conver(nniter, good, good_cut, M_f, rho_aop, fock_aop, Xmat,   &
                       Ymat, 1)
#       endif

        call g2g_timer_sum_pause('SCF acceleration')
!------------------------------------------------------------------------------!
!  Fock(ON) diagonalization
        if ( allocated(morb_coefon) ) deallocate(morb_coefon)
        allocate( morb_coefon(M_f,M_f) )
        call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
        call fock_aop%Diagon_datamat( morb_coefon, morb_energy )
        call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')
!
!
!------------------------------------------------------------------------------!
!  Base change of coeficients ( (X^-1)*C ) and construction of new density
!  matrix
        call g2g_timer_sum_start('SCF - MOC base change (sum)')
#       ifdef CUBLAS
           call cumxp_r( morb_coefon, dev_Xmat, morb_coefat, M_f)
#       else
           morb_coefat = matmul( Xmat, morb_coefon )
#       endif
        call standard_coefs( morb_coefat )
        call g2g_timer_sum_pause('SCF - MOC base change (sum)')

        if ( allocated(morb_coefon) ) deallocate(morb_coefon)
        call rho_aop%Dens_build(M_f, NCOa_f, ocupF, morb_coefat)
        call rho_aop%Gets_data_AO(rho_a)
        call messup_densmat( rho_a )

        Eorbs      = morb_energy
        MO_coef_at = morb_coefat
        
    if (OPEN) then
!%%%%%%%%%%%%%%%%%%%%
!OPEN SHELL OPTION  |
!%%%%%%%%%%%%%%%%%%%%
        call g2g_timer_sum_start('SCF acceleration')
#       ifdef CUBLAS
           call conver(nniter, good, good_cut, M_f, rho_bop, fock_bop,         &
                       dev_Xmat, dev_Ymat, 2)
#       else
           call conver(nniter, good, good_cut, M_f, rho_bop, fock_bop, Xmat,     &
                       Ymat, 2)
#       endif

        call g2g_timer_sum_pause('SCF acceleration')

!------------------------------------------------------------------------------!
!  Fock(ON) diagonalization
        if ( allocated(morb_coefon) ) deallocate(morb_coefon)
        allocate( morb_coefon(M_f,M_f) )

        call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
        call fock_bop%Diagon_datamat( morb_coefon, morb_energy )
        call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')
!
!
!------------------------------------------------------------------------------!
!  Base change of coeficients ( (X^-1)*C ) and construction of new density
!  matrix
        call g2g_timer_sum_start('SCF - MOC base change (sum)')
#       ifdef CUBLAS
           call cumxp_r( morb_coefon, dev_Xmat, morb_coefat, M_f)
#       else
           morb_coefat = matmul( Xmat, morb_coefon )
#       endif
        call standard_coefs( morb_coefat )
        call g2g_timer_sum_pause('SCF - MOC base change (sum)')

        if ( allocated(morb_coefon) ) deallocate(morb_coefon)
        call rho_bop%Dens_build(M_f, NCOb_f, ocupF, morb_coefat)
        call rho_bop%Gets_data_AO(rho_b)
        call messup_densmat( rho_b )

        Eorbs_b      = morb_energy
        MO_coef_at_b = morb_coefat
        
    end if!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------!
!carlos: storing matrices
!------------------------------------------------------------------------------!
!  We are not sure how to translate the sumation over molecular orbitals
!  and energies when changing from TBDFT system to DFT subsystem. Forces
!  may be broken due to this.
!  This should not be affecting normal DFT calculations.

!carlos: this is not working with openshell. morb_coefat is not depending of
!        the spin factor.
        do ii=1,M
          do jj=1,M
            X( ii, M2+jj ) = morb_coefat( MTB+ii, jj )
          enddo
        enddo

!------------------------------------------------------------------------------!
! carlos: added to separate from rho the DFT part
!
! TODO: again, this should be handled differently...
! TODO: make xnano go away...only remains here
!
        allocate ( xnano(M,M) )

        if (tbdft_calc) then
          rhoa_TBDFT = rho_a
          call extract_rhoDFT(M, rho_a, rho_a0)

          if (OPEN) then
              rhob_TBDFT = rho_b
              call extract_rhoDFT(M, rho_b, rho_b0)
              call sprepack('L',M,rhoalpha,rho_a0)
              call sprepack('L',M,rhobeta,rho_b0)
              xnano=rho_a0+rho_b0
          else
              xnano=rho_a0
          end if

        else

          if (OPEN) then
             call sprepack('L',M,rhoalpha,rho_a)
             call sprepack('L',M,rhobeta,rho_b)
             xnano=rho_a+rho_b
          else
              xnano=rho_a
          end if
        end if
!------------------------------------------------------------------------------!
! Convergence criteria and lineal search in P

       IF (OPEN) call P_conver(nniter, En, E1, E2, Ex, good, xnano, rho_a, rho_b)
       IF (.not. OPEN) call P_conver(nniter, En, E1, E2, Ex, good, xnano, rho_a, rho_a)
!------------------------------------------------------------------------------!

      deallocate ( xnano )

!------------------------------------------------------------------------------!
! TODO: finalization of the loop is a little bit messy. Also: "999 continue"??
!       I think it is time we regularized this loop...

        ! Damping factor update
        DAMP=DAMP0
        E=E1+E2+En

!       write energy at every step
        if (niter.eq.NMAX) then
           if (Rho_LS .eq.0) then
              write(6,*) 'NO CONVERGENCE AT ',NMAX,' ITERATIONS'
              write(6,*) 'trying Lineal search'
              Rho_LS=1
              NMAX=2*NMAX
              changed_to_LS=.true.
              call P_linearsearch_init()
           end if
        end if

        Egood=abs(E+Ex-Evieja)
        Evieja=E+Ex

        ! Write energy at every step
        call write_energy_convergence(niter, Evieja, good, told, egood, etold)

        call g2g_timer_stop('Total iter')
        call g2g_timer_sum_pause('Iteration')

 999  continue

      call g2g_timer_sum_start('Finalize SCF')

!------------------------------------------------------------------------------!
!     Checks of convergence
!
      if (niter.ge.NMAX) then
         call write_final_convergence(.false., NMAX, Evieja)
         noconverge = noconverge + 1
         converge   = 0
      else
         call write_final_convergence(.true., niter, Evieja)
         noconverge = 0
         converge   = converge + 1
      endif

      if (changed_to_LS) then
         changed_to_LS=.false.
         NMAX=NMAX/2
         Rho_LS=0
      end if

      if (noconverge.gt.4) then
         write(6,'(A)')  'FATAL ERROR - No convergence achieved 4 times.'
         stop
      endif
!------------------------------------------------------------------------------!

!TBDFT: Mulliken analysis of TB part
   if (tbdft_calc) call tbdft_scf_output(M,OPEN)

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

!       Part of the QM/MM contrubution are in E1
        E=E1+E2+En+Ens+Exc+E_restrain



!       Write Energy Contributions
        if (npas.eq.1) npasw = 0

        if (npas.gt.npasw) then
           call ECP_energy( MM, Pmat_vec, Eecp, Es )
           call write_energies(E1, E2, En, Ens, Eecp, Exc, ecpmode, E_restrain,&
                               number_restr, nsol)
           npasw=npas+10
        end if
      endif ! npas



      ! Calculation of energy weighted density matrix
      call g2g_timer_sum_start('energy-weighted density')
      kkk = 0
      Pmat_en_wgt = 0.0D0
      if (.not. OPEN) then
         ! Closed shell
         do jj = 1 , M_f
            kkk = kkk +1
            do kk = 1, NCOa_f
               Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 2.0D0 * Eorbs(kk) * &
                                  MO_coef_at(jj,kk) * MO_coef_at(jj,kk)
            enddo

            do ii = jj+1, M_f
               kkk = kkk +1
               do kk = 1, NCOa_f
                  Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 4.0D0 * Eorbs(kk) * &
                                     MO_coef_at(ii,kk) * MO_coef_at(jj,kk)
               enddo
            enddo
         enddo
         
      else
         ! Open shell
         do jj = 1 , M_f
            kkk = kkk +1
            do kk = 1, NCOa_f
               Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - Eorbs(kk) * &
                                  MO_coef_at(jj,kk) * MO_coef_at(jj,kk)
            enddo
            do kk = 1, NCOb_f
               Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - Eorbs_b(kk) * &
                                  MO_coef_at_b(jj,kk) * MO_coef_at_b(jj,kk)
            enddo

            do ii = jj+1, M_f
               kkk = kkk +1
               do kk=1, NCOa_f
                  Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 2.0D0 * Eorbs(kk) * &
                                     MO_coef_at(ii,kk) * MO_coef_at(jj,kk)
               enddo
               do kk=1, NCOb_f
                  Pmat_en_wgt(kkk) = Pmat_en_wgt(kkk) - 2.0D0 * Eorbs_b(kk) * &
                                     MO_coef_at_b(ii,kk) * MO_coef_at_b(jj,kk)
               enddo
            enddo
         enddo
      endif

      call g2g_timer_sum_stop('energy-weighted density')

      call cubegen_matin( M, X )

   if (gaussian_convert) then       ! Density matrix translation from Gaussian09
      allocate(rho_exc(M,M))
      call translation(M,rho_exc)   ! Reorganizes Rho to LIO format.

      do jj=1,M                     ! Stores matrix in vector form.
         Pmat_vec(jj + (M2-jj)*(jj-1)/2) = rho_exc(jj,jj)
         do kk = jj+1, M
            Pmat_vec( kk + (M2-jj)*(jj-1)/2) = rho_exc(jj,kk) * 2.0D0
         enddo
      enddo

      deallocate(rho_exc)
   endif                            ! End of translation

   if (lresp) then
     if (OPEN) then
       print*, "LINEAR RESPONSE ONLY WORKS WITH CLOSED SHELL"
     else
       call linear_response(morb_coefat,morb_energy)
     endif
   endif

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
! TODO: Why is TD being called again?? In any case have TD call SCF as a first
!       step to obtain the first density...
!
      if (timedep.eq.1) then
        call g2g_timer_sum_start('TD')
        if (OPEN) then
           call TD(fock_aop, rho_aop, fock_bop, rho_bop)
        else
           call TD(fock_aop, rho_aop)
        end if
        call g2g_timer_sum_stop('TD')
      endif


!------------------------------------------------------------------------------!
! TODO: Cublas should be handled differently. Hidden behind SOP type and an
!       interface or general endstep call that takes care of cublas_shutdown
!       and all other similar stuff.
!
      call cublas_release( dev_Xmat )
      call cublas_release( dev_Ymat )
      call cublas_release( )

      call g2g_timer_stop('SCF')
      call g2g_timer_sum_stop('Finalize SCF')
      call g2g_timer_sum_stop('SCF')
      call g2g_timer_stop('SCF_full')
      end subroutine SCF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
