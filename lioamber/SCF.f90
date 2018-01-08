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
   use ehrensubs, only: ehrensetup
   use garcha_mod, only : M,Md, NCO,natom,Nang, number_restr, hybrid_converg, MEMO, &
   npas, verbose, RMM, X, SHFT, GRAD, npasw, igrid, energy_freq, converge,          &
   noconverge, cubegen_only, cube_dens, cube_orb, cube_elec, VCINP, Nunp, GOLD,     &
   igrid2, predcoef, nsol, r, pc, timedep, tdrestart, DIIS, told, Etold, Enucl,     &
   Eorbs, kkind,kkinds,cool,cools,NMAX,Dbug, idip, Iz, epsilon, nuc,                &
   doing_ehrenfest, first_step, RealRho, tdstep, total_time, field, Fx, Fy, Fz, a0, &
   MO_coef_at, Smat, lowdin, good_cut, ndiis
   use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC, &
   FOCK_ECP_read,FOCK_ECP_write,IzECP
   use transport, only : generate_rho0
   use time_dependent, only : TD
   use faint_cpu77, only: int1, int2, intsol, int3mem, int3lu, intfld
   use dftb_data, only : dftb_calc, MDFTB, MTB
   use dftb_subs, only : dftb_init, getXY_DFTB, find_neighbors, build_chimera,      &
                            extract_rhoDFT

   use typedef_sop   , only: sop              ! Testing SOP
   use atompot_subs  , only: atompot_oldinit  ! Testing SOP
   use tmpaux_SCF    , only: neighbor_list_2e, starting_guess
   use liosubs_dens  , only: builds_densmat, messup_densmat
   use linear_algebra, only: matrix_diagon
   use converger_subs, only: converger_init, conver


   implicit none
   integer :: Nel
   integer :: niter
   real*8  :: sq2
   real*8  :: good
   real*8  :: del
   real*8  :: factor !estan en una parte del codigo q no se usa. consultar para sacarlas
   integer :: IDAMP !estan en una parte del codigo q no se usa. consultar para sacarlas
   real*8  :: DAMP0
   real*8  :: DAMP
   integer :: igpu

   real*8, allocatable :: rho_test(:,:)
   real*8, allocatable :: rho_0(:,:)
   real*8, allocatable :: fock_0(:,:)
   real*8, allocatable :: rho(:,:)
   real*8, allocatable :: fock(:,:)
   real*8, allocatable :: morb_coefon(:,:)
   real*8, allocatable :: morb_coefat(:,:)
   real*8, allocatable :: morb_energy(:)
   integer             :: i0, ii, jj, kk, kkk

!------------------------------------------------------------------------------!
!carlos: variable mas comoda para inputs
   integer :: M_in
   integer :: NCO_in


!------------------------------------------------------------------------------!
! FFR - vvterm
   logical             :: dovv
   real*8              :: weight
   real*8, allocatable :: fockbias(:,:)
   real*8, allocatable :: Xmat(:,:)
   real*8, allocatable :: Ymat(:,:)
   real*8, allocatable :: Dvec(:)
   type(sop)           :: overop      ! Testing SOP
   real*8, allocatable :: sqsmat(:,:) ! Testing SOP

! FFR - ehrenfest (temp)
   real*8 :: dipxyz(3)
   real*8 :: dipole_norm

!FIELD variables (maybe temporary)
   real*8 :: Qc, Qc2, g


!------------------------------------------------------------------------------!
! Energy contributions and convergence
   real*8, intent(inout) :: E ! Total SCF energy 

   real*8 :: E1          ! kinetic + nuclear attraction + e-/MM charge
                         !    interaction + effective core potetial
   real*8 :: E1s         ! kinetic + nuclear attraction + effective core
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
#ifdef  CUBLAS
   integer sizeof_real
   parameter(sizeof_real=8)
   integer stat
   integer*8 devPtrX, devPtrY
   external CUBLAS_INIT, CUBLAS_SET_MATRIX
   external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_FREE
   integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX 
#endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! TODO : Variables to eliminate...
   real*8, allocatable :: xnano(:,:)
   integer :: MM, MM2, MMd, Md2
   integer :: M1, M2, M3, M5, M7, M9, M11, M13, M15, M17, M18, M19, M20

   real*8, allocatable :: Y(:,:)
   real*8, allocatable :: Ytrans(:,:)
   real*8, allocatable :: Xtrans(:,:)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


   call g2g_timer_start('SCF_full')



!--------------------------------------------------------------------!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!carlos: init para TB

   allocate (fock_0(M,M), rho_0(M,M))

   if (dftb_calc) then
      call dftb_init(M)
      allocate(fock(MDFTB,MDFTB), rho(MDFTB,MDFTB))
      allocate(morb_energy(MDFTB), morb_coefat(MDFTB,MDFTB))
   else
      allocate(fock(M,M), rho(M,M))
      allocate(morb_energy(M), morb_coefat(M,M))
   end if

   M_in = M
   if (dftb_calc) M_in=MDFTB

   if ( allocated(morb_coefon) ) deallocate(morb_coefon)
   if ( allocated(morb_coefat) ) deallocate(morb_coefat)
   if ( allocated(morb_energy) ) deallocate(morb_energy)

   allocate( morb_coefon(M_in,M_in), morb_coefat(M_in,M_in) )
   allocate( morb_energy(M_in) )


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%    Effective Core Potential Fock    %%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! TODO: unnecessary detail here in SCF, maybe replace by a subroutine
!       ecp_init() with all of these insides? Or even better, have
!       ecp module realize if it has to initialize stuff or not and
!       do it inside the sub that adds to fock.
        if (ecpmode) then
         if (FOCK_ECP_read) then
            call intECP(0) ! alocatea variables comunes y las lee del archivo ECP_restart
         else
            call g2g_timer_start('ECP Routines')
            call intECP(1) !alocatea variables, calcula variables comunes, y calcula terminos de 1 centro
            call intECP(2) !calcula terminos de 2 centros
            call intECP(3) !calcula terminos de 3 centros
              call g2g_timer_stop('ECP Routines')
         end if
         if (FOCK_ECP_write) call WRITE_ECP()
         call WRITE_POST(1)
       end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      call g2g_timer_start('SCF')
      call g2g_timer_sum_start('SCF')
      call g2g_timer_sum_start('Initialize SCF')

      npas=npas+1
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


!------------------------------------------------------------------------------!
! TODO: RMM should no longer exist. As a first step, maybe we should
!       put all these pointers inside the maskrmm module instead of
!       having them around all subroutines...
      sq2=sqrt(2.D0)
      MM=M*(M+1)/2
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M

      M1=1 ! first P
      M3=M1+MM ! now Pnew
      M5=M3+MM! now S, F also uses the same position after S was used
      M7=M5+MM! now G
      M9=M7+MMd ! now Gm
      M11=M9+MMd! now H
      M13=M11+MM! W ( eigenvalues ), also this space is used in least squares
      M15=M13+M! aux ( vector for ESSl)
      M17=M15+MM! Least squares
      M18=M17+MMd! vectors of MO
      M19=M18+M*NCO! weights (in case of using option )
      M20 = M19 + natom*50*Nang ! RAM storage of two-electron integrals (if MEMO=T)


!------------------------------------------------------------------------------!
! TODO: This seem like it should go in a separated subroutine. Also, please lets
!       not open timers outside conditionals and close them inside of them
!       (I know it is because of the "return" clause; even so this is not clear
!       programming)
       if (cubegen_only.and.(cube_dens.or.cube_orb.or.cube_elec)) then
          if (.not.VCINP) then
             write(*,*) "cubegen_only CAN ONLY BE USED WITH VCINP"
             stop
          endif

          kkk=0
          do kk=1,NCO
          do ii=1,M
             kkk=kkk+1
             morb_coefat(kk,ii) = MO_coef_at(kkk)
          enddo
          enddo

          call g2g_timer_sum_start('cube gen')
          call cubegen(M15,morb_coefat)
          call g2g_timer_sum_stop('cube gen')

          call g2g_timer_sum_stop('Initialize SCF')
          call g2g_timer_sum_stop('SCF')
          return
       endif


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
! FFR - vvterm : Variable Allocation
! TODO: this is not the proper way to do this, SCF should not have to handle
!       the bias, it should be added directly into fock via a subroutine.
       allocate(fockbias(M,M))
       dovv=.false.


!------------------------------------------------------------------------------!
! TODO: this whole part which calculates the non-electron depending terms of
!       fock and the overlap matrix should probably be in a separated sub.
!       (diagonalization of overlap, starting guess, the call to TD, should be taken out)
!
! Reformat from here...

! Para hacer lineal la integral de 2 electrones con lista de vecinos. Nano
!
      call neighbor_list_2e()

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

      if (predcoef.and.npas.gt.3) then
        write(*,*) 'no devería estar aca!'
      endif

! Calculate 1e part of F here (kinetic/nuc in int1, MM point charges
! in intsol)
!
      call g2g_timer_sum_start('1-e Fock')
      call g2g_timer_sum_start('Nuclear attraction')
      call int1(En)

! Effective Core Potential Add
!
! TODO: Maybe condense this in a single subroutine that receives only
!       RMM and returns it modified as needed?
!
      if (ecpmode ) then
          write(*,*) "Modifying Fock Matrix with ECP terms"
          do kk=1,MM
               term1e(kk)=RMM(M11+kk-1) !backup of 1e terms
               RMM(M11+kk-1)=RMM(M11+kk-1)+VAAA(kk)+VAAB(kk)+VBAC(kk) !add EC
          enddo
      end if

! Other terms
!
      call g2g_timer_sum_stop('Nuclear attraction')
      if(nsol.gt.0.or.igpu.ge.4) then
          call g2g_timer_sum_start('QM/MM')
       if (igpu.le.1) then
          call g2g_timer_start('intsol')
          call intsol(E1s,Ens,.true.)
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
! TODO: remove or sistematize
!
      E1=0.D0
      do kk=1,MM
        E1=E1+RMM(kk)*RMM(M11+kk-1)
      enddo
      call g2g_timer_sum_stop('1-e Fock')


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! OVERLAP DIAGONALIZATION
! TODO: Simplify, this has too much stuff going on...
! (maybe trans mats are not even necessary?)
!
        if ( allocated(Xmat) ) deallocate(Xmat)
        if ( allocated(Ymat) ) deallocate(Ymat)
        allocate( Xmat(M,M), Ymat(M,M) )

        if ( allocated(Y) )      deallocate(Y)
        if ( allocated(Ytrans) ) deallocate(Ytrans)
        if ( allocated(Xtrans) ) deallocate(Xtrans)
        allocate( Y(M,M), Ytrans(M,M), Xtrans(M,M) )

        call overop%Sets_smat( Smat )
        if (lowdin) then
!          TODO: LOWDIN CURRENTLY NOT WORKING (NOR CANONICAL)
           call overop%Gets_orthog_4m( 2, 0.0d0, Xmat, Ymat, Xtrans, Ytrans )
        else
           call overop%Gets_orthog_4m( 1, 0.0d0, Xmat, Ymat, Xtrans, Ytrans )
        end if

! TODO: initializations related to atompot should be dealt with differently
        if (dovv) then
           if ( allocated(fockbias) ) deallocate(fockbias)
           if ( allocated(sqsmat) )   deallocate(sqsmat)
           allocate( fockbias(M,M), sqsmat(M,M) )
           call atompot_oldinit( natom, nuc, sqsmat, fockbias )
        end if

! TODO: replace X,Y,Xtrans,Ytrans with Xmat, Ymat, Xtrp, Ytrp
        do ii=1,M
        do jj=1,M
           X(ii,jj)      = Xmat(ii,jj)
           Y(ii,jj)      = Ymat(ii,jj)
!          Xtrans(ii,jj) = Xtrp(ii,jj)
!          Ytrans(ii,jj) = Ytrp(ii,jj)
        end do
        end do

        if ( allocated(Dvec) ) deallocate(Dvec)
        allocate( Dvec(M) )
        call overop%Gets_eigens_v( 0.0d0, Dvec )
        do kk = 1, M
           RMM(M13+kk-1) = Dvec(kk)
        end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!carlosDFTB: allocation of X and Y matrix for TB calculations
!
! TODO: this is nasty, a temporary solution would be to have a Msize variable
!       be assigned M (or, even better, "basis_size") or MDFTB
!       ("basis_size_dftb") according to the case

       if (dftb_calc) then

          M_in=MDFTB

          if ( allocated(Xmat) ) deallocate(Xmat)
          allocate( Xmat(MDFTB,MDFTB) )

          if ( allocated(Ymat) ) deallocate(Ymat)
          allocate( Ymat(MDFTB,MDFTB) )
       
          call getXY_DFTB(M,x,y,xmat,ymat)
       
       else

          M_in=M
          
          if ( allocated(Xmat) ) deallocate(Xmat)
          allocate( Xmat(M,M) )

          if ( allocated(Ymat) ) deallocate(Ymat)
          allocate( Ymat(M,M) )

          do jj = 1, M
          do ii = 1, M
             Xmat(ii,jj) = x(ii,jj)
             Ymat(ii,jj) = y(ii,jj)
          enddo
          enddo
       
      end if 


! TODO: are these comments relevant?
!
!      MDFTB must be declared before
!      call :
!      agrandar y modificar las Xmat/Ymat
!      modificar MDFTB

! CUBLAS
! TODO: these should be inside of the overlap type
!
#ifdef CUBLAS
            call CUBLAS_INIT()
            stat = CUBLAS_ALLOC(M_in*M_in, sizeof_real, devPtrX)
            stat = CUBLAS_ALLOC(M_in*M_in, sizeof_real, devPtrY)
            if (stat.NE.0) then
              write(*,*) "X and/or Y memory allocation failed"
              call CUBLAS_SHUTDOWN
              stop
            endif
            stat = CUBLAS_SET_MATRIX(M_in,M_in,sizeof_real,Xmat,M_in,devPtrX,M_in)
            stat = CUBLAS_SET_MATRIX(M_in,M_in,sizeof_real,Ymat,M_in,devPtrY,M_in)
            if (stat.NE.0) then
              write(*,*) "X and/or Y setting failed"
              call CUBLAS_SHUTDOWN
              stop
            endif 
#endif


!##############################################################################!
! FFR: Currently working here.
!
      write(666,*) "X on input...", size(X,1), size(X,2)
      do ii=1,size(X,1)
      do jj=1,size(X,2)
         if ((ii>M).or.(jj>M)) then
            write(666,*) ii , jj , X(ii,jj)
         endif
      enddo
      enddo
      write(666,*)
      call starting_guess( morb_coefat )
      write(666,*) "morb_coefat on output..."
      do ii=1,M
      do jj=1,M
         write(666,*) ii , jj , morb_coefat(ii,jj)
      enddo
      enddo
!      stop
!
! FFR: When finished, uncomment the following starting guess...
!##############################################################################!
!------------------------------------------------------------------------------!
! TODO: remove from here...
!
!      call starting_guess( morb_coefat )

      if ((timedep.eq.1).and.(tdrestart)) then
        call g2g_timer_sum_start('TD')
        call TD()
        call g2g_timer_sum_stop('TD')
        return
      endif


!----------------------------------------------------------!
! Precalculate two-index (density basis) "G" matrix used in density fitting
! here (S_ij in Dunlap, et al JCP 71(8) 1979) into RMM(M7)
! Also, pre-calculate G^-1 if G is not ill-conditioned into RMM(M9)
!
      call g2g_timer_sum_start('Coulomb G matrix')
      call int2()
      call g2g_timer_sum_stop('Coulomb G matrix')
!
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
         call int3mem()
!        Small elements of t_i put into single-precision cools here
!        call int3mems()
         call g2g_timer_stop('int3mem')
         call g2g_timer_sum_stop('Coulomb precalc')
      endif
!
!
! TODO: ...to here
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

        if (verbose) call WRITE_CONV_STATUS(GOOD,TOLD,EGOOD,ETOLD)
!       Escribe los criterios de convergencia y el valor del paso de dinamica

        call g2g_timer_start('Total iter')
        call g2g_timer_sum_start('Iteration')
        call g2g_timer_sum_start('Fock integrals')
        niter=niter+1
        E1=0.0D0


!------------------------------------------------------------------------------!
!       Fit density basis to current MO coeff and calculate Coulomb F elements
!
! TODO: Calculation of fock terms should be separated. Also, maybe it would be
!       convenient to add the NaN checks and the timers inside the subroutine
!       calls is a more systematic way.

!       Test for NaN
        if (Dbug) call SEEK_NaN(RMM,1,MM,"RHO Start")
        if (Dbug) call SEEK_NaN(RMM,M5-1,M5-1+MM,"FOCK Start")

!       Computes Coulomb part of Fock, and energy on E2
        call g2g_timer_sum_start('Coulomb fit + Fock')
        call int3lu(E2)
        call g2g_timer_sum_pause('Coulomb fit + Fock')

!       Test for NaN
        if (Dbug) call SEEK_NaN(RMM,1,MM,"RHO Coulomb")
        if (Dbug) call SEEK_NaN(RMM,M5-1,M5-1+MM,"FOCK Coulomb")

!       XC integration / Fock elements
        call g2g_timer_sum_start('Exchange-correlation Fock')
        call g2g_solve_groups(0,Ex,0)
        call g2g_timer_sum_pause('Exchange-correlation Fock')

!       Test for NaN
        if (Dbug) call SEEK_NaN(RMM,1,MM,"RHO Ex-Corr")
        if (Dbug) call SEEK_NaN(RMM,M5-1,M5-1+MM,"FOCK Ex-Corr")


!------------------------------------------------------------------------------!
! REACTION FIELD CASE
!
! TODO: what is reaction field?? in any case, it should be either in its own
!       module and only appear here as a subroutine, or be included in the
!       separated fock calculation subroutine...
!
        if ( field .and. generate_rho0 ) then
           dipxyz(:)=0.0D0
           call dip(dipxyz)
           g=1.0D0
           factor=2.54D0
           call intfld(g,Fx,Fy,Fz)
           E1=-1.00D0*g*(Fx*dipxyz(1)+Fy*dipxyz(2)+Fz*dipxyz(3))/factor- &
           0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0

           do kk=1,MM
               E1=E1+RMM(kk)*RMM(M11+kk-1)
           enddo

        else
!          E1 includes solvent 1 electron contributions
           do kk=1,MM
              E1=E1+RMM(kk)*RMM(M11+kk-1)
           enddo

        endif
        call g2g_timer_start('actualiza rmm')
        call g2g_timer_sum_pause('Fock integrals')


!------------------------------------------------------------------------------!
!carlos: extraemos rho y fock antes
!
! TODO: extraction of fock an rho via subroutines from maskrmm as a first step,
!       total removal once rmm is gone.

        do jj=1,M

          do kk=1,jj-1
            fock_0(jj,kk) = RMM(M5+jj+(M2-kk)*(kk-1)/2-1)
            rho_0(jj,kk)  = (RMM(jj+(M2-kk)*(kk-1)/2))/2
          enddo

          fock_0(jj,jj) = RMM(M5+jj+(M2-jj)*(jj-1)/2-1)
          rho_0(jj,jj)  = RMM(jj+(M2-jj)*(jj-1)/2)

          do kk=jj+1,M
            fock_0(jj,kk) = RMM(M5+kk+(M2-jj)*(jj-1)/2-1)
            rho_0(jj,kk)  = RMM(kk+(M2-jj)*(jj-1)/2)/2
          enddo

        enddo


!------------------------------------------------------------------------------!
!carlos: armamos la fock posta
!
! TODO: this should be wrapped inside a single dftb subroutine. Also, two
!       consecutive dftb_calc switches? really?
!
        if (dftb_calc) then
          call find_neighbors(M_in,Nuc,natom)
          call build_chimera (M, fock_0, fock, natom)

          if (niter==1) then
            rho=0.0D0
            do ii=1, MTB
              rho(ii,ii)=1.0D0
              rho(MTB+M+ii,MTB+M+ii)=1.0D0 
            end do
          end if
    
          rho(MTB+1:MTB+M, MTB+1:MTB+M)=rho_0(:,:)

        else
          fock = fock_0
          rho  = rho_0
        endif

        if (dftb_calc) then
          NCO_in = NCO+MTB
        else
          NCO_in = NCO
        end if


!------------------------------------------------------------------------------!
        call g2g_timer_sum_start('SCF acceleration')
        if (niter==1) then
           call converger_init( M_in, ndiis, DAMP, DIIS, hybrid_converg )
        end if
#       ifdef CUBLAS
           call conver(niter, good, good_cut, M_in, rho, fock, devPtrX, devPtrY)
#       else
           call conver(niter, good, good_cut, M_in, rho, fock, Xmat, Ymat)
#       endif
        call g2g_timer_sum_pause('SCF acceleration')


!------------------------------------------------------------------------------!
!  Fock(ON) diagonalization
        if ( allocated(morb_coefon) ) deallocate(morb_coefon)
        allocate( morb_coefon(M_in,M_in) )

        call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
        call matrix_diagon( fock, morb_coefon, morb_energy )
        call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')


!------------------------------------------------------------------------------!
!  Base change of coeficients ( (X^-1)*C ) and construction of new density
!  matrix
        call g2g_timer_sum_start('SCF - MOC base change (sum)')
#       ifdef CUBLAS
           call cumxp_r( morb_coefon, devPtrX, morb_coefat, M_in)
#       else
           morb_coefat = matmul( Xmat, morb_coefon )
#       endif
        call g2g_timer_sum_pause('SCF - MOC base change (sum)')

        if ( allocated(morb_coefon) ) deallocate(morb_coefon)

        call builds_densmat( M_in, NCO_in, 2.0d0, morb_coefat, rho)
        call messup_densmat( rho )


!------------------------------------------------------------------------------!
!  We are not sure how to translate the sumation over molecular orbitals
!  and energies when changing from DFTB system to DFT subsystem. Forces
!  may be broken due to this.
!  This should not be affecting normal DFT calculations.
        i0 = 0
        if (dftb_calc) i0=MTB

        do kk=1,M
          RMM(M13+kk-1) = morb_energy(kk)
        end do

        kkk = 0
        do kk=1,NCO
        do ii=1,M
          kkk = kkk+1
          MO_coef_at(kkk) = morb_coefat( i0+ii, kk )
        enddo
        enddo

        do ii=1,M
!         do jj=1,NCO
          do jj=1,M
            X( ii, M2+jj ) = morb_coefat( i0+ii, jj )
          enddo
!         do jj=NCO,M
!           X( ii, M2+jj ) = morb_coefat( i0+ii, jj+2*MTB )
!         enddo
        enddo


!------------------------------------------------------------------------------!
!carlos: agregado para separar de rho la parte DFT
!
! TODO: again, this should be handled differently...
! TODO: make xnano go away...only remains here
!
        allocate ( xnano(M,M) )

        if (dftb_calc) then
          call extract_rhoDFT(M, rho, xnano)
        else
          xnano=rho
        end if


!------------------------------------------------------------------------------!
! TODO: convergence criteria should be a separated subroutine...
        good = 0.0d0
        do jj=1,M
        do kk=jj,M
          del=xnano(jj,kk)-(RMM(kk+(M2-jj)*(jj-1)/2))
          del=del*sq2
          good=good+del**2
          RMM(kk+(M2-jj)*(jj-1)/2)=xnano(jj,kk)
        enddo
        enddo
        good=sqrt(good)/float(M)
        deallocate ( xnano )


! TODO: what is this doing here???
        call g2g_timer_stop('dens_GPU')

!------------------------------------------------------------------------------!
! TODO: finalization of the loop is a little bit messy. Also: "999 continue"??
!       I think it is time we regularized this loop...

! Damping factor update
        DAMP=DAMP0

        E=E1+E2+En
!        E=E+Es
!
        call g2g_timer_stop('otras cosas')
        call g2g_timer_sum_pause('new density')

        if(verbose) call WRITE_E_STEP(niter, E+Ex) !escribe energia en cada paso

        Egood=abs(E+Ex-Evieja)
        Evieja=E+Ex
!
        call g2g_timer_stop('Total iter')
        call g2g_timer_sum_pause('Iteration')

 999  continue



!------------------------------------------------------------------------------!
      call g2g_timer_sum_start('Finalize SCF')

      if (niter.ge.NMAX) then
        write(6,*) 'NO CONVERGENCE AT ',NMAX,' ITERATIONS'
        noconverge=noconverge + 1
        converge=0
      else
        write(6,*) 'CONVERGED AT',niter,'ITERATIONS'
        noconverge = 0
        converge=converge+1
      endif

      if(noconverge.gt.4) then
        write(6,*)  'stop for not convergion 4 times'
        stop
      endif

!------------------------------------------------------------------------------!
! TODO: Comments about a comented sections? Shouldn't it all of this go away?
!
!    CH - Why call intsol again here? with the .false. parameter,
!    E1s is not recalculated, which would be the only reason to do
!    this again; Ens isn't changed from before...
! -- SOLVENT CASE --------------------------------------
!      if (sol) then
!      call g2g_timer_sum_start('intsol 2')
!      if(nsol.gt.0) then
!        call intsol(E1s,Ens,.false.)
!        write(*,*) 'cosillas',E1s,Ens
!        call g2g_timer_sum_stop('intsol 2')
!      endif
!      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
      Es=Es+E1s+Ens
!     endif


!------------------------------------------------------------------------------!
      if (MOD(npas,energy_freq).eq.0) then
      if (GRAD) then
!       call g2g_timer_sum_start('exchnum')
        call g2g_timer_sum_start('Exchange-correlation energy') 
!       Resolve with last density to get XC energy
        call g2g_new_grid(igrid)
        call g2g_solve_groups(1, Exc, 0)
        call g2g_timer_sum_stop('Exchange-correlation energy')

!----------- COmputing the QM/MM contribution to total energy

!       NucleusQM-CHarges MM
        Es=Ens
        
!       One electron Kinetic (with aint >3) or Kinetic + Nuc-elec (aint >=3)
        call int1(En)

!       Computing the E1-fock without the MM atoms
        if (nsol.gt.0.and.igpu.ge.1) then
          call aint_qmmm_init(0,r,pc)
          call aint_qmmm_fock(E1s,Etrash)
          call aint_qmmm_init(nsol,r,pc)
        endif

!       E1s (here) is the 1e-energy without the MM contribution
        E1s=0.D0
        do kk=1,MM
          E1s=E1s+RMM(kk)*RMM(M11+kk-1)
        enddo

!       Es is the QM/MM energy computated as total 1e - E1s + QMnuc-MMcharges
        Es=Es+E1-E1s
 
! -------------------------------------------------
!       Total SCF energy =
!       E1 - kinetic + nuclear attraction + QM/MM interaction + effective core potential
!       E2 - Coulomb
!       En - nuclear-nuclear repulsion
!       Ens - MM point charge - nuclear interaction
!       Exc - exchange-correlation
!       Eecp - Efective core potential
!       E_restrain - distance restrain
! -------------------------------------------------

!       Part of the QM/MM contrubution are in E1
        E=E1+E2+En+Ens+Exc+E_restrain


!%%%%%%%%%%%%%%   Write Energy Contributions   %%%%%%%%%%%%%%
        if (npas.eq.1) npasw = 0

        if (npas.gt.npasw) then  
          if (ecpmode) then
            Eecp=0.d0
            do kk=1,MM
              Eecp=Eecp+RMM(kk)*(VAAA(kk)+VAAB(kk)+VBAC(kk))
            enddo
            Es=Es-Eecp  ! 
          end if
          call WriteEnergies(E1,E2,En,Ens,Eecp,Exc,ecpmode,E_restrain)
          npasw=npas+10
        endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      endif ! GRAD
      endif ! npas




!------------------------------------------------------------------------------!
! calculation of energy weighted density matrix
!
! TODO: this should be a separated subroutine...
!
      call g2g_timer_sum_start('energy-weighted density')
      kkk=0
      do jj=1,M
      do ii=jj,M
         kkk=kkk+1
         RMM(M15+kkk-1)=0.D0
!
         if (ii.eq.jj) then
            factor=2.D0
         else
            factor=4.D0
         endif

         do kk=1,NCO
            RMM(M15+kkk-1)= &
            RMM(M15+kkk-1)-RMM(M13+kk-1)*factor*X(ii,M2+kk)*X(jj,M2+kk)
         enddo
      enddo
      enddo

      call g2g_timer_sum_stop('energy-weighted density')

!     Variables needed for further calculations (Populations, Dip, etc).
      Enucl = En
      do kkk=1, M
          Eorbs(kkk) = RMM(M13+kkk-1)
      enddo


!------------------------------------------------------------------------------!
! FFR - Ehrenfest

! TODO: have ehrendyn call SCF and have SCF always save the resulting rho in
!       a module so that ehrendyn can retrieve it afterwards. Remove all this.
      if (doing_ehrenfest) then
         call spunpack('L',M,RMM(M1),RealRho)
         call fixrho(M,RealRho)
         call ehrensetup(natom, M, RealRho)
      endif

! TODO: have a separate module handle the dipole moment
       if (idip .eq. 1) then
           if (first_step) then
             call write_dipole(dipxyz, 0, 134, .true.)
             total_time=0.0d0
           endif
           dipxyz(:)=0.0D0
           call dip(dipxyz)
           dipole_norm = sqrt(dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)
           call write_dipole(dipxyz, dipole_norm, 134, .false.)
       
           print*,''
           print*,' Timer: ',total_time
           print*,''
           total_time=total_time+tdstep*0.0241888
       endif

 901  format(F15.9,2x,F15.9)


!------------------------------------------------------------------------------!
! Performs orbital/density plots.
!
! TODO: this should be a separated subroutine...
!
        if (cube_dens.or.cube_orb.or.cube_elec) then
          call g2g_timer_sum_start('cube gen')
          kk=0

          do kk=1,M
          do ii=1,M
             kkk=kkk+1 ! TODO: unnecessary?
             morb_coefat(kk,ii)  = X(ii,M2+kk)
          enddo
          enddo

          call cubegen(M15,morb_coefat)
          call g2g_timer_sum_stop('cube gen')
        endif


!------------------------------------------------------------------------------!
! TODO: Deallocation of variables that should be removed
! TODO: MEMO should be handled differently...
!
      deallocate (Y)
      deallocate (Ytrans, Xtrans )


      if (MEMO) then
        deallocate (kkind,kkinds)
        deallocate(cool,cools)
      endif

!------------------------------------------------------------------------------!
! TODO: Why is TD being called again?? remove from here
!
      if(timedep.eq.1) then
        call g2g_timer_sum_start('TD')
        call TD()
        call g2g_timer_sum_stop('TD')
      endif


!------------------------------------------------------------------------------!
! TODO: Cublas should be handled differently. Hidden behind SOP type and an
!       interface or general endstep call that takes care of cublas_shutdown
!       and all other similar stuff.
!
#ifdef  CUBLAS
      call CUBLAS_FREE(devPtrX)
      call CUBLAS_FREE(devPtrY)
      call CUBLAS_SHUTDOWN 
#endif


!--------------------------------------------------------------------!
      call g2g_timer_stop('SCF')
      call g2g_timer_sum_stop('Finalize SCF')
      call g2g_timer_sum_stop('SCF')
      call g2g_timer_stop('SCF_full')
      end subroutine SCF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
