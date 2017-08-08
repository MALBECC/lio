!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! SCF subroutine
! DIRECT VERSION
! Calls all integrals generator subroutines : 1 el integrals,
! 2 el integrals, exchange fitting , so it gets S matrix, F matrix
! and P matrix in lower storage mode ( symmetric matrices)
!
! Modified to f90, 2017
! Dario Estrin, 1992
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      subroutine SCF(E)
!      use linear_algebra
      use ehrenfest, only: ehrensetup
      use garcha_mod, only : M,Md, NCO,natom,Nang, number_restr, hybrid_converg, MEMO, &
      npas, verbose, RMM, X, SHFT, GRAD, npasw, igrid, energy_freq, converge,          &
      noconverge, cubegen_only, cube_dens, cube_orb, cube_elec, VCINP, Nunp, GOLD,     &
      igrid2, predcoef, nsol, r, pc, timedep, tdrestart, DIIS, told, Etold, Enucl,     &
      Eorbs, kkind,kkinds,cool,cools,NMAX,Dbug, idip                                   &
      , do_ehrenfest, first_step, RealRho, tdstep, total_time 
!      use mathsubs
      use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC, &
       FOCK_ECP_read,FOCK_ECP_write,IzECP
!      use general_module 
!#ifdef  CUBLAS
!      use cublasmath 
!#endif
	IMPLICIT NONE
	integer :: M1,M2,M3, M5, M7, M9, M11, M13, M15, M17, M18,M19, M20, MM,MM2,MMd,  &
        Md2 !temporales hasta q rompamos RMM
	integer :: i,j,k, kk !Auxiliares
	integer :: Nel, niter
	real*8 :: sq2, ff
	REAL*8 :: good, del
	REAL*8 :: D1, D2, E0, factor !estan en una parte del codigo q no se usa. consultar para sacarlas
	INTEGER :: IDAMP !estan en una parte del codigo q no se usa. consultar para sacarlas
	REAL*8 :: DAMP0, DAMP
	INTEGER :: igpu
      integer:: l
       real*8, dimension (:,:), ALLOCATABLE::xnano,znano
!       real*8, dimension (:,:), ALLOCATABLE::scratch1, scratch
       real*8, dimension (:), ALLOCATABLE :: rmm5,rmm15,rmm13
      real*8, dimension (:,:), allocatable :: Y,Ytrans,Xtrans
       logical  hagodiis, ematalloc
        INTEGER :: ErrID,iii,jjj
        LOGICAL :: docholesky
        INTEGER            :: LWORK2
        INTEGER, ALLOCATABLE :: IPIV(:)
        logical :: just_int3n,ematalloct

! FFR - vvterm
!----------------------------------------------------------!
       logical             :: dovv
       real*8              :: weight, softness
       real*8,dimension(:,:),allocatable :: fockbias

! FFR - ehrenfest (temp)
!----------------------------------------------------------!
       real*8 :: dipxyz(3), dipole_norm

!----------------------------------------------------------!


! Energy and contributions
       real*8 :: E, E1, E1s, E2,Eecp, En, Ens, Es, E_restrain, Ex, Exc,Etrash
!--------------------------------------------------------------------!

 
#ifdef  CUBLAS
        integer sizeof_real
        parameter(sizeof_real=8)
        integer stat
        integer*8 devPtrX, devPtrY
        external CUBLAS_INIT, CUBLAS_SET_MATRIX
        external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_FREE
        integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX 
#endif

!------------------------------------------------------

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        double precision :: Egood, Evieja !variables para criterio de co
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!------------------------------------------------------

      call g2g_timer_start('SCF_full')
!--------------------------------------------------------------------!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%    Effective Core Potential Fock    %%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
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
      ematalloc=.false.
      hagodiis=.false.

      allocate (znano(M,M),xnano(M,M))

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

      if (cubegen_only.and.(cube_dens.or.cube_orb.or.cube_elec)) then
        if (.not.VCINP) then
          write(*,*) "cubegen_only CAN ONLY BE USED WITH VCINP"
          stop
        endif
        kk=0
        do k=1,NCO
          do i=1,M
            kk=kk+1
            Xnano(k,i) = RMM(M18+kk-1)
          enddo
        enddo

        call g2g_timer_sum_start('cube gen')
        call cubegen(M15,Xnano)
        call g2g_timer_sum_stop('cube gen')

        deallocate (znano,xnano)
!,scratch,scratch1)
        call g2g_timer_sum_stop('Initialize SCF')
        call g2g_timer_sum_stop('SCF')
        return
      endif
!
      Nel=2*NCO+Nunp
!
      allocate(rmm5(MM),rmm15(mm))
!
      good=1.00D0
      Egood=1.00D0
      Evieja=0.d0

      niter=0
      D1=1.D0
      D2=1.D0
      DAMP0=GOLD
      DAMP=DAMP0

! FFR - vvterm : Variable Allocation
!--------------------------------------------------------------------!
       allocate(fockbias(M,M))
       dovv=.false.
!----------------------------------------
      call neighbor_list_2e() ! Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano

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

!
! Calculate 1e part of F here (kinetic/nuc in int1, MM point charges
! in intsol)
!
      call g2g_timer_sum_start('1-e Fock')
      call g2g_timer_sum_start('Nuclear attraction')
      call int1(En)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%    Effective Core Potential Add    %%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      if (ecpmode ) then
          write(*,*) "Modifying Fock Matrix with ECP terms"
          do k=1,MM
               term1e(k)=RMM(M11+k-1) !backup of 1e terms
               RMM(M11+k-1)=RMM(M11+k-1)+VAAA(k)+VAAB(k)+VBAC(k) !add EC
          enddo
      end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

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

!
! test ---------------------------------------------------------
      E1=0.D0
      do k=1,MM
        E1=E1+RMM(k)*RMM(M11+k-1)
      enddo
      call g2g_timer_sum_stop('1-e Fock')

!#########################################################################################
!#########################################################################################
	allocate (Y(M,M),Ytrans(M,M),Xtrans(M,M))
! Diagonalization of S matrix, after this is not needed anymore
	call overlap_diag(fockbias, dovv,Y,Ytrans,Xtrans)
!#########################################################################################
!#########################################################################################

!! CUBLAS ---------------------------------------------------------------------!
#ifdef CUBLAS
            call CUBLAS_INIT()
            stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
            stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrY)
            if (stat.NE.0) then
              write(*,*) "X and/or Y memory allocation failed"
              call CUBLAS_SHUTDOWN
              stop
            endif
            stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,X,M,devPtrX,M)
            stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,Y,M,devPtrY,M)
            if (stat.NE.0) then
              write(*,*) "X and/or Y setting failed"
              call CUBLAS_SHUTDOWN
              stop
            endif 
#endif
!------------------------------------------------------------------------------!

!############################################################################
!############################################################################
	call starting_guess(xnano)
!############################################################################
!############################################################################

!
      if ((timedep.eq.1).and.(tdrestart)) then
        call g2g_timer_sum_start('TD')
        call TD()
        call g2g_timer_sum_stop('TD')
        return
      endif
!
! Precalculate two-index (density basis) "G" matrix used in density fitting
! here (S_ij in Dunlap, et al JCP 71(8) 1979) into RMM(M7)
! Also, pre-calculate G^-1 if G is not ill-conditioned into RMM(M9)
!
      call g2g_timer_sum_start('Coulomb G matrix')
      call int2()
      call g2g_timer_sum_stop('Coulomb G matrix')
!
!*
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
! Large elements of t_i put into double-precision cool here
! Size criteria based on size of pre-factor in Gaussian Product Theorem
! (applied to MO basis indices)
         call int3mem()
! Small elements of t_i put into single-precision cools here
!         call int3mems()
         call g2g_timer_stop('int3mem')
         call g2g_timer_sum_stop('Coulomb precalc')
      endif
!***
!---------------------------------------------------------------------
! Now, damping is performed on the density matrix
! The first 4 iterations ( it may be changed, if necessary)
! when the density is evaluated on the grid, the density
! matrix is used ( slower), after that it is calculated
! using the vectors . Since the vectors are not damped,
! only at the end of the SCF, the density matrix and the
! vectors are 'coherent'
!---------------------------------------------------------------
! LEVEL SHIFT CASE, contruction of initial vectors ------------------
!
      if (SHFT) then
        write(*,*) 'Level SHIFT is not suported'
        stop
      endif
!

      if (hybrid_converg) DIIS=.true. ! cambio para convergencia damping-diis

      call g2g_timer_sum_stop('Initialize SCF')

!-------------------------------------------------------------------
!-------------------------------------------------------------------
      do 999 while ((good.ge.told.or.Egood.ge.Etold).and.niter.le.NMAX)

      if (verbose) call WRITE_CONV_STATUS(GOOD,TOLD,EGOOD,ETOLD)
! Escribe los criterios de convergencia y el valor del paso de dinamica

        call g2g_timer_start('Total iter')
        call g2g_timer_sum_start('Iteration')
        call g2g_timer_sum_start('Fock integrals')
        niter=niter+1

! Fit density basis to current MO coeff and calculate Coulomb F elements
!

!-------------------------------------------------------------------
! Test for NaN
        if (Dbug) call SEEK_NaN(RMM,1,MM,"RHO Start")
        if (Dbug) call SEEK_NaN(RMM,M5-1,M5-1+MM,"FOCK Start")
!-------------------------------------------------------------------
            call g2g_timer_sum_start('Coulomb fit + Fock')
            call int3lu(E2) ! Computes Coulomb part of Fock, and energy on E2
            call g2g_timer_sum_pause('Coulomb fit + Fock')
!-------------------------------------------------------------------
! Test for NaN
        if (Dbug) call SEEK_NaN(RMM,1,MM,"RHO Coulomb")
        if (Dbug) call SEEK_NaN(RMM,M5-1,M5-1+MM,"FOCK Coulomb")
!-------------------------------------------------------------------
            call g2g_timer_sum_start('Exchange-correlation Fock')
            call g2g_solve_groups(0,Ex,0) ! XC integration / Fock elements
            call g2g_timer_sum_pause('Exchange-correlation Fock')
!-------------------------------------------------------------------
! Test for NaN
        if (Dbug) call SEEK_NaN(RMM,1,MM,"RHO Ex-Corr")
        if (Dbug) call SEEK_NaN(RMM,M5-1,M5-1+MM,"FOCK Ex-Corr")
!-------------------------------------------------------------------



!-------------------------------------------------------
        E1=0.0D0
!
! REACTION FIELD CASE --------------------------------------------
!
        call g2g_timer_start('actualiza rmm')
!----------------------------------------------------------------
! E1 includes solvent 1 electron contributions
        do k=1,MM
          E1=E1+RMM(k)*RMM(M11+k-1)
        enddo
        call g2g_timer_sum_pause('Fock integrals')
        call g2g_timer_sum_start('SCF acceleration')


#ifdef CUBLAS
        call obtain_new_P(niter, DAMP, dovv, fockbias, good, xnano, znano, devPtrX, devPtrY)
#else
        call obtain_new_P(niter, DAMP, dovv, fockbias, good, xnano, znano, Y)
#endif



      do j=1,M
         do k=j,M
               del=xnano(j,k)-(RMM(k+(M2-j)*(j-1)/2))
               del=del*sq2
               good=good+del**2
               RMM(k+(M2-j)*(j-1)/2)=xnano(j,k)
         enddo
      enddo
      good=sqrt(good)/float(M)

      call g2g_timer_stop('dens_GPU')

       if (SHFT) then
! Level Shifting
         do i=1,M
           do j=1,M
             X(i,j)=X(i,M2+j)
           enddo
         enddo
       endif
!
       if (SHFT) then
         DAMP=0.0D0
       endif
!
!--- Damping factor update -
       DAMP=DAMP0




!#################################################################################
!###################    PREGUNTAR    #############################################
!#################################################################################
       IDAMP=0
       if (IDAMP.EQ.1) then !este if es siempre false. sirve para algo??, Nick
         DAMP=DAMP0
         if (abs(D1).lt.1.D-5) then
           factor=dmax1(0.90D0,abs(D1/D2))
           factor=dmin1(factor,1.1D0)
           DAMP=DAMP0*factor
         endif
!
         E=E1+E2+En
         E=E+Es
!
         D2=D1
         D1=(E-E0)
!
         E0=E
         DAMP0=DAMP
       endif
!#################################################################################
!#################################################################################
!#################################################################################

!
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
!-------------------------------------------------------------------
!
!-------------------------------------------------------------------
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
!--------------------------------------------------------------

      if (MOD(npas,energy_freq).eq.0) then
      if (GRAD) then
!       call g2g_timer_sum_start('exchnum')
        call g2g_timer_sum_start('Exchange-correlation energy') 
#ifdef  G2G 
#ifdef  ULTIMA_CPU
        call exchnum(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM, &
                    M18,NCO,Exc,nopt)
#else
      ! Resolve with last density to get XC energy
        call g2g_new_grid(igrid)
        call g2g_solve_groups(1, Exc, 0)
#endif
#else 
#ifdef  ULTIMA_G2G
        call g2g_new_grid(igrid)
        call g2g_solve_groups(1, Exc, 0)
#else 
#endif 
#endif
        call g2g_timer_sum_stop('Exchange-correlation energy')
!----------- COmputing the QM/MM contribution to total energy

         Es=Ens   ! NucleusQM-CHarges MM

       call int1(En) ! One electron Kinetic (with aint >3) or Kinetic + Nuc-elec (aint >=3)

       if(nsol.gt.0.and.igpu.ge.1) then ! Computing the E1-fock without the MM atoms
          call aint_qmmm_init(0,r,pc)
          call aint_qmmm_fock(E1s,Etrash)
          call aint_qmmm_init(nsol,r,pc)
      endif
      E1s=0.D0
      do k=1,MM
        E1s=E1s+RMM(k)*RMM(M11+k-1)  ! E1s (here) is the 1e-energy without the MM contribution 
      enddo

      Es=Es+E1-E1s ! Es is the QM/MM energy computated as total 1e - E1s + QMnuc-MMcharges
 
        ! -------------------------------------------------
        ! Total SCF energy =
        ! E1 - kinetic+nuclear attraction+QM/MM interaction + effective core potential
        ! E2 - Coulomb
        ! En - nuclear-nuclear repulsion
        ! Ens - MM point charge-nuclear interaction
        ! Exc - exchange-correlation
        ! Eecp - Efective core potential
        ! E_restrain - distance restrain
        ! -------------------------------------------------

        E=E1+E2+En+Ens+Exc+E_restrain ! Part of the QM/MM contrubution are in E1


        if (npas.eq.1) npasw = 0

!%%%%%%%%%%%%%%   Write Energy Contributions   %%%%%%%%%%%%%%
        if (npas.gt.npasw) then  
          if (ecpmode) then
            Eecp=0.d0
            do k=1,MM
              Eecp=Eecp+RMM(k)*(VAAA(k)+VAAB(k)+VBAC(k))
            enddo
            Es=Es-Eecp  ! 
          end if
          call WriteEnergies(E1,E2,En,Ens,Eecp,Exc,ecpmode,E_restrain)
          npasw=npas+10
        endif
      endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      endif
! calculation of energy weighted density matrix
!


      call g2g_timer_sum_start('energy-weighted density')
      kk=0
      do j=1,M
        do i=j,M
          kk=kk+1
          RMM(M15+kk-1)=0.D0
!
          if(i.eq.j) then
            ff=2.D0
          else
            ff=4.D0
          endif
          do k=1,NCO
            RMM(M15+kk-1)= &
            RMM(M15+kk-1)-RMM(M13+k-1)*ff*X(i,M2+k)*X(j,M2+k)
          enddo
        enddo
      enddo

      call g2g_timer_sum_stop('energy-weighted density')

!     Variables needed for further calculations (Populations, Dip, etc).
      Enucl = En
      do kk=1, M
          Eorbs(kk) = RMM(M13+kk-1)
      enddo

!==============================================================================!
! FFR - Ehrenfest

! TODO: have ehrendyn call SCF and have SCF always save the resulting rho in
!       a module so that ehrendyn can retrieve it afterwards. Remove all this.
      if (do_ehrenfest) then
         call spunpack('L',M,RMM(M1),RealRho)
         call fixrho(M,RealRho)
         call ehrensetup(M,RealRho)
      endif

! TODO: have a separate module handle the dipole moment
       if (idip .eq. 1) then
           if (first_step) then
             call write_dipole(dipxyz, 0, 134, .true.)
             total_time=0.0d0
           endif
           call dip(dipxyz)
           dipole_norm = sqrt(dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)
           call write_dipole(dipxyz, dipole_norm, 134, .false.)
       
           print*,''
           print*,' Timer: ',total_time
           print*,''
           total_time=total_time+tdstep*0.0241888
       endif

 901  format(F15.9,2x,F15.9)
!==============================================================================!


!     Performs orbital/density plots.
        if (cube_dens.or.cube_orb.or.cube_elec) then
          call g2g_timer_sum_start('cube gen')
          kk=0

          do k=1,M
            do i=1,M
              kk=kk+1
              xnano(k,i)  = X(i,M2+k)
            enddo
          enddo

          call cubegen(M15,Xnano)
          call g2g_timer_sum_stop('cube gen')
        endif


      deallocate (Y)
      deallocate (Ytrans, Xtrans, znano)
      deallocate (xnano, rmm5, rmm15)


      if (MEMO) then
        deallocate (kkind,kkinds)
        deallocate(cool,cools)
      endif

      if(timedep.eq.1) then
        call g2g_timer_sum_start('TD')
        call TD()
        call g2g_timer_sum_stop('TD')
      endif
 
#ifdef  CUBLAS
      call CUBLAS_FREE(devPtrX)
      call CUBLAS_FREE(devPtrY)
      call CUBLAS_SHUTDOWN 
#endif

!
!--------------------------------------------------------------------!
      call g2g_timer_stop('SCF')
      call g2g_timer_sum_stop('Finalize SCF')
      call g2g_timer_sum_stop('SCF')
      call g2g_timer_stop('SCF_full')
      return
      END SUBROUTINE SCF
!  -------------------------




        SUBROUTINE neighbor_list_2e()
! Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano
        USE garcha_mod, ONLY : natom, natomc, r, d, jatc, rmax, nshell, atmin, nnps, nnpp, nnpd, M, nuc
        IMPLICIT NONE
        INTEGER :: i,j, iij, iik, iikk
        REAL*8 :: zij, ti, tj, alf, rexp
          do i=1,natom
          natomc(i)=0
            do j=1,natom
              d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+ &
                 (r(i,3)-r(j,3))**2

              zij=atmin(i)+atmin(j)
              ti=atmin(i)/zij
              tj=atmin(j)/zij
              alf=atmin(i)*tj
              rexp=alf*d(i,j)
              if (rexp.lt.rmax) then
                natomc(i)=natomc(i)+1
                jatc(natomc(i),i)=j
              endif
            enddo
          enddo

          do iij=nshell(0),1,-1
            nnps(nuc(iij))=iij
          enddo

          do iik=nshell(0)+nshell(1),nshell(0)+1,-1
            nnpp(nuc(iik))=iik
          enddo

          do iikk=M,nshell(0)+nshell(1)+1,-1
            nnpd(nuc(iikk))=iikk
          enddo
        END SUBROUTINE neighbor_list_2e



	SUBROUTINE overlap_diag(fockbias, dovv,Y,Ytrans,Xtrans)
	use garcha_mod, ONLY: lowdin,M,Md,RMM,Smat,X,sqsm, nuc, natom
	USE general_module, ONLY :  sdcmp_cholesky, sdiag_canonical, vector_selection, read_list, atmorb
	IMPLICIT NONE
	LOGICAL :: docholesky
	real*8,dimension(:),  allocatable :: Dvec
	real*8,dimension(:,:),allocatable :: Vmat
        real*8,dimension(:,:),allocatable :: Ymat, Xmat, Xtrp, Ytrp
        real*8,dimension(M,M),intent(inout) :: Y,Ytrans,Xtrans
	real*8,dimension(M,M),intent(inout) :: fockbias
	logical, intent(in) :: dovv
       integer,allocatable :: atom_group(:)
       integer,allocatable :: orb_group(:)
       integer,allocatable :: orb_selection(:)
	real*8 :: weight
	INTEGER :: M13, MM, MMd !temporal hasta que rompamos RMM
	INTEGER :: iii, jjj, kk
	allocate(Vmat(M,M),Dvec(M))
        allocate(Ymat(M,M),Xmat(M,M),Xtrp(M,M),Ytrp(M,M))

        MM=M*(M+1)/2
        MMd=Md*(Md+1)/2
	M13=1 + 4*MM + 2*MMd !temporal hasta que rompamos RMM

       if (dovv.eqv..true.) then
        if (.not.allocated(atom_group)) then
          allocate(atom_group(natom))
          call read_list('atomgroup',atom_group)
        endif
        if (.not.allocated(orb_group)) then
          allocate(orb_group(M))
          call atmorb(atom_group,nuc,orb_group)
        endif
        if (.not.allocated(orb_selection)) then
          allocate(orb_selection(M))
        endif
       endif
!
! Diagonalization of S matrix, after this is not needed anymore
! S = YY^T ; X = (Y^-1)^T
! => (X^T)SX = 1
!
      docholesky=.true.
      if (lowdin) docholesky=.false.
      call g2g_timer_start('cholesky')
      call g2g_timer_sum_start('Overlap decomposition')
      IF (docholesky) THEN
         call sdcmp_cholesky(Smat,Dvec,Vmat,Ymat,Xtrp,Ytrp,Xmat)
         do iii=1,M
         do jjj=1,M
           X(iii,jjj)=Xmat(iii,jjj)
         enddo
         enddo
         Y=Ymat
         Xtrans=Xtrp
         Ytrans=Ytrp
         do kk=1,M
           RMM(M13+kk-1)=0.0d0
         enddo
      ELSE

! FFR: Canonical Diagonalization of Overlap
!--------------------------------------------------------------------!
! I am keeping Y,Ytrans and Xtrans but they should be replaced
! by the much nicer Ymat,Ytrp,Xtrp (and X by Xmat). Also, copy
! into RMM.
!
         call sdiag_canonical(Smat,Dvec,Vmat,Xmat,Xtrp,Ymat,Ytrp)
         sqsm=matmul(Vmat,Ytrp)

         if (dovv.eqv..true.) then
          fockbias=0.0d0

          weight=0.195d0
          call vector_selection(1,orb_group,orb_selection)
          call fterm_biaspot(M,sqsm,orb_selection,weight,fockbias)

          weight=-weight
          call vector_selection(2,orb_group,orb_selection)
          call fterm_biaspot(M,sqsm,orb_selection,weight,fockbias)
         endif

         do iii=1,M
         do jjj=1,M
           X(iii,jjj)=Xmat(iii,jjj)
         enddo
         enddo
         Y=Ymat
         Xtrans=Xtrp
         Ytrans=Ytrp
         do kk=1,M
           RMM(M13+kk-1)=Dvec(kk)
         enddo
      ENDIF
      call g2g_timer_stop('cholesky')
	END SUBROUTINE overlap_diag



	SUBROUTINE starting_guess(xnano)
	use garcha_mod, ONLY: RMM, ATRHO, VCINP, primera, M, X, Md, NCO
	IMPLICIT NONE
	integer :: info
	real*8, dimension (M,M), intent(inout)::xnano
	real*8, dimension (:), ALLOCATABLE :: rmm5,rmm15
	integer :: M1,M2,M3, M5, M7, M9, M11, M13, M15, M17, M18, MM, MMd !temporales hasta q rompamos RMM
	integer :: i,j,k,kk !auxiliares
	real*8 :: ff
      call g2g_timer_start('initial guess')
      call g2g_timer_sum_stop('Overlap decomposition')

      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      allocate(rmm5(MM),rmm15(mm))

      M1=1 ! first P
      M2=2*M
      M3=M1+MM ! now Pnew
      M5=M3+MM! now S, F also uses the same position after S was used
      M7=M5+MM! now G
      M9=M7+MMd ! now Gm
      M11=M9+MMd! now H
      M13=M11+MM! W ( eigenvalues ), also this space is used in least squares
      M15=M13+M! aux ( vector for ESSl)
      M17=M15+MM! Least squares
      M18=M17+MMd! vectors of MO

!
! CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
! FCe = SCe; (X^T)SX = 1
! F' = (X^T)FX
! => (X^-1*C)^-1 * F' * (X^-1*C) = e
!
! Calculate F' in RMM(M5)
      if((.not.ATRHO).and.(.not.VCINP).and.primera) then
        call g2g_timer_sum_start('initial guess')
        primera=.false.
        do i=1,M
          do j=1,M
            X(i,M+j)=0.D0
            do k=1,j
              X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+j+(M2-k)*(k-1)/2-1)
            enddo
            do k=j+1,M
              X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+k+(M2-j)*(j-1)/2-1)
            enddo
          enddo

        enddo

        kk=0
        do j=1,M
          do i=j,M
            kk=kk+1
            RMM(M5+kk-1)=0.D0
            do k=1,j
              RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
            enddo
          enddo
        enddo
!
! F' diagonalization now
! xnano will contain (X^-1)*C
!
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M13+i-1)=0.D0
        enddo
!
! ESSL OPTION
        do i=1,MM
          rmm5(i)=RMM(M5+i-1)
        enddo
        rmm15=0
        xnano=0
#ifdef  essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
! LAPACK OPTION -----------------------------------------
#ifdef pack
!
        call dspev('V','L',M,RMM5,RMM(M13),Xnano,M,RMM15,info)
#endif

        do i =1,M
          do j=1,M
            X(i,M+j)=xnano(i,j)
          enddo
        enddo
!-----------------------------------------------------------
! Recover C from (X^-1)*C
        do i=1,MM
          RMM(M5+i-1)=rmm5(i)
        enddo

        do i=1,M
          do j=1,M
            X(i,M2+j)=0.D0
            do k=1,M
              X(i,M2+j)=X(i,M2+j)+X(i,k)*X(k,M+j)
            enddo
          enddo
        enddo
      call g2g_timer_stop('initial guess')

!
! Density Matrix
!
        kk=0
        do k=1,NCO
          do i=1,M
            kk=kk+1
            RMM(M18+kk-1)=X(i,M2+k)
          enddo
        enddo
!
        kk=0
        do j=1,M
          do i=j,M
            kk=kk+1
            RMM(kk)=0.D0
!
! one factor of 2 for alpha+beta
            if(i.eq.j) then
             ff=2.D0
! another factor of 2 for direct triangular sum (j>i) w/ real basis
            else
             ff=4.D0
            endif
!
            do k=1,NCO
              RMM(kk)=RMM(kk)+ff*X(i,M2+k)*X(j,M2+k)
            enddo
          enddo
        enddo
!
        call g2g_timer_sum_stop('initial guess')
      endif
	deallocate(rmm5,rmm15)
! End of Starting guess (No MO , AO known)-------------------------------
	END SUBROUTINE starting_guess



#ifdef CUBLAS
	SUBROUTINE obtain_new_P(niter, DAMP, dovv, fockbias, good, xnano, znano, devPtrX, devPtrY)
        use cublasmath, only : cumxp_r, cumfx, cumxtf, cu_calc_fock_commuts
#else
        SUBROUTINE obtain_new_P(niter, DAMP, dovv, fockbias, good, xnano, znano, Y)
	use mathsubs, only : basechange_gemm
#endif
	use garcha_mod, ONLY: RMM, DIIS, M, Md, NCO, ndiis, hybrid_converg, SHFT, good_cut, X
	use linear_algebra, only :matrix_diagon
	IMPLICIT NONE
	INTEGER :: i,ii,j,jnuevo,k,kk,kknueva,kk2, l !auxiliares
	INTEGER, INTENT(IN) :: niter
	REAL*8, INTENT(INOUT) :: DAMP, good
	logical, intent(in) :: dovv
	real*8,dimension(M,M),intent(in) :: fockbias
	real*8, dimension (:,:), allocatable :: fock, rho, rho1, EMAT
        real*8, dimension (:,:), allocatable, save :: fockm, FP_PFm, EMAT2
	real*8, dimension (M,M), intent(inout) :: xnano, znano
       real*8, dimension (:), ALLOCATABLE, save :: bcoef
       real*8, dimension (:), ALLOCATABLE :: suma
	REAL*8, dimension(1000) :: WORK, IWORK
        REAL*8,ALLOCATABLE :: WORK2(:)
        INTEGER, ALLOCATABLE :: IWORK2(:)
	INTEGER :: LIWORK  !hay algo mal con LIWORK e IWORK, Nick

	integer :: M1,M2,M3, M5, M7, M9, M11, M13, M15, M17, M18, MM, MMd !temporales hasta q rompamos RMM
	integer :: ndiist
	logical, save :: hagodiis
	INTEGER :: INFO, LWORK
        real*8, dimension (M,M) :: scratch
!FFR!
       real*8,allocatable :: eigen_vecs(:,:),eigen_vals(:)

#ifdef  CUBLAS
        integer*8 :: devPtrX, devPtrY !ver intent
!        real*8, dimension (M,M) :: scratch
#else
	real*8, dimension (M,M) :: Y
        real*8, dimension (M,M) :: scratch1
#endif

        if(niter.le.ndiis) then
          ndiist=niter
        else
          ndiist=ndiis
        endif

	allocate (fock(M,M), rho(M,M), rho1(M,M))

	MM=M*(M+1)/2
	MMd=Md*(Md+1)/2
      M1=1 ! first P
      M2=2*M
      M3=M1+MM ! now Pnew
      M5=M3+MM! now S, F also uses the same position after S was used
      M7=M5+MM! now G
      M9=M7+MMd ! now Gm
      M11=M9+MMd! now H
      M13=M11+MM! W ( eigenvalues ), also this space is used in least squares
      M15=M13+M! aux ( vector for ESSl)
      M17=M15+MM! Least squares
      M18=M17+MMd! vectors of MO


!
!
! now, we know S matrix, and F matrix, and E for a given P
! 1) diagonalize S, get X=U s^(-1/2)
! 2) get U F Ut
! 3) diagonalize F
! 4) get vec ( coeff) ---->  P new
! 5) iterate
! call diagonalization routine for S , get after U s^(-1/2)
! where U matrix with eigenvectors of S , and s is vector with
! eigenvalues
!
! here in RMM(M5) it is stored the new Fock matrix
! test damping on Fock matrix
!
!      DAMP=gold
        if (niter.eq.1) then
          DAMP=0.0D0
	  hagodiis=.false.
	  if (DIIS ) then ! agregado para cambio de damping a diis, Nick
	    if (.not. allocated(fockm)) allocate(fockm(MM,ndiis))
            if (.not. allocated(FP_PFm)) allocate(FP_PFm(MM,ndiis))
            if (.not. allocated(bcoef)) allocate(bcoef(ndiis+1))
          endif
       endif

	if (DIIS) then ! agregado para cambio de damping a diis, Nick
            allocate(suma(MM))
        endif



!-----------------------------------------------------------------------------------------
! If DIIS is turned on, update fockm with the current transformed F' (into ON
! basis) and update FP_PFm with the current transformed [F',P']
!-----------------------------------------------------------------------------------------
        if (DIIS) then! agregado para cambio de damping a diis, Nick
          call g2g_timer_sum_start('DIIS')
          call g2g_timer_sum_start('DIIS prep')
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Expand F into square form
! (for better memory access patterns in coming multiplications)
!-----------------------------------------------------------------------------------------
          do j=1,M
            do k=1,j
              fock(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
            enddo
            do k=j+1,M
              fock(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
            enddo
          enddo

! FFR: Van Voorhis Term for DIIS
!--------------------------------------------------------------------!
         if (dovv.eqv..true.) fock=fock+fockbias

!-----------------------------------------------------------------------------------------
! Expand density matrix into full square form (before, density matrix was set up for triangular sums
! (sum j>=i) so off-diagonal elements need to be divided by 2 to get square-form numbers)
! (for better memory access patterns in coming multiplications)
!-----------------------------------------------------------------------------------------
          do j=1,M
            do k=1,j-1
              rho(j,k)=(RMM(j+(M2-k)*(k-1)/2))/2
            enddo
            rho(j,j)=RMM(j+(M2-j)*(j-1)/2)
            do k=j+1,M
              rho(j,k)=RMM(k+(M2-j)*(j-1)/2)/2
            enddo
          enddo

          ! Calculate F' and [F',P']
          call g2g_timer_start('commutators + basechange')
#ifdef  CUBLAS
          call cu_calc_fock_commuts(fock,rho,devPtrX,devPtrY,scratch,M)
#else
          call calc_fock_commuts(fock,rho,X,Y,scratch,scratch1,M)
#endif
          call g2g_timer_stop('commutators + basechange')
          ! update fockm with F'
          do j=ndiis-(ndiist-1),ndiis-1
            do i=1,MM
              fockm(i,j)=fockm(i,j+1)
            enddo
          enddo
          do k=1,M
            do j=k,M
              i=j+(M2-k)*(k-1)/2
              fockm(i,ndiis)=fock(j,k)
            enddo
          enddo
!-----------------------------------------------------------------------------------------
! now, scratch = A = F' * P'; scratch1 = A^T
! [F',P'] = A - A^T
!-----------------------------------------------------------------------------------------
          ! update FP_PFm with [F',P']
          do j=ndiis-(ndiist-1),ndiis-1
            do i=1,MM
              FP_PFm(i,j)=FP_PFm(i,j+1)
            enddo
          enddo
#ifdef  CUBLAS
          do k=1,M
            do j=k,M
              i=j+(M2-k)*(k-1)/2
              FP_PFm(i,ndiis)=scratch(j,k)
            enddo
          enddo
#else
          do k=1,M
            do j=k,M
              i=j+(M2-k)*(k-1)/2
              FP_PFm(i,ndiis)=scratch(j,k)-scratch1(j,k)
            enddo
          enddo
#endif
          call g2g_timer_sum_pause('DIIS prep')
          call g2g_timer_sum_pause('DIIS')
        endif

!
!-------------Decidiendo cual critero de convergencia usar-----------
!-----------------------------------------------------------------------------------------
! iF DIIS=T
! Do simple damping 2nd iteration; DIIS afterwards
!-----------------------------------------------------------------------------------------
! IF DIIS=F and hybrid_converg T
! Do Do simple damping until GOOD<=GOOD_CUT. DIIS afterwards
!-----------------------------------------------------------------------------------------
! IF DIIS=F and hybrid_converg F
! Always do damping (after first iteration)
!-----------------------------------------------------------------------------------------

      if ( .not. hybrid_converg ) then
              if (niter.gt.2.and.(DIIS)) then
                hagodiis=.true.
            end if
      else ! change damping to diis, Nick
            if (good .lt. good_cut ) then
              if ( .not. hagodiis ) then
                write(6,*)
                write(6,8503)
                write(6,8504) niter
                write(6,8505)
                write(6,*)
              end if
              hagodiis=.true.
            end if
      end if

!-----------------------------------------------------------------------------------------
! If we are not doing diis this iteration, apply damping to F, save this
! F in RMM(M3) for next iteration's damping and put F' = X^T * F * X in RMM(M5)
!----------------------------------------------------------------------------------------
        if(.not.hagodiis) then
          call g2g_timer_start('Fock damping')


          if(niter.ge.2) then
            do k=1,MM
              kk=M5+k-1
              kk2=M3+k-1
              RMM(kk)=(RMM(kk)+DAMP*RMM(kk2))/(1.D0+DAMP)
            enddo
          endif
! the newly constructed damped matrix is stored, for next iteration
! in RMM(M3)
!
         do k=1,MM
            kk=M5+k-1
            kk2=M3+k-1
            RMM(kk2)=RMM(kk)
          enddo
!-------------------------------------------------------------!
            fock=0
            do j=1,M
              do k=1,j
                 fock(k,j)=RMM(M5+j+(M2-k)*(k-1)/2-1)
              enddo
              do k=j+1,M
                 fock(k,j)=RMM(M5+k+(M2-j)*(j-1)/2-1)
              enddo
            enddo

! FFR: Van Voorhis Term for not DIIS
!--------------------------------------------------------------------!
         if (dovv.eqv..true.) fock=fock+fockbias

#ifdef  CUBLAS
            call cumxtf(fock,devPtrX,fock,M)
            call cumfx(fock,DevPtrX,fock,M)
#else
            fock=basechange_gemm(M,fock,x)
#endif
          do j=1,M
             do k=1,j
                RMM(M5+j+(M2-k)*(k-1)/2-1)=fock(j,k)
             enddo
             do k=j+1,M
                RMM(M5+k+(M2-j)*(j-1)/2-1)=fock(j,k)
             enddo
          enddo
          call g2g_timer_stop('Fock damping')
        endif


!
! now F contains transformed F
! diagonalization now
!
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M13+i-1)=0.D0
        enddo

!---- LEVEL SHIFT
!
        if (SHFT) then
!       shi=shi*0.99
! adition of level shifts
! constant to diagonal (virtual) elements
          do i=NCO+1,M
            ii=i+(i-1)*(M2-i)/2
            RMM(M5+ii-1)=RMM(M5+ii-1)! +shi
! Este shi no esta definido en ningun lado, sirve para algo, Nick
          enddo
        endif

        call g2g_timer_stop('actualiza rmm')


!----------Si hagodiis(ver mas arriba) es true entonces sigo-----------------------
        call g2g_timer_start('diis')
!--------Pasar columnas de FP_PFm a matrices y multiplicarlas y escribir EMAT-------


        if(DIIS) then
          call g2g_timer_sum_start('DIIS')
          call g2g_timer_sum_start('DIIS Fock update')
!          deallocate(EMAT)
          allocate(EMAT(ndiist+1,ndiist+1))
! Before ndiis iterations, we just start from the old EMAT
          if(niter.gt.1.and.niter.le.ndiis) then
            EMAT=0
            do k=1,ndiist-1
              do kk=1,ndiist-1
                EMAT(K,KK)=EMAT2(K,KK)
              enddo
            enddo
            deallocate (EMAT2)
! After ndiis iterations, we start shifting out the oldest iteration stored
          elseif(niter.gt.ndiis) then
            do k=1,ndiist-1
              do kk=1,ndiist-1
                EMAT(K,KK)=EMAT2(K+1,KK+1)
              enddo
            enddo
            deallocate (EMAT2)
          endif

          k=ndiist
          do kk=1,ndiist
            kknueva=kk+(ndiis-ndiist)
!-------Escribimos en xnano y znano dos conmutadores de distintas iteraciones------
            do i=1,M
              do j=1,i
                xnano(i,j)=FP_PFm(i+(M2-j)*(j-1)/2,ndiis)
                znano(i,j)=FP_PFm(i+(M2-j)*(j-1)/2,kknueva)
              enddo
              do j=i+1,M
                xnano(i,j)=FP_PFm(j+(M2-i)*(i-1)/2,ndiis)
                znano(i,j)=FP_PFm(j+(M2-i)*(i-1)/2,kknueva)
              enddo
            enddo
!#ifdef CUBLAS (Only diagonal elements must be computed, so the entire multiplication is a waste..)
!                   call cumatmul_r(xnano,znano,rho1,M)
!#else
            call matmuldiag(xnano,znano,rho1,M)
!#endif

            EMAT(ndiist,kk)=0.
            if(kk.ne.ndiist) EMAT(kk,ndiist)=0.
            do l=1,M
              EMAT(ndiist,kk)=EMAT(ndiist,kk)+rho1(l,l)
              if (kk.ne.ndiist) then
                EMAT(kk,ndiist)=EMAT(ndiist,kk)
              endif
            enddo
          enddo

          do i=1,ndiist
            EMAT(i,ndiist+1)= -1.0
            EMAT(ndiist+1,i)= -1.0
          enddo
          EMAT(ndiist+1, ndiist+1)= 0.0


          if (.not. allocated(EMAT2)) allocate(EMAT2(ndiist+1,ndiist+1))


          EMAT2=EMAT

!        ematalloct=.true.
!********************************************************************
!   THE MATRIX EMAT SHOULD HAVE FORM
!
!      |<E(1)*E(1)>  <E(1)*E(2)> ...   -1.0|
!      |<E(2)*E(1)>  <E(2)*E(2)> ...   -1.0|
!      |<E(3)*E(1)>  <E(3)*E(2)> ...   -1.0|
!      |<E(4)*E(1)>  <E(4)*E(2)> ...   -1.0|
!      |     .            .      ...     . |
!      |   -1.0         -1.0     ...    0. |
!
!   WHERE <E(I)*E(J)> IS THE SCALAR PRODUCT OF [F*P] FOR ITERATION I
!   TIMES [F*P] FOR ITERATION J.
!
!********************************************************************
!-----Pasamos a resolver con DGELS el problema EMAT*ci=bcoef-------
          if (hagodiis) then
            do i=1,ndiist
              bcoef(i)=0
            enddo
            bcoef(ndiist+1)=-1

!----------Cálculo de parámetro optimo para DGELS-------------------
!
!
        LWORK = -1
      CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
      ndiist+1, bcoef, ndiist+1, WORK, LWORK, INFO )
        LWORK = MIN( 1000, INT( WORK( 1 ) ) )


!-----Resuelve la ecuación A*X = B. (EMAT*ci=bcoef). La solución la escribe en bcoef------

      CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
       ndiist+1, bcoef, ndiist+1, WORK, LWORK, INFO )
!--------Construccion de la "nueva" matriz de fock como cl de las anteriores--------------
!--------Eventualmente se puede probar con la matriz densidad-----------------------------
            suma=0
            do j=1,ndiist
              jnuevo=j+(ndiis-ndiist)
              do i=1,MM
                suma(i)=suma(i)+bcoef(j)*fockm(i,jnuevo)
              enddo
            enddo
            do i=1,MM
              RMM(M5+i-1)=suma(i)
            enddo
            call g2g_timer_stop('diis')
        endif
        call g2g_timer_sum_pause('DIIS Fock update')
        call g2g_timer_sum_pause('DIIS')
      endif

!
! F' diagonalization now
! X(1,M) will contain (X^-1)*C
!

       call g2g_timer_start('dspev')
       call g2g_timer_sum_pause('SCF acceleration')
       call g2g_timer_sum_start('diagonalization')
! ESSL OPTION ---------------------------------------------------
#ifdef essl
       call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif

       if(.not.allocated(fock)) allocate (fock(M,M))
       fock=0
       do j=1,M
       do k=1,j
         i=j+(M2-k)*(k-1)/2
         fock(j,k)=RMM(M5+i-1)
       enddo
       enddo

!

      if(.not.allocated(fock)) allocate (fock(M,M))
      fock=0
      do j=1,M
        do k=1,j
         i=j+(M2-k)*(k-1)/2
         fock(j,k)=RMM(M5+i-1)
        enddo
      enddo
!
! FFR: Fock Diagonaliation is now external
!--------------------------------------------------------------------!
       call g2g_timer_start('scf - fock diagonalization')

       if ( allocated(eigen_vecs) ) deallocate(eigen_vecs)
       allocate( eigen_vecs(M,M) )
       if ( allocated(eigen_vals) ) deallocate(eigen_vals)
       allocate( eigen_vals(M) )

       call matrix_diagon( fock, eigen_vecs, eigen_vals )
       fock=eigen_vecs
       do kk=1,M
         RMM(M13+kk-1) = eigen_vals(kk)
       end do


       call g2g_timer_stop('scf - fock diagonalization')
!
!
! FFR: Old version - kept until functionality is assured
!--------------------------------------------------------------------!

!######################################################################################
!####################    PREGUNTAR  ###################################################
!######################################################################################


      if (.false.) then   !Esto lo seguimos dejando o ya se puede sacar???, Nick
       LWORK=-1
#ifdef  magma
       call magmaf_dsyevd('V','L',M,fock,M,RMM(M13),WORK,LWORK &
       ,IWORK,LWORK,info)
#else
       call dsyevd('V','L',M,fock,M,RMM(M13),WORK,LWORK &
       ,IWORK,LWORK,info)
#endif
       LWORK=work(1)          
       LIWORK=IWORK(1)
       if(allocated(WORK2)) deallocate (WORK2)
       if(allocated(IWORK2)) deallocate (IWORK2)
       allocate (WORK2(LWORK),IWORK2(LIWORK))

#ifdef  magma
       call magmaf_dsyevd('V','L',M,fock,M,RMM(M13),WORK2,LWORK &
       ,IWORK2,LIWORK,info)
#else
       call dsyevd('V','L',M,fock,M,RMM(M13),WORK2,LWORK &
       ,IWORK2,LIWORK,info)
#endif
      end if

!######################################################################################
!######################################################################################
!######################################################################################

!
!--------------------------------------------------------------------!

       call g2g_timer_stop('dspev')
       call g2g_timer_sum_pause('diagonalization')

       call g2g_timer_start('coeff')
       call g2g_timer_sum_start('MO coefficients')

!-----------------------------------------------------------
!
! diagonalization now
!-----------------------------------------------------------
! Recover C from (X^-1)*C; put into xnano
!-----------------------------------------------------------
#ifdef CUBLAS
        call cumxp_r(fock,devPtrX,xnano,M)
        do i=1,M
          do j=1,M
             X(i,M2+j)=xnano(i,j)
          enddo
       enddo
#else
! new coefficients
        do i=1,M
           do j=1,M
              xnano(i,j)=x(i,j)
           enddo
        enddo
        xnano=matmul(xnano,fock)
        do i=1,M
          do j=1,M
             X(i,M2+j)=xnano(i,j)
          enddo
       enddo
#endif
      call g2g_timer_stop('coeff')
      call g2g_timer_start('otras cosas')
      call g2g_timer_sum_pause('MO coefficients')
      call g2g_timer_sum_start('new density')
!
! --- For the first iteration, damping on density matrix
! Important for the case of strating guess of AO
!
      kk=0
      do k=1,NCO
        do i=1,M
          kk=kk+1
          RMM(M18+kk-1)=X(i,M2+k)
          xnano(k,i)  = X(i,M2+k)
        enddo
      enddo
!
! Construction of new density matrix and comparison with old one
      kk=0

      good=0.

      call g2g_timer_start('dens_GPU')
      call density(M,NCO,X,xnano)


      if(DIIS) deallocate (fock, rho, rho1)


!%%%%%%%%%%%%%%%%%%%%%%%%% Nuevos Formatos, Nick %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 8503 FORMAT(4x,"╔════════════════", &
      "══╦══════╦═══════╗")
 8504 FORMAT(4x,"║ Changing to DIIS ║ step ║",2x,i4,1x,"║")
 8505 FORMAT(4x,"╚════════════════", &
      "══╩══════╩═══════╝")
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	END SUBROUTINE obtain_new_P



	SUBROUTINE COPY_VEC(VEC,DIM_VEC,POINTER_RMM)
!subrutina temporal para empezar a romper RMM
!copia el vector VEC a RMM posicion POINTER_RMM
	use garcha_mod, ONLY: RMM
	IMPLICIT NONE
	integer, intent(in) :: DIM_VEC,POINTER_RMM
	real*8, dimension(DIM_VEC), intent(in) :: VEC
	integer :: i
	do i=1, DIM_VEC
	   RMM(POINTER_RMM+i-1)=VEC(i)
	end do
	END SUBROUTINE COPY_VEC


