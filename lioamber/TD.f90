!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
         SUBROUTINE TD()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! REAL TIME-TDDFT
!
! Dario Estrin, 1992
! Nano, Dario, Uriel, Damian 2012
!
!  This subrutine takes the converged density matrix from an SCF calculation
!  and evolves it in time. In the input file the total number of propagation
!  steps is specified (nstep) as well as the time of each evolution step
!  (tdstep).
!  This implementation has two alternatives to evolve the density in time. The
!  first one (propagator=1) is the Verlet algorithm that uses a convination of
!  Liouville von Newmann expresion for the time derivative of the density matrix
!  and a first order Taylor expansion of the density matrix. The second one
!  (propagator=2) is the Magnus propagation scheme that uses Backer Campbell
!  Hausdorff (BCH) formula. For this reason when Magnus is used the number of
!  total conmutators in the BCH espansion has to be specified (NBCH, default=10).
!  A narrow gaussian type electric field can be introduced during the time
!  evolution in order to excite all electronic frequencies with the same intensity.
!  Once this perturbation is turned on (Field=t, exter=t) each component of the
!  external electric field has to be specified in the input file (Fx,Fy,Fz).
!  In each step of the propagation the cartesian components of the sistems dipole
!  are stored in files x.dip, y.dip, z.dip.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!       USE latom
        USE garcha_mod, only : M, Md, NBCH, propagator, tdstep, idip,          &
                               tdrestart, exists, RMM, NCO, nang, Iz, natomc,  &
                               r, d, atmin, rmax, jatc, nshell, nnps, nnpp,    &
                               nnpd, igrid2, Nuc, predcoef, npas, nsol, pc, X, &
                               Smat, MEMO, ntdstep, field, exter, epsilon,     &
                               writedens, a0, sol, kkind, kkinds, cool, cools, &
                               GRAD, natom, sqsm, Fx, Fy, Fz
       use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC
       use mathsubs
       use transport
       use general_module, ONLY :  atmorb
       use faint_cpu77, only: int1, intsol, int3mem, intfld, int3lu
#ifdef CUBLAS
        use cublasmath
#endif
       implicit none
       !IMPLICIT REAL*8 (a-h,o-z)
       real*8 :: dipxyz(3), dipole_norm,Qc,Qc2,zij,ti,tj,alf,&
                 rexp,E1s,Ens,Ex,g,factor,fxx,fyy,fzz,Es,ff,t0
       INTEGER :: istep,nopt,M1,M2,M3,i,j,M5,M7,MM,MM2,&
                  MMd,Md2,k,M13,M15,M17,M9,M20,M18,M19,M11,Nel,&
                  Nunp,igpu,info,kk,n,unit1,unit2
       REAL*8 :: t,E2,E,En,E1
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: &
        xnano2,xmm,xtrans,ytrans,Y,fock,&
        F1a,F1b,overlap,rhoscratch
       real*8, dimension (:,:), ALLOCATABLE :: elmu

#ifdef TD_SIMPLE
       COMPLEX*8 :: Im,Ix
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) :: &
        rho,rhonew,rhold,xnano,rho1
#else
       COMPLEX*16 :: Im,Ix
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) :: &
        rho,rhonew,rhold,xnano,rho1
#endif
       REAL*8,DIMENSION(natom) ::  q
       REAL*8,dimension(:),ALLOCATABLE :: factorial
       INTEGER            :: LWORK,ii,jj
       REAL*8,ALLOCATABLE :: WORK(:)
!!------------------------------------!!
!! FFR ADD
       INTEGER :: &
        pert_steps,lpfrg_steps,chkpntF1a,chkpntF1b
       REAL*8 :: &
        dt_magnus,dt_lpfrg
        logical :: ematalloct
!! CUBLAS
#ifdef CUBLAS
      integer sizeof_real
      parameter(sizeof_real=8)
      integer sizeof_complex
#ifdef TD_SIMPLE
      parameter(sizeof_complex=8)
#else
      parameter(sizeof_complex=16)
#endif
      integer stat
      integer*8 devPtrX, devPtrY,devPtrXc
      external CUBLAS_INIT, CUBLAS_SET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_GET_MATRIX
      external CUBLAS_FREE
     ! integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX,CUBLAS_GET_MATRIX  sacar el coment
#endif
!!   GROUP OF CHARGES
       INTEGER             :: ngroup
       INTEGER,ALLOCATABLE :: group(:)
       REAL*8,ALLOCATABLE  :: qgr(:)
       REAL*8 :: tiempo1000

!FFR
       logical             :: dovv
       integer,allocatable :: orb_group(:)
       integer,allocatable :: orb_selection(:)

       real*8,dimension(:,:),allocatable :: Vmat
       real*8,dimension(:),allocatable   :: Dvec


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       call g2g_timer_start('TD')
       call g2g_timer_start('inicio')
       open(unit = 134, file = "dipole_moment_td")
       ALLOCATE(factorial(NBCH))
!!------------------------------------!!
        nopt=0

#ifdef CUBLAS
       write(*,*) 'USING CUBLAS'
!        las dos alloc eran con stat
       call CUBLAS_INIT()
       call CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
       call CUBLAS_ALLOC(M*M, sizeof_real, devPtrY)
       if (stat.NE.0) then
           write(*,*) "initialization failed -TD"
           call CUBLAS_SHUTDOWN
           stop
       endif
#endif
#ifdef TD_SIMPLE
        write(*,*) 'simple presition complex'
#else
        write(*,*) 'double presition complex'
#endif
       if(propagator.eq.2) then
          dt_magnus=tdstep
          dt_lpfrg=tdstep*0.10D0
          factorial(1)=1.0D0
#ifdef CUBLAS
          DO ii=1,NBCH
             factorial(ii)=1.0D0/ii
          ENDDO
#else
       DO ii=2,NBCH
         factorial(ii)=factorial(ii-1)/ii
       ENDDO
#endif
       endif
       if(propagator.eq.1) then
          dt_lpfrg=tdstep
       endif
!!------------------------------------!!
!! FFR ADD:
       pert_steps=100
       lpfrg_steps=200
       chkpntF1a=185
       chkpntF1b=195
!--------------------------------------------------------------------!
! Pointers -
       !Ndens=1
       E=0.0D0
       E1=0.0D0
       En=0.0D0
       E2=0.0D0
       idip=1
       !ngeo=ngeo+1
       Im=(0.0D0,2.0D0)
       !sq2=sqrt(2.D0)
       MM=M*(M+1)/2
       MM2=M**2
       MMd=Md*(Md+1)/2
       Md2=2*Md
       M2=2*M
!
       ALLOCATE(xnano(M,M),xnano2(M,M),fock(M,M),rhonew(M,M),&
        rhold(M,M),rho(M,M),xmm(M,M),xtrans(M,M),Y(M,M),ytrans(M,M),&
        rho1(M,M), rhofirst(M,M))
!
!LOWDIN POPULATION: temporarily off
      lpop=.false.


!TRANSPORT------------------------------------------------------------!

      IF (TRANSPORT_CALC) THEN
         ALLOCATE(mapmat(M,M))
      ENDIF

      if (transport_calc) then
         inquire(file='atomgroup',exist=exists)
         if (.not.exists) then
            write(*,*) 'ERROR CANNOT FIND atomgroup file'
            stop
         endif
            open(unit=678,file='atomgroup')
            allocate(group(natom))
            ngroup=0
            do n=1,natom
               read(678,*) kk
               group(n)=kk
               if (kk.gt.ngroup) ngroup=kk
            enddo
            allocate(qgr(ngroup))
            if(ngroup.gt.3) write(*,*) 'if the number of group &
      is greater than 3 then group 1 should be the donor and 2 the &
      acceptor'
            close(unit=678)
      endif
      IF(TRANSPORT_CALC) THEN
         call mat_map(group, mapmat, Nuc, M, natom)
      ENDIF
!--------------------------------------------------------------------!
      if(propagator.eq.2) allocate (F1a(M,M),F1b(M,M))
!--------------------------------------------------------------------!
      if (tdrestart) then
         inquire(file='rho.restart',exist=exists)
         if (.not.exists) then
             write(*,*) 'ERROR CANNOT FIND rho.restart'
             write(*,*) '(if you are not restarting a previous &
      run set tdrestart= false)'
             stop
         endif
         open(unit=1544,file='rho.restart',status='old')
         do j=1,M
            do k=1,M
               read(1544,*) rho(j,k)
            enddo
         enddo
         do j=1,M
            do k=j,M
               if(j.eq.k) then
                  RMM(k+(M2-j)*(j-1)/2)=REAL(rho(j,k))
               else
                  RMM(k+(M2-j)*(j-1)/2)=(REAL(rho(j,k)))*2
               endif
            enddo
         enddo
         if (propagator .eq. 2) then
            inquire(file='F1a.restart',exist=exists)
            if (.not.exists) then
               write(*,*) 'ERROR CANNOT FIND F1a.restart'
               write(*,*) '(if you are not restarting a &
      previous run set tdrestart= false)'
               stop
            endif
            inquire(file='F1b.restart',exist=exists)
            if (.not.exists) then
               write(*,*) 'ERROR CANNOT FIND F1b.restart'
               write(*,*) '(if you are not restarting a &
      previous run set tdrestart= false)'
               stop
            endif
            open(unit=7777,file='F1a.restart',status='old')
            do i=1,M
               do j=1,M
                  read(7777,*) F1a(i,j)
               enddo
            enddo
            open(unit=7399,file='F1b.restart',status='old')
            do i=1,M
               do j=1,M
                  read(7399,*) F1b(i,j)
               enddo
            enddo
         endif
!--------------------------------------------------------------------!
! We read the density matrix stored in RMM(1,2,3,...,MM) and it is copied in rho matrix.
         else
!          do j=1,M
!             do k=1,j-1
!                rho(j,k)=RMM(j+(M2-k)*(k-1)/2)/2
!             enddo
!             rho(j,j)=RMM(j+(M2-k)*(j-1)/2)
!             do k=j+1,M
!                rho(j,k)=RMM(k+(M2-j)*(j-1)/2)/2
!             enddo
!           enddo
            call spunpack_rtc('L',M,RMM,rho)
         endif
!------------------------------------------------------------------------------!
! first i
            M1=1
! now Fold
            M3=M1+MM
! now S, F also uses the same position after S was used
            M5=M3+MM
! now G
            M7=M5+MM
! now Gm
            M9=M7+MMd
! now H
            M11=M9+MMd
! W ( eigenvalues ), also this space is used in least squares
            M13=M11+MM
! aux ( vector for ESSl)
            M15=M13+M
! Least squares
            M17=M15+MM
! vectors of MO
            M18=M17+MMd
! weights (in case of using option )
            M19=M18+M*NCO
! RAM storage of two-electron integrals (if MEMO=T)
            M20 = M19 + natom*50*nang
!
            Nel=2*NCO+Nunp
! Initializations/Defaults
       write(*,*) ' TD CALCULATION  '
!--------------------------------------!
           !niter=0
           !D1=1.D0
           !D2=1.D0
!--------------------------------------!
           Qc=0.0D0
           do i=1,natom
             Qc=Qc+Iz(i)
           enddo
           Qc=Qc-Nel
           Qc2=Qc**2

! FFR: Variable Allocation
!--------------------------------------------------------------------!
       if ( allocated(Vmat) ) deallocate(Vmat)
       if ( allocated(Dvec) ) deallocate(Dvec)
       if ( allocated(sqsm) ) deallocate(sqsm)

       allocate(Vmat(M,M),Dvec(M))
       allocate(sqsm(M,M))

       if (transport_calc) then
        if (.not.allocated(orb_group)) then
          allocate(orb_group(M))
          call atmorb(group,nuc,orb_group)
        endif
        if (.not.allocated(orb_selection)) then
          allocate(orb_selection(M))
        endif
       endif

!------------------------------------------------------------------------------!
! Two electron integral with neighbor list.
!
            do i=1,natom
              natomc(i)=0
!
              do j=1,natom
                d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+&
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
            do ii=nshell(0),1,-1
              nnps(Nuc(ii))=ii
            enddo
            do ii=nshell(0)+nshell(1),nshell(0)+1,-1
              nnpp(Nuc(ii))=ii
            enddo
            do ii=M,nshell(0)+nshell(1)+1,-1
              nnpd(Nuc(ii))=ii
            enddo
!------------------------------------------------------------------!
!
! Create integration grid for XC here
! Assign points to groups (spheres/cubes)
! Assign significant functions to groups
! -Calculate point weights
!
      call g2g_timer_sum_start('Exchange-correlation grid setup')
      call g2g_reload_atom_positions(igrid2)
      call g2g_timer_sum_stop('Exchange-correlation grid setup')

      call aint_query_gpu_level(igpu)
      if (igpu.gt.1) call aint_new_step()

      if (predcoef.and.npas.gt.3) then

!        if (.not.OPEN) then
!          if(verbose) write(*,*) 'prediciendo densidad'
!          do i=1,MM
!            RMM(i)=(3*old1(i))-(3*old2(i))+(old3(i))
!          enddo
!         endif
       endif
!------------------------------------------------------------------------------!
! H H core, 1 electron matrix elements
      call g2g_timer_sum_start('1-e Fock')
      call g2g_timer_sum_start('Nuclear attraction')
      call int1(En)

      if (ecpmode) then
          write(*,*) "agrego terminos AAA,AAB a los de 1e"
          do k=1,MM
               term1e(k)=RMM(M11+k-1)
!copia los terminos de 1e
!               write(89,*) RMM(M11+k-1),VAAA(k),VAAB(k)
               RMM(M11+k-1)=RMM(M11+k-1)+VAAA(k)+VAAB(k)+VBAC(k)
!               write(89,*) RMM(M11+k-1)
!agrega el ECP AAA a los terminos de 1 e
          enddo
      end if

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
!--------------------------------------!
            E1=0.D0
            do k=1,MM
              E1=E1+RMM(k)*RMM(M11+k-1)
            enddo
            call g2g_timer_sum_stop('1-e Fock')

            if( transport_calc )then
              allocate(overlap(M,M),rhoscratch(M,M))
              call spunpack('L',M,RMM(M5),overlap)
            endif
!--------------------------------------!
! Diagonalization of S matrix, after this is not needed anymore
! s is in RMM(M13,M13+1,M13+2,...,M13+MM)
!--------------------------------------!
! ESSL OPTION
#ifdef essl
       call DSPEV(1,RMM(M5),RMM(M13),X,M,M,RMM(M15),M2)
#endif
!--------------------------------------!
! LAPACK OPTION
#ifdef pack
       do ii=1,M; do jj=1,M
         X(ii,jj)=Smat(ii,jj)
       enddo; enddo
       if (allocated(WORK)) deallocate(WORK); allocate(WORK(1))
       call dsyev('V','L',M,X,M,RMM(M13),WORK,-1,info)
       LWORK=int(WORK(1));  deallocate(WORK); allocate(WORK(LWORK))
       call dsyev('V','L',M,X,M,RMM(M13),WORK,LWORK,info)
#endif
!--------------------------------------!
! Here, we obtain the transformation matrices X and Y for converting
! from the atomic orbital to a molecular orbital basis (truncated
! during linear dependency elimination).
! S is the overlap matrix
! s is the diagonal eigenvalue matrix of S
! U is the eigenvector matrix of S
! X=U s^(-1/2)
! matrix X's dimension is M*3M. In the first M*M terms it contains
! the transformation matrices and in the other M*2M terms it contains
! auxiliar matrices.
            call g2g_timer_start('inicio1')
            do j=1,M
              if (RMM(M13+j-1).lt.1.0D-06) then
                write(*,*) 'LINEAR DEPENDENCY DETECTED'
                do i=1,M
                  X(i,j)=0.0D0
                  Y(i,j)=0.0D0
                enddo
              else
                do i=1,M
                  Y(i,j)=X(i,j)*sqrt(RMM(M13+j-1))
                  X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))
                enddo
              endif
            enddo
!------------------------------------------------------------------------------!
! Here rho1 is used as an auxiliar matrix. Then we need to restore its original value
#ifdef CUBLAS
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=cmplx(X(i,j),0.0D0)
               ENDDO
            ENDDO
!           estas tres tenian stat
            call CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
            call CUBLAS_ALLOC(M*M, sizeof_complex, devPtrXc)
            call CUBLAS_ALLOC(M*M, sizeof_complex, devPtrY)
            if (stat.NE.0) then
            write(*,*) "X and/or Y memory allocation failed"
            call CUBLAS_SHUTDOWN()
            stop
            endif
!           dos con stat
            call CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrXc,M)
            call CUBLAS_SET_MATRIX(M,M,sizeof_real,x,M,devPtrX,M)
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=cmplx(Y(i,j),0.0D0)
               ENDDO
            ENDDO
!           ponerle stat
            call CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrY,M)
            if (stat.NE.0) then
            write(*,*) "X and/or Y setting failed"
            call CUBLAS_SHUTDOWN
            stop
            endif
            rho1=0D0
#endif
!------------------------------------------------------------------------------!
! the transformation matrices is copied in xmm
!
            do i=1,M
               do j=1,M
                  xmm(i,j)=X(i,j)
               enddo
            enddo
! the tranposed matrixes are calculated
            do i=1,M
               do j=1,M
                 xtrans(j,i)=X(i,j)
                 ytrans(j,i)=Y(i,j)
               enddo
            enddo
!------------------------------------------------------------------------------!
! External Electric Field components
!
!       write(*,*) 'fx =', fx
!       write(*,*) 'fy =', fy
!       write(*,*) 'fz =', fz
!------------------------------------------------------------------------------!

!Transport --------------------------------------------------------------------!
      IF(TRANSPORT_CALC) THEN
         open(unit=100000,file='rhofirst')
         IF(generate_rho0) then
            DO i=1,M
               DO j=1,M
                  write(100000,*) rho(i,j)
               ENDDO
            ENDDO
            rhofirst=rho
            write(*,*) 'HE ESCRITO RHOFIRST'
         else
            DO i=1,M
               DO j=1,M
                  read(100000,*) rhofirst(i,j)
               ENDDO
            ENDDO
            write(*,*) ' HE LEIDO RHOFIRST '
         ENDIF
      ENDIF
!-----------------------------------------------------------------------------!

! Rho is transformed to the orthonormal basis
! with matmul:
#ifdef CUBLAS
       call g2g_timer_start('complex_rho_ao_to_on-cu')
       rho1=basechange_cublas(M,rho,devPtrY,'dir')
       rho=rho1
       call g2g_timer_stop('complex_rho_ao_to_on-cu')
#else
       rho=matmul(ytrans,rho)
       rho=matmul(rho,y)
#endif

!Transport: trace writing---------------------------------------------------!

      IF (TRANSPORT_CALC) THEN
         DO i=1,M
            traza0=traza0+Real(rho(i,i))
         ENDDO
         write(*,*) 'traza0 =', traza0
      ENDIF
!---------------------------------------------------------------------------!

! with matmulnanoc
! (NO LONGER AVAILABLE; USE BASECHANGE INSTEAD)
!            call matmulnanoc(rho,Y,rho1,M)
!            rho=rho1
!            rho=basechange(M,Ytrans,rho,Y)
!--------------------------------------!
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
#ifdef CUBLAS
          if(.not.TRANSPORT_CALC) call CUBLAS_FREE(devPtrY)
#endif

            call g2g_timer_stop('inicio')

!Tansport: initializing Driging rate to the propagation------------------------!
      if(TRANSPORT_CALC) then
         GammaMagnus=driving_rate
         GammaVerlet=GammaMagnus*0.1D0
         write(*,*) ' Driving Rate =', GammaMagnus
      end if
!------------------------------------------------------------------------------!

!##############################################################################!
! HERE STARTS THE TIME EVOLUTION
!##############################################################################!
            write(*,*) 'PROPAGATION'
            do 999 istep=1, ntdstep
!--------------------------------------!
              call g2g_timer_start('TD step')
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)&
                       .and.(.not.tdrestart)) then

                  t=(istep-1)*tdstep*0.1

              else
                 t=20*tdstep
                 t=t+(istep-200)*tdstep
              endif
              if (propagator.eq.1) then
                 t=(istep-1)*tdstep
              endif
              t=t*0.024190D0
              write(*,*) 'evolution time (fs)  =', t
!--------------------------------------!

            if ((propagator.eq.1).or. &
            (((propagator.eq.2).and.(istep.lt.lpfrg_steps)) &
            .and. (.not.tdrestart))) then
               call int3lu(E2)
               call g2g_solve_groups(0,Ex,0)
            endif

              write(*,*) '! step & energy', istep,E+Ex
              E1=0.0D0
! ELECTRIC FIELD CASE - Type=gaussian (ON)
         fxx=0.0D0
         fyy=0.0D0
         fzz=0.0D0
         if(.not.TRANSPORT_CALC) then
            if(istep.lt.pert_steps) then
               if (field) then
                 call dip(dipxyz)
                 if (exter) then
                   g=1.0D0
                   factor=2.54D0
                   fxx=fx*exp(-0.2D0*(real(istep-50))**2)
                   fyy=fy*exp(-0.2D0*(real(istep-50))**2)
                   fzz=fz*exp(-0.2D0*(real(istep-50))**2)
                   write(*,*) fxx,fyy,fzz
!
                 else
                   g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
                   Fx=dipxyz(1)/2.54D0
                   Fy=dipxyz(2)/2.54D0
                   Fz=dipxyz(3)/2.54D0
                   factor=(2.54D0*2.00D0)
!
                 endif
                 write(*,*) 'epsilon =', epsilon
                 call intfld(g,Fxx,Fyy,Fzz)
!TODO: shouldn't E1 use Fxx instead of Fx?? (idem y, z)
                 E1=-1.00D0*g*(Fx*dipxyz(1)+Fy*dipxyz(2)+Fz*dipxyz(3))&
                          /factor - 0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
              endif
            else
               field=.false.
            endif
         endif
!------------------------------------------------------------------------------!
! E1 includes solvent 1 electron contributions
            do k=1,MM
              E1=E1+RMM(k)*RMM(M11+k-1)
            enddo
!        write(*,*) '1 electron contribution',E1
!------------------------------------------------------------------------------!
! Here we obtain the fock matrix in the molecular orbital (MO) basis.
! where U matrix with eigenvectors of S , and s is vector with
! eigenvalues

            if ((propagator.eq.1).or. &
            (((propagator.eq.2).and.(istep.lt.lpfrg_steps)) &
            .and. (.not.tdrestart))) then
              call g2g_timer_start('fock')
              call spunpack('L',M,RMM(M5),fock)
#ifdef CUBLAS
!            call cumxtf(fock,devPtrX,fock,M)
!            call cumfx(fock,DevPtrX,fock,M)
!             call fock_ao_to_on(fock,devPtrX,fock,M)
              xnano2 = basechange_cublas(M,fock,devPtrX,'dir')
              fock=xnano2
#else
              xnano2=matmul(xtrans,fock)
              fock=matmul(xnano2,xmm)
!             call fock_ao_to_on(fock,x,fock,M)
#endif

!            do j=1,M
!              do k=1,j
!                 xnano(k,j)=RMM(M5+j+(M2-k)*(k-1)/2-1)
!              enddo
!              do k=j+1,M
!                 xnano(k,j)=RMM(M5+k+(M2-j)*(j-1)/2-1)
!              enddo
!            enddo
!            do i=1,M
!              do j=1,M
!                 X(i,M+j)=0.D0
!                 xnano2(i,j)=X(j,i)
!              enddo
!           enddo
!           call g2g_timer_start('actualiza rmm1')
!           do j=1,M
!              do k=1,M
!                 do i=1,M
!                 X(i,M+j)=X(i,M+j)+X(k,i)*xnano(k,j)
!                 enddo
!              enddo
!           enddo
!           call g2g_timer_stop('actualiza rmm1')
!           kk=0
!           do i=1,M
!              do k=1,M
!              xnano(k,i)=X(i,M+k)
!              enddo
!           enddo
!
!           do j=1,M
!              do i=j,M
!                 kk=kk+1
!                 RMM(M5+kk-1)=0.D0
!                 do k=1,M
!                    RMM(M5+kk-1)=RMM(M5+kk-1)+Xnano(k,i)*X(k,j)
!                  enddo
!               enddo
!            enddo
!
! Fock triangular matrix contained in RMM(M5,M5+1,M5+2,...,M5+MM) is copied to square matrix fock.
!            do j=1,M
!               do k=1,j
!                  fock(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
!               enddo
!               do k=j+1,M
!                  fock(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
!               enddo
!            enddo
              call sprepack('L',M,RMM(M5),fock)
              call g2g_timer_stop('fock')
           endif
! Now fock is stored in molecular orbital basis.
!
!  stores F1a and F1b for magnus propagation
            if((propagator.eq.2) .and. (.not.tdrestart)) then
               if(istep.eq.chkpntF1a) then
                  F1a=fock
               endif
               if(istep.eq.chkpntF1b) then
                  F1b=fock
               endif
            endif
!  stores F1a and F1b checkpoints to restart the dynamics
            if(writedens .and. propagator.eq.2) then
               kk=istep+5
               ii=istep+15
            if(mod (kk,500) == 0) then
               open(unit=7624,file='F1b.restart')
               rewind 7624
               do i=1,M
                  do j=1,M
                     write(7624,*) fock(i,j)
                  enddo
               enddo
               endif
               if(mod (ii,500) == 0) then
                 open(unit=7625,file='F1a.restart')
                 rewind 7625
                 do i=1,M
                    do j=1,M
                       write(7625,*) fock(i,j)
                    enddo
                 enddo
               endif
            endif
            E=E1+E2+En
            if (sol) then
                E=E+Es
            endif
!--------------------------------------------------------------------!
            if ((propagator.eq.1).or.(((propagator.eq.2)&
                 .and.(istep.lt.lpfrg_steps))&
                   .and.(.not.tdrestart))) then
                write(*,*) 'Verlet'
! In the first step of the propagation we extrapolate rho back in time
! using Verlet algorithm to calculate rhold.
! using matmul
!           if(istep.eq.1) then
!             rhold=rho+(dt_lpfrg*Im*(matmul(fock,rho)))
!             rhold=rhold-(dt_lpfrg*Im*(matmul(rho,fock)))
!           endif
! using commutator
              if(istep.eq.1) then
#ifdef CUBLAS
                call g2g_timer_start('cuconmut')
                rhold=commutator_cublas(fock,rho)
                rhold=rho+dt_lpfrg*(Im*rhold)
                call g2g_timer_stop('cuconmut')
#else
                call g2g_timer_start('conmutc')
                 rhold=commutator(fock,rho)
                 rhold=rho+dt_lpfrg*(Im*rhold)
                call g2g_timer_stop('conmutc')
#endif
              endif
!####################################################################!
! DENSITY MATRIX PROPAGATION USING VERLET ALGORITHM
! using matmul:
!           rhonew=rhold-(dt_lpfrg*Im*(matmul(fock,rho)))
!           rhonew=rhonew+(dt_lpfrg*Im*(matmul(rho,fock)))
!--------------------------------------c

!Transport: propagation in transport---------------------------------!

!####### Modificacion de transport
      if(TRANSPORT_CALC) then
         call g2g_timer_start('TRANSPORT - b Verlet -')
         if(istep.eq.1) then
            unit1 = 55555
            if ( Pop_Drive == 1 ) then
              open(unit=unit1,file='DriveMul')
            elseif ( Pop_Drive == 2 ) then
              open(unit=unit1,file='DriveLowd')
            endif
         endif
         if(istep.ge.3) then
           scratchgamma=GammaVerlet*exp(-0.0001D0*(dble(istep-1000))**2)
           call ELECTROSTAT(rho1,mapmat,overlap,rhofirst,scratchgamma,M)
           if ( Pop_Drive == 1 .or. Pop_Drive == 2 ) then
              if (mod(istep-1,save_charge_freq*10)==0) then
               rhoscratch=REAL(rho1)
               q(:) = 0.0D0
               call Drive_Population(Pop_Drive,ngroup,&
                    rhoscratch,overlap,group,sqsm,q,unit1)
              endif
           endif
         endif
!#####################################
#ifdef CUBLAS
            call g2g_timer_start('complex_rho_ao_to_on-cu')
            rho1=basechange_cublas(M,rho1,devPtrY,'dir')
           call g2g_timer_stop('complex_rho_ao_to_on-cu')
#endif
            call g2g_timer_stop('TRANSPORT - b Verlet -')

      endif  !end of transport propagation


! using commutator:
            call g2g_timer_start('commutator')
#ifdef CUBLAS
              rhonew=commutator_cublas(fock,rho)
              rhonew=rhold-dt_lpfrg*(Im*rhonew)
#else
              rhonew=commutator(fock,rho)
              rhonew=rhold-dt_lpfrg*(Im*rhonew)
#endif
            call  g2g_timer_stop('commutator')

!Transport: Add the driving term to the propagation------------------!

            if((istep.ge.3).and.(TRANSPORT_CALC)) then
               write(*,*) 'adding driving term to the density'
               rhonew=rhonew-rho1
            endif
!--------------------------------------------------------------------!

! Density update (rhold-->rho, rho-->rhonew)
              do i=1,M
                 do j=1,M
                    rhold(i,j)=rho(i,j)
                    rho(i,j)=rhonew(i,j)
                 enddo
              enddo
! END OF VERLET PROPAGATOR
!####################################################################!
              else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! DENSITY MATRIX PROPAGATION USING MAGNUS ALGORITHM
                 write(*,*) 'Magnus'

!################# MODIFICACION de TRANSPORT
         if(TRANSPORT_CALC) then
             call g2g_timer_start('TRANSPORT - b magnus -')
             if(istep.le.1000) then
             scratchgamma=GammaMagnus*exp(-0.0001D0*(dble(istep-1000))**2)
             else
             scratchgamma=GammaMagnus
             endif
             call ELECTROSTAT(rho1,mapmat,overlap,rhofirst,scratchgamma,M)
             if (Pop_Drive == 1 .or. Pop_Drive == 2) then
                if (mod(istep-1,save_charge_freq)==0) then
                 rhoscratch=REAL(rho1)
                 q(:) = 0.0D0
                 call Drive_Population(Pop_Drive,ngroup,&
                      rhoscratch,overlap,group,sqsm,q,unit1)
                endif
             endif

#ifdef CUBLAS
              call g2g_timer_start('complex_rho_ao_to_on-cu')
              rho1=basechange_cublas(M,rho1,devPtrY,'dir')
              call g2g_timer_stop('complex_rho_ao_to_on-cu')
#endif
              call g2g_timer_stop('TRANSPORT - b magnus -')

         endif ! end of transport propagation
!##############################################

#ifdef CUBLAS
                call g2g_timer_start('cupredictor')
                call cupredictor(F1a,F1b,fock,rho,devPtrX,factorial,&
                                 fxx,fyy,fzz,g,devPtrXc)
                call g2g_timer_stop('cupredictor')
                call g2g_timer_start('cumagnus')
                call cumagnusfac(fock,rho,rhonew,M,NBCH,dt_magnus,&
                                 factorial)
                call g2g_timer_stop('cumagnus')



!                rhold=rhonew
!                call g2g_timer_start('MAGNUS_MODIFIED')
!                call magnus_cublas(fock,rho,rhonew,M,NBCH,dt_magnus,
!                  factorial)
!                call g2g_timer_stop('MAGNUS_MODIFIED')
!                rhold=rhonew-rhold
!                write(22222222,*) rhold
!                stop 'hemos escrito rhold'
#else
                call g2g_timer_start('predictor')
                call predictor(F1a,F1b,fock,rho,factorial,fxx,fyy,fzz,g)
                call g2g_timer_stop('predictor')
                call g2g_timer_start('magnus')
                call magnus(fock,rho,rhonew,M,NBCH,dt_magnus,factorial)
                call g2g_timer_stop('magnus')
#endif

! Add the driving term to the propagation
         if(TRANSPORT_CALC .eqv. .true.) then
            write(*,*) 'adding driving term to the density'
            rhonew=rhonew-rho1
         endif

! density update and fock storage

                   F1a=F1b
                   F1b=fock
                   rho=rhonew
! END OF MAGNUS PROPAGATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
              endif
!####################################################################!
! Here we transform the density to the atomic orbital basis and take the real part of it. The imaginary part of the density
! can be descarted since for a basis set of purely real functions the fock matrix is real and symetric and depends only on
! the real part of the complex density matrix. (This wont be true in the case of hybrid functionals)
! with matmul:
#ifdef CUBLAS
             call g2g_timer_start('complex_rho_on_to_ao-cu')
!             call cumxp(rho,devPtrX,rho1,M)
!             call cumpxt(rho1,devPtrX,rho1,M)
             rho1=basechange_cublas(M,rho,devPtrXc,'inv')
             call g2g_timer_stop('complex_rho_on_to_ao-cu')
#else
             call g2g_timer_start('complex_rho_on_to_ao')
!---------------------------------------------------------------!
!FFR: I used to get an error that had to do with the size of X.
!     I'm not sure the problem persists, I'm temporarily leaving
!     here my implemented solution.
!
!        This should go in declarations:
!            real*8, allocatable :: xmat(:,:)
!
!        This should go here:
!             if (.not.allocated(xmat)) allocate(xmat(M,M))
!             do ii=1,M
!             do jj=1,M
!               xmat(ii,jj) = x(ii,jj)
!             end do
!             end do
!             rho1=matmul(xmat,rho)
!             rho1=matmul(rho1,xtrans)
!---------------------------------------------------------------!
             rho1=matmul(x(1:M,1:M),rho)
!             rho1=matmul(x,rho)
             rho1=matmul(rho1,xtrans)
             call g2g_timer_stop('complex_rho_on_to_ao')
#endif
!       rho1=REAL(rho1)
! with matmulnanoc:
! (NO LONGER AVAILABLE; USE BASECHANGE INSTEAD)
!          call matmulnanoc(rho,xtrans,rho1,M)
!          rho1=basechange(M,X,rho,Xtrans)
!          rho1 = REAL(rho1)
! The real part of the density matrix in the atomic orbital basis is copied in RMM(1,2,3,...,MM) to compute the corresponding fock matrix.
              do j=1,M
                  do k=j,M
                      if(j.eq.k) then
                        RMM(k+(M2-j)*(j-1)/2)=REAL(rho1(j,k))
                      else
                        RMM(k+(M2-j)*(j-1)/2)=(REAL(rho1(j,k)))*2
                      endif
                  enddo
              enddo
! Stores the density matrix each 500 steps to be able to restart the dynamics
              if(writedens) then
                 if(mod (istep,500) == 0) then
                     open(unit=5374,file='rho.restart')
                     rewind 5374
                     do j=1,M
                        do k=1,M
                           write(5374,*) rho1(j,k)
                        enddo
                     enddo
                     open(unit=7624,file='F1b.restart')
                     open(unit=7625,file='F1a.restart')
                     rewind 7624
                     rewind 7625
                     do i=1,M
                        do j=1,M
                           write(7624,*) F1b(i,j)
                           write(7625,*) F1a(i,j)
                        enddo
                     enddo
                  endif
! In the last step density matrix is stored
                  if (istep.eq.ntdstep) then
                    open(unit=44,file='rholast')
                    do j=1,M
                       do k=1,M
                          write(44,*) rho1(j,k)
                       enddo
                    enddo
                  endif
              endif

!Compute the trace of the density matrix for population analysis

         if(TRANSPORT_CALC) then
            traza=dcmplx(0.0D0,0.0D0)
            DO i=1,M
               traza=traza+rho(i,i)
            ENDDO
            write(*,*) 'TRAZA =', real(traza)
         endif

!###################################################################!
!# DIPOLE MOMENT CALCULATION
         if(istep.eq.1) then
            call write_dipole_td_header(tdstep, Fx, Fy, Fz, 134)
         endif
         if ((propagator.eq.2).and.(istep.lt.lpfrg_steps).and.(.not.tdrestart))&
         then
            if(mod ((istep-1),10) == 0) then
               call g2g_timer_start('DIPOLE_TD')
               call dip(dipxyz)
               call g2g_timer_stop('DIPOLE_TD')
               call write_dipole_td(dipxyz, t, 134)
            endif
         else
            call g2g_timer_start('DIPOLE_TD')
            call dip(dipxyz)
            call g2g_timer_stop('DIPOLE_TD')
            call write_dipole_td(dipxyz, t, 134)
         endif
! u in Debyes
!# END OF DIPOLE MOMENT CALCULATION
!------------------------------------------------------------------------------------
!-------------------------MULLIKEN CHARGES-----------------------------------------------!

!################ MODIFICACION TRANSPORT
            if (transport_calc) then
                if (istep .eq. 1) then
                    unit2 = 678
                    if ( Pop_Drive == 1 ) then
                       open(unit=unit2,file="MullikenGroup")
                    elseif ( Pop_Drive == 2 ) then
                       open(unit=unit2,file="LowdinGroup")
                    endif
                endif
                call g2g_timer_start("Mulliken Population")
                if (Pop_Drive == 1 .or. Pop_Drive == 2) then
                   do n=1,natom
                     q(n)=Iz(n)
                   enddo
                   rhoscratch=REAL(rho1)
                   if ((propagator.eq.2).and.(istep.lt.lpfrg_steps) &
                               .and. (.not.tdrestart)) then
                       if(mod((istep-1),save_charge_freq*10)==0) then
                          call Drive_Population(Pop_Drive,ngroup,&
                          rhoscratch,overlap,group,sqsm,q,unit2)
                       endif
                    elseif(mod((istep-1),save_charge_freq)==0) then
                   call Drive_Population(Pop_Drive,ngroup,&
                        rhoscratch,overlap,group,sqsm,q,unit2)
                    endif
                endif
            endif
!########################################

!!       if (nopt.ne.3) then
!       write(*,300) niter,DAMP,E
!       endif
!
!      write(*,*) 'Coulomb E',E2-Ex,Ex
               call g2g_timer_stop('TD step')
               write(*,*)
           if (istep .eq. 1000) then
             call g2g_timer_start('corrida 1000')
              tiempo1000=t
           elseif (istep .eq. 2000) then
               call g2g_timer_stop('corrida 1000')
               write(*,*) t-tiempo1000
           end if

!Transport: stops TD and indicate the generation of RHO0

         if((istep.ge.1).and.(generate_rho0)) then
            print*, "RHO0 GENERATED"
            exit
         endif

 999           continue
!
!##############################################################################!
! HERE FINISHES THE PROPAGATION
!##############################################################################!

 995   continue


         if (memo) then
            deallocate (kkind,kkinds)
            deallocate(cool,cools)
         endif
         if(propagator.eq.2) then
           deallocate (F1a,F1b)
         endif
         if (GRAD) then
            if(nopt.eq.0) then
              write(*,*)
              write(*,600)
              write(*,610)
              write(*,620) E1,E2-Ex,En
              if (sol) then
                 write(*,615)
                 write(*,625) Es
              endif
            endif
            write(*,*)
            write(*,450) E
         else
            E=E-Ex
         endif
! calculation of energy weighted density matrix

          kk=0
          do 307 j=1,M
             do 307 i=j,M
                kk=kk+1
                RMM(M15+kk-1)=0.D0
                if(i.eq.j) then
                    ff=2.D0
                else
                    ff=4.D0
                endif
                do 309 k=1,NCO
                   RMM(M15+kk-1)=RMM(M15+kk-1)-RMM(M13+k-1)&
       *ff*X(i,M2+k)*X(j,M2+k)
 309  continue
 307   continue

          if (nopt.eq.0) then
! calculates Mulliken poputations
                call int1(En)
                do n=1,natom
                   q(n)=Iz(n)
                enddo
                do i=1,M
                   do j=1,i-1
                      kk=i+(M2-j)*(j-1)/2
                      t0=RMM(kk)*RMM(M5+kk-1)/2.D0
                      q(Nuc(i))=q(Nuc(i))-t0
                   enddo
                   kk=i+(M2-i)*(i-1)/2
                   t0=RMM(kk)*RMM(M5+kk-1)
                   q(Nuc(i))=q(Nuc(i))-t0
                   do j=i+1,M
                      kk=j+(M2-i)*(i-1)/2
                      t0=RMM(kk)*RMM(M5+kk-1)/2.D0
                      q(Nuc(i))=q(Nuc(i))-t0
                   enddo
                 enddo
                 write(*,*) 'MULLIKEN POPULATION ANALYSIS'
                 write(*,770)
                 do n=1,natom
                    write(*,760) n,Iz(n),q(n)
                 enddo
                 write(*,*)
          endif
! ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION

!        if (icharge.eq.1) then
!          Q1=-(2*NCO+Nunp)
!         do n=1,natom
!          Q1=Q1+Iz(n)
!         enddo
!         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
!     >            c,a,RMM,map,Q1)
!        endif
!
! outputs final  MO ---------------------
!      do l=1,M
!      do n=1,NCO+3
!      do n=1,M
!        X(indexii(l),M+n)=X(l,M2+n)
!      enddo
!      enddo
!
!      do l=1,M
!        write(2,400) (X(l,M+n),n=1,NCO)
!      enddo
!--------------------------------------!
! Writes down MO coefficients and orbital energies
!       write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
!       do n=1,NCO
!         write(29,850) n,RMM(M13+n-1)
!         write(29,400) (X(l,M+n),l=1,M)
!       enddo
!       do n=NCO+1,M
!         write(29,851) n,RMM(M13+n-1)
!         write(29,400) (X(l,M+n),l=1,M)
!       enddo
!       close(29)
!--------------------------------------!
#ifdef CUBLAS
            call CUBLAS_FREE(devPtrX)
            call CUBLAS_FREE(devPtrXc)
            call CUBLAS_FREE(devPtrY)
            call CUBLAS_SHUTDOWN()
#endif


!---- DEBUGGINGS
!      write(*,*) 'Exc, integrated and calculated',Exc,Ex
!      write(*,*) 'Coulomb energy',E2-Ex

       close(134)
       call g2g_timer_stop('TD')
       DEALLOCATE(xnano,xnano2,fock,rhonew,rhold, &
                  rho,xmm,xtrans,Y,ytrans,rho1)
       DEALLOCATE(factorial)
!------------------------------------------------------------------------------!
 500  format('SCF TIME ',I6,' sec')
 450  format ('FINAL ENERGY = ',F19.12)
 400  format(4(E14.7E2,2x))
 300  format(I3,E14.6,2x,F14.7)
 600  format('  ENERGY CONTRIBUTIONS IN A.U.')
 610  format(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
 615  format(2x,'SOLVENT')
 620  format(F14.7,4x,F14.7,4x,F14.7)
 625  format(F14.7)
 760  format(I3,9x,I3,6x,F10.4)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
 850  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
 851  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7,&
         '(NON OCC.)')
 900  format(F15.9,2x,3(F15.9,2x),2x,F15.9)
 901  format(F15.9,2x,F15.9)
 777  format(4(F8.4,2x))
 776  format (3(F8.4,2x))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
