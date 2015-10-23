!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE TD()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c REAL TIME-TDDFT
c 
c Dario Estrin, 1992
c Nano, Dario, Uriel, Damian 2012
c
c  This subrutine takes the converged density matrix from an SCF calculation
c  and evolves it in time. In the input file the total number of propagation
c  steps is specified (nstep) as well as the time of each evolution step 
c  (tdstep). 
c  This implementation has two alternatives to evolve the density in time. The 
c  first one (propagator=1) is the Verlet algorithm that uses a convination of 
c  Liouville von Newmann expresion for the time derivative of the density matrix 
c  and a first order Taylor expansion of the density matrix. The second one 
c  (propagator=2) is the Magnus propagation scheme that uses Backer Campbell
c  Hausdorff (BCH) formula. For this reason when Magnus is used the number of 
c  total conmutators in the BCH espansion has to be specified (NBCH, default=10). 
c  A narrow gaussian type electric field can be introduced during the time 
c  evolution in order to excite all electronic frequencies with the same intensity.
c  Once this perturbation is turned on (Field=t, exter=t) each component of the
c  external electric field has to be specified in the input file (Fx,Fy,Fz).
c  In each step of the propagation the cartesian components of the sistems dipole
c  are stored in files x.dip, y.dip, z.dip.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c       USE latom
       USE garcha_mod
       use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC
       use mathsubs
#ifdef CUBLAS
        use cublasmath
#endif
       IMPLICIT REAL*8 (a-h,o-z)
       INTEGER :: istep
       REAL*8 :: t,E2
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: 
     >   xnano2,xmm,xtrans,ytrans,Y,fock,
     >   F1a,F1b,overlap,rhoscratch
       real*8, dimension (:,:), ALLOCATABLE :: elmu
#ifdef TD_SIMPLE
       COMPLEX*8 :: Im,Ix
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) ::
     >   rho,rhonew,rhold,xnano,rho1
#else
       COMPLEX*16 :: Im,Ix
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) ::
     >   rho,rhonew,rhold,xnano,rho1
#endif
       DIMENSION q(natom)
       REAL*8,dimension(:),ALLOCATABLE :: factorial
       INTEGER            :: LWORK,ii,jj
       REAL*8,ALLOCATABLE :: WORK(:)
!!------------------------------------!!
!! FFR ADD
       INTEGER ::
     >   pert_steps,lpfrg_steps,chkpntF1a,chkpntF1b
       REAL*8 ::
     >   dt_magnus,dt_lpfrg
        logical :: just_int3n,ematalloct
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
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX,CUBLAS_GET_MATRIX
#endif
!!   GROUP OF CHARGES
       LOGICAL             :: groupcharge
       INTEGER             :: ngroup
       INTEGER,ALLOCATABLE :: group(:)
       REAL*8,ALLOCATABLE  :: qgr(:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       call g2g_timer_start('TD')
       call g2g_timer_start('inicio')
       just_int3n = .false.
       ALLOCATE(factorial(NBCH))
!!------------------------------------!!
! Mulliken
       ipop=1
! Group of charges
       groupcharge=.false.
!!------------------------------------!!
#ifdef CUBLAS
       write(*,*) 'USING CUBLAS'
       stat=CUBLAS_INIT()
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
       Ndens=1
       E=0.0D0
       E1=0.0D0
       En=0.0D0
       E2=0.0D0
       idip=1
       ngeo=ngeo+1
       Im=(0.0D0,2.0D0)
       sq2=sqrt(2.D0)
       MM=M*(M+1)/2 
       MM2=M**2
       MMd=Md*(Md+1)/2
       Md2=2*Md
       M2=2*M
!
       ALLOCATE(xnano(M,M),xnano2(M,M),fock(M,M),rhonew(M,M),
     >   rhold(M,M),rho(M,M),xmm(M,M),xtrans(M,M),Y(M,M),ytrans(M,M),
     >   rho1(M,M))
!
      if(propagator.eq.2) allocate (F1a(M,M),F1b(M,M))
!--------------------------------------------------------------------!
      if (tdrestart) then
         inquire(file='rho.restart',exist=exists)
         if (.not.exists) then
             write(*,*) 'ERROR CANNOT FIND rho.restart'
             write(*,*) '(if you are not restarting a previous 
     > run set tdrestart= false)'
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
               write(*,*) '(if you are not restarting a 
     > previous run set tdrestart= false)'
               stop
            endif
            inquire(file='F1b.restart',exist=exists)
            if (.not.exists) then
               write(*,*) 'ERROR CANNOT FIND F1b.restart'
               write(*,*) '(if you are not restarting a
     > previous run set tdrestart= false)'
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
c first i
            M1=1
c now Fold
            M3=M1+MM
c now S, F also uses the same position after S was used
            M5=M3+MM
c now G
            M7=M5+MM
c now Gm
            M9=M7+MMd
c now H
            M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
            M13=M11+MM
c aux ( vector for ESSl)
            M15=M13+M
c Least squares
            M17=M15+MM
c vectors of MO
            M18=M17+MMd
c weights (in case of using option )
            M19=M18+M*NCO
c RAM storage of two-electron integrals (if MEMO=T)
            M20 = M19 + natom*50*Nang   
c
            Nel=2*NCO+Nunp
c Initializations/Defaults
       write(*,*) ' TD CALCULATION  '
!--------------------------------------!
           niter=0
           D1=1.D0
           D2=1.D0
!--------------------------------------!
           Qc=0.0D0
           do i=1,natom
             Qc=Qc+Iz(i)
           enddo
           Qc=Qc-Nel
           Qc2=Qc**2
!------------------------------------------------------------------------------!
! Two electron integral with neighbor list.
!
            do i=1,natom
              natomc(i)=0
!
              do j=1,natom
                d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >            (r(i,3)-r(j,3))**2
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
              nnps(nuc(ii))=ii
            enddo
            do ii=nshell(0)+nshell(1),nshell(0)+1,-1
              nnpp(nuc(ii))=ii
            enddo
            do ii=M,nshell(0)+nshell(1)+1,-1
              nnpd(nuc(ii))=ii
            enddo
!------------------------------------------------------------------!
c
c Create integration grid for XC here
c Assign points to groups (spheres/cubes)
c Assign significant functions to groups
c -Calculate point weights
c
      call g2g_timer_sum_start('Exchange-correlation grid setup')
      call g2g_reload_atom_positions(igrid2)
      call g2g_timer_sum_stop('Exchange-correlation grid setup')

      call aint_query_gpu_level(igpu)
      if (igpu.gt.1) call aint_new_step()

      if (predcoef.and.npas.gt.3) then

c        if (.not.OPEN) then
c          if(verbose) write(*,*) 'prediciendo densidad'
c          do i=1,MM
c            RMM(i)=(3*old1(i))-(3*old2(i))+(old3(i))
c          enddo
c         endif
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
            if(ipop.eq.1)then
              allocate(overlap(M,M),rhoscratch(M,M))
              call spunpack('L',M,RMM(M5),overlap)
            endif
!--------------------------------------!
c Diagonalization of S matrix, after this is not needed anymore
c s is in RMM(M13,M13+1,M13+2,...,M13+MM)
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
#ifdef CUBLAS
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=cmplx(X(i,j),0.0D0)
               ENDDO
            ENDDO
            stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
            stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrXc)
            stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrY)
            if (stat.NE.0) then
            write(*,*) "X and/or Y memory allocation failed"
            call CUBLAS_SHUTDOWN
            stop
            endif
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrXc,M)
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_real,x,M,devPtrX,M)
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=cmplx(Y(i,j),0.0D0)
               ENDDO
            ENDDO
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrY,M)
            if (stat.NE.0) then
            write(*,*) "X and/or Y setting failed"
            call CUBLAS_SHUTDOWN
            stop
            endif
            rho1=0
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
! with matmulnanoc
! (NO LONGER AVAILABLE; USE BASECHANGE INSTEAD)
!            call matmulnanoc(rho,Y,rho1,M)
!            rho=rho1
!            rho=basechange(M,Ytrans,rho,Y)
!--------------------------------------!
c Precalculate three-index (two in MO basis, one in density basis) matrix
c used in density fitting / Coulomb F element calculation here
c (t_i in Dunlap)
c
      call aint_query_gpu_level(igpu)
      if (igpu.gt.2) then
        call aint_coulomb_init()
      endif
      if (igpu.eq.5) MEMO = .false.
      !MEMO=.true.
      if (MEMO) then
         call g2g_timer_start('int3mem')
         call g2g_timer_sum_start('Coulomb precalc')
c Large elements of t_i put into double-precision cool here
c Size criteria based on size of pre-factor in Gaussian Product Theorem
c (applied to MO basis indices)
         call int3mem()
c Small elements of t_i put into single-precision cools here
c         call int3mems()
         call g2g_timer_stop('int3mem')
         call g2g_timer_sum_stop('Coulomb precalc')
      endif
#ifdef CUBLAS
            call CUBLAS_FREE(devPtrY)
#endif

            call g2g_timer_stop('inicio')
!##############################################################################!
! HERE STARTS THE TIME EVOLUTION
!##############################################################################!
            write(*,*) 'PROPAGATION'
            do 999 istep=1, ntdstep
!--------------------------------------!
              call g2g_timer_start('TD step')
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >      .and. (.not.tdrestart)) then
                 t=(istep-1)*tdstep*0.1
              else
                 t=20*tdstep
                 t=t+(istep-200)*tdstep
              endif
              if (propagator.eq.1) then
                 t=(istep-1)*tdstep
              endif
              t=t*0.02419
              write(*,*) 'evolution time (fs)  =', t
!--------------------------------------!
              call int3lu(E2)
              call g2g_solve_groups(0,Ex,0)
              write(*,*) '! step & energy', istep,E
              E1=0.0D0
c ELECTRIC FIELD CASE - Type=gaussian (ON)
            if(istep.lt.pert_steps) then
               if (field) then
                 call dip(ux,uy,uz)
                 if (exter) then
                   g=1.0D0
                   factor=2.54D0
                   fxx=fx*exp(-0.2*(real(istep-50))**2)
                   fyy=fy*exp(-0.2*(real(istep-50))**2)
                   fzz=fz*exp(-0.2*(real(istep-50))**2)
                   write(*,*) fxx,fyy,fzz
!
                 else
                   g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
                   Fx=ux/2.54D0
                   Fy=uy/2.54D0
                   Fz=uz/2.54D0
                   factor=(2.54D0*2.00D0)
!
                 endif
                 write(*,*) 'epsilon =', epsilon
                 call intfld(g,Fxx,Fyy,Fzz)
                 E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/factor -
     >        0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
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
c Fock triangular matrix contained in RMM(M5,M5+1,M5+2,...,M5+MM) is copied to square matrix fock.
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
c Now fock is stored in molecular orbital basis.
c
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
            if ((propagator.eq.1).or.
     >      (((propagator.eq.2).and.(istep.lt.lpfrg_steps))
     >      .and. (.not.tdrestart))) then
           write(*,*) 'Verlet'
c In the first step of the propagation we extrapolate rho back in time
c using Verlet algorithm to calculate rhold.
c using matmul 
!           if(istep.eq.1) then
!             rhold=rho+(dt_lpfrg*Im*(matmul(fock,rho)))
!             rhold=rhold-(dt_lpfrg*Im*(matmul(rho,fock)))
!           endif
c using commutator
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
c--------------------------------------c
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
c Density update (rhold-->rho, rho-->rhonew)
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
#ifdef CUBLAS
                call g2g_timer_start('cupredictor')
                call cupredictor(F1a,F1b,fock,rho,devPtrX,factorial,
     > fxx,fyy,fzz,g,devPtrXc) 
                call g2g_timer_stop('cupredictor')
                call g2g_timer_start('cumagnus')
                call cumagnusfac(fock,rho,rhonew,M,NBCH,dt_magnus,
     >factorial)
                call g2g_timer_stop('cumagnus')
!                rhold=rhonew
!                call g2g_timer_start('MAGNUS_MODIFIED')
!                call magnus_cublas(fock,rho,rhonew,M,NBCH,dt_magnus,
!     >factorial) 
!                call g2g_timer_stop('MAGNUS_MODIFIED')
!                rhold=rhonew-rhold
!                write(22222222,*) rhold
!                stop 'hemos escrito rhold'
#else
                call g2g_timer_start('predictor')
                call predictor(F1a,F1b,fock,rho,factorial,
     > fxx,fyy,fzz,g)
                call g2g_timer_stop('predictor')
                call g2g_timer_start('magnus')
                call magnus(fock,rho,rhonew,M,NBCH,dt_magnus,factorial)
                call g2g_timer_stop('magnus')
#endif
                 F1a=F1b
                 F1b=fock
                 rho=rhonew
! END OF MAGNUS PROPAGATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
              endif
!####################################################################!
c Here we transform the density to the atomic orbital basis and take the real part of it. The imaginary part of the density 
c can be descarted since for a basis set of purely real functions the fock matrix is real and symetric and depends only on 
c the real part of the complex density matrix. (This wont be true in the case of hybrid functionals)
c with matmul:
#ifdef CUBLAS
             call g2g_timer_start('complex_rho_on_to_ao-cu')
!             call cumxp(rho,devPtrX,rho1,M)
!             call cumpxt(rho1,devPtrX,rho1,M)
             rho1=basechange_cublas(M,rho,devPtrXc,'inv')
             call g2g_timer_stop('complex_rho_on_to_ao-cu')
#else
             call g2g_timer_start('complex_rho_on_to_ao')
             rho1=matmul(x,rho)
             rho1=matmul(rho1,xtrans)
             call g2g_timer_stop('complex_rho_on_to_ao')
#endif
!       rho1=REAL(rho1)
c with matmulnanoc:
c (NO LONGER AVAILABLE; USE BASECHANGE INSTEAD)
c          call matmulnanoc(rho,xtrans,rho1,M)
c          rho1=basechange(M,X,rho,Xtrans)
c          rho1 = REAL(rho1)
c The real part of the density matrix in the atomic orbital basis is copied in RMM(1,2,3,...,MM) to compute the corresponding fock matrix.
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
!###################################################################!
!# DIPOLE MOMENT CALCULATION
              if(istep.eq.1) then
                open(unit=134,file='x.dip')
                open(unit=135,file='y.dip')
                open(unit=136,file='z.dip')
                open(unit=13600,file='abs.dip')
        write(134,*) '#Time (fs) vs DIPOLE MOMENT, X COMPONENT (DEBYES)'
        write(135,*) '#Time (fs) vs DIPOLE MOMENT, Y COMPONENT (DEBYES)'
        write(136,*) '#Time (fs) vs DIPOLE MOMENT, Z COMPONENT (DEBYES)'
        write(13600,*) '#Time (fs) vs DIPOLE MOMENT (DEBYES)'
              endif
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >      .and. (.not.tdrestart)) then
                  if(mod ((istep-1),10) == 0) then
                     call g2g_timer_start('DIPOLE')
                     call dip(ux,uy,uz)
                     call g2g_timer_stop('DIPOLE')
                     write(134,901) t,ux
                     write(135,901) t,uy
                     write(136,901) t,uz
                  endif
              else
                  call g2g_timer_start('DIPOLE')
                  call dip(ux,uy,uz)
                  call g2g_timer_stop('DIPOLE')
                  write(134,901) t,ux
                  write(135,901) t,uy
                  write(136,901) t,uz
              endif
c u in Debyes
!# END OF DIPOLE MOMENT CALCULATION
c------------------------------------------------------------------------------------
!       if (nopt.ne.3) then
!       write(*,300) niter,DAMP,E
!       endif
c
c      write(*,*) 'Coulomb E',E2-Ex,Ex
               call g2g_timer_stop('TD step')
               write(*,*)
 999           continue
!
!##############################################################################!
! HERE FINISHES THE PROPAGATION
!##############################################################################!

 995   continue
c
c
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
c calculation of energy weighted density matrix
c
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
                   RMM(M15+kk-1)=RMM(M15+kk-1)-RMM(M13+k-1)
     >  *ff*X(i,M2+k)*X(j,M2+k)
 309  continue
 307   continue
c
          if (nopt.eq.0) then
c calculates Mulliken poputations
             if (ipop.eq.1) then
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
c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
c
c        if (icharge.eq.1) then
c          Q1=-(2*NCO+Nunp)
c         do n=1,natom
c          Q1=Q1+Iz(n)
c         enddo
c         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c     >            c,a,RMM,map,Q1)
c        endif
c
c outputs final  MO ---------------------
!      do l=1,M
c      do n=1,NCO+3
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
      endif
#ifdef CUBLAS
            call CUBLAS_FREE(devPtrX)
            call CUBLAS_FREE(devPtrXc)
            call CUBLAS_FREE(devPtrY)
            call CUBLAS_SHUTDOWN()
#endif
c
c
c---- DEBUGGINGS
c      write(*,*) 'Exc, integrated and calculated',Exc,Ex
c      write(*,*) 'Coulomb energy',E2-Ex
c
       call g2g_timer_stop('TD')
       deallocate(xnano,fock,rho)
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
 851  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7,
     >    '(NON OCC.)')
 900  format(F15.9,2x,3(F15.9,2x),2x,F15.9)
 901  format(F15.9,2x,F15.9)
 777  format(4(F8.4,2x))
 776  format (3(F8.4,2x))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
