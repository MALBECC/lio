!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c SCF subroutine
c DIRECT VERSION
c Calls all integrals generator subroutines : 1 el integrals,
c 2 el integrals, exchange fitting , so it gets S matrix, F matrix
c and P matrix in lower storage mode ( symmetric matrices)
c
c Dario Estrin, 1992
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      subroutine SCF(E,dipxyz)
      use ehrenfest
      use garcha_mod
      use mathsubs
      use general_module
      use faint_cpu77, only: int2, intsol, int3mem, int3lu
      use faint_cpu, only: int1
#ifdef CUBLAS
      use cublasmath
#endif
c      use qmmm_module, only : qmmm_struct, qmmm_nml
c
      implicit real*8 (a-h,o-z)
      integer:: l
       dimension q(natom),work(1000),IWORK(1000)
       REAL*8 , intent(inout)  :: dipxyz(3)
       real*8, dimension (:,:), allocatable::xnano,znano,scratch
       real*8, dimension (:,:), allocatable::scratch1
       real*8, dimension (:), allocatable :: rmm5,rmm15,rmm13,
     >   bcoef, suma
      real*8, dimension (:,:), allocatable :: fock,fockm,rho,!,FP_PF,
     >   FP_PFm,EMAT,Y,Ytrans,Xtrans,rho1,EMAT2,Xcpy
c
       integer ndiist
c       dimension d(natom,natom)
       logical  hagodiis,alloqueo, ematalloc
c       REAL*8 , intent(in)  :: qmcoords(3,natom)
c       REAL*8 , intent(in)  :: clcoords(4,nsolin)
        INTEGER :: ErrID,iii,jjj
        LOGICAL :: docholesky
        INTEGER            :: LWORK2
        REAL*8,ALLOCATABLE :: WORK2(:)
        INTEGER, ALLOCATABLE :: IWORK2(:),IPIV(:)
        logical :: just_int3n,ematalloct

!FFR!
!      NewForce
       real*8     :: Sinv(M,M),Fmtx(M,M),Bmat(M,M)
       complex*16 :: Pmtx(M,M)

!      vvterm
       logical             :: dovv
       real*8              :: weight
       integer,allocatable :: atom_group(:)
       integer,allocatable :: orb_group(:)
       integer,allocatable :: orb_selection(:)

       real*8,dimension(:,:),allocatable :: fockbias
       real*8,dimension(:,:),allocatable :: Xmat,Xtrp,Ymat,Ytrp
       real*8,dimension(:,:),allocatable :: sqsm
       real*8,dimension(:,:),allocatable :: Vmat
       real*8,dimension(:),allocatable   :: Dvec


#ifdef CUBLAS
        integer sizeof_real
        parameter(sizeof_real=8)
        integer stat
        integer*8 devPtrX, devPtrY
        external CUBLAS_INIT, CUBLAS_SET_MATRIX
        external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_FREE
        integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX
#endif
!--------------------------------------------------------------------!


#ifdef magma
       call magmaf_init()
#endif

      call g2g_timer_start('SCF')
      call g2g_timer_sum_start('SCF')
      call g2g_timer_sum_start('Initialize SCF')
c      just_int3n = .false.
      alloqueo = .true.
      ematalloc=.false.
      hagodiis=.false.
c     if(verbose)  write(6,*) 'ntatom',ntatom,nsol,natom

c------------------------------------------------------------------
c
c Pointers
c
c chequeo -----------
c
      Ndens=1
c---------------------
c       write(*,*) 'M=',M
      allocate (znano(M,M),xnano(M,M),scratch(M,M),scratch1(M,M),
     > fock(M,M))

      npas=npas+1
      E=0.0D0
      E1=0.0D0
      En=0.0D0
      E2=0.0D0
      Es=0.0D0
      Ens=0.0D0

      ngeo=ngeo+1

      sq2=sqrt(2.D0)
      MM=M*(M+1)/2
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M

c first P
      M1=1
c now Pnew
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
c
* RAM storage of two-electron integrals (if MEMO=T)
      M20 = M19 + natom*50*Nang

      if (cubegen_only.and.(cube_dens.or.cube_orb.or.cube_elec)) then
        if (.not.VCINP) then
          write(*,*) "cubegen_only CAN ONLY BE USED WITH VCINP"
          stop
        endif
        kk=0
        do k=1,NCO
          do i=1,M
            kk=kk+1
            Xnano(k,i) = MO_coef_at(kk)
          enddo
        enddo

        call g2g_timer_sum_start('cube gen')
        call cubegen(M15,Xnano)
        call g2g_timer_sum_stop('cube gen')

        deallocate (znano,xnano,scratch,scratch1)
        call g2g_timer_sum_stop('Initialize SCF')
        call g2g_timer_sum_stop('SCF')
        return
      endif
c
      Nel=2*NCO+Nunp
c
      allocate(rmm5(MM),rmm15(mm))
c
      good=1.00D0
      niter=0
      D1=1.D0
      D2=1.D0
      DAMP0=GOLD
      DAMP=DAMP0
c
      Qc=0.0D0
      do i=1,natom
       Qc=Qc+Iz(i)
      enddo
      Qc=Qc-Nel
      Qc2=Qc**2


! FFR: Variable Allocation
!--------------------------------------------------------------------!
       allocate(Xmat(M,M),Xtrp(M,M),Ymat(M,M),Ytrp(M,M))
       allocate(Vmat(M,M),Dvec(M))
       allocate(sqsm(M,M))
       allocate(fockbias(M,M))

       dovv=.false.
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


C----------------------------------------
c Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano

      do i=1,natom
        natomc(i)=0
        do j=1,natom
          d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
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

      ! get MM pointers in g2g
      ! call g2g_mm_init(nsol,r,pc)

c
c -Create integration grid for XC here
c -Assign points to groups (spheres/cubes)
c -Assign significant functions to groups
c -Calculate point weights
c
      call g2g_timer_sum_start('Exchange-correlation grid setup')
      call g2g_reload_atom_positions(igrid2)
      call g2g_timer_sum_stop('Exchange-correlation grid setup')

      call aint_query_gpu_level(igpu)
      if (igpu.gt.1) call aint_new_step()

      if (predcoef.and.npas.gt.3) then
        write(*,*) 'no devería estar aca!'

c        if (.not.OPEN) then
c          if(verbose) write(*,*) 'prediciendo densidad'
c          do i=1,MM
c            RMM(i)=(3*old1(i))-(3*old2(i))+(old3(i))
c          enddo
c         endif
       endif

c
c Calculate 1e part of F here (kinetic/nuc in int1, MM point charges
c in intsol)
c
      call g2g_timer_sum_start('1-e Fock')
      call g2g_timer_sum_start('Nuclear attraction')
      call int1(En, RMM, Smat, Nuc, a, c, d, r, Iz, ncont, NORM, natom, M, Md, nshell,ntatom)
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
c
c test ---------------------------------------------------------
      E1=0.D0
      do k=1,MM
        E1=E1+RMM(k)*RMM(M11+k-1)
      enddo
      call g2g_timer_sum_stop('1-e Fock')
c
c Diagonalization of S matrix, after this is not needed anymore
c S = YY^T ; X = (Y^-1)^T
c => (X^T)SX = 1
c
      docholesky=.true.
      call g2g_timer_start('cholesky')
      call g2g_timer_sum_start('Overlap decomposition')
      IF (docholesky) THEN
#ifdef magma
        ! ESTO SIGUE USANDO Smat EN RMM(M5)
        ! CAMBIARLO CUANDO SE SIGA PROBANDO MAGMA
        PRINT*,'DOING MAGMA-CHOLESKY'

        ALLOCATE(Y(M,M),Ytrans(M,M))
        DO iii=1,M;DO jjj=1,M
          Y(iii,jjj)=0
          IF (jjj.LE.iii) THEN
            iiindex=iii+(2*M-jjj)*(jjj-1)/2
            Y(iii,jjj)=RMM(M5+iiindex-1)
          ENDIF
        ENDDO;ENDDO

        CALL MAGMAF_DPOTRF('L',M,Y,M,ErrID)

        Ytrans= transpose(Y)
        ALLOCATE(Xtrans(M,M))
        Xtrans=Y
        CALL MAGMAF_DTRTRI('L','N',M,Xtrans,M,ErrID)
        if(ErrID.ne.0) STOP ('Error in cholesky decomp.')
        xnano= transpose(Xtrans)
        do i=1,M;doj=1,M
          X(i,j)=Xnano(i,j)
        enddo;enddo
        PRINT*,'CHOLESKY MAGMA'

#else

! FFR: Cholesky Decomposition of Overlap
!--------------------------------------------------------------------!
! I am keeping Y,Ytrans and Xtrans but they should be replaced
! by the much nicer Ymat,Ytrp,Xtrp (and X by Xmat). The outputs
! Dvec and Vmat don't have the same meaning as in the diagona-
! lization (they are not eigenvalues or eigenvectors. No S1/2
! matrix can be obtained.
!
! Magma option should be introduced INSIDE of che call
!
         call sdcmp_cholesky(Smat,Dvec,Vmat,Ymat,Xtrp,Ytrp,Xmat)

         allocate (Y(M,M),Ytrans(M,M),Xtrans(M,M))
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

#endif
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

         allocate (Y(M,M),Ytrans(M,M),Xtrans(M,M))
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
!! CUBLAS ---------------------------------------------------------------------!
#ifdef CUBLAS
            stat=CUBLAS_INIT()
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
      call g2g_timer_start('initial guess')
      call g2g_timer_sum_stop('Overlap decomposition')

c
c CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
c FCe = SCe; (X^T)SX = 1
c F' = (X^T)FX
c => (X^-1*C)^-1 * F' * (X^-1*C) = e
c

c Calculate F' in RMM(M5)
      if((.not.VCINP).and.primera) then
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
c
c F' diagonalization now
c xnano will contain (X^-1)*C
c
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M13+i-1)=0.D0
        enddo
c
c ESSL OPTION
        do i=1,MM
          rmm5(i)=RMM(M5+i-1)
        enddo
        rmm15=0
        xnano=0
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c LAPACK OPTION -----------------------------------------
#ifdef pack
c
        call dspev('V','L',M,RMM5,RMM(M13),Xnano,M,RMM15,info)
#endif
        do i =1,M
          do j=1,M
            X(i,M+j)=xnano(i,j)
          enddo
        enddo
c-----------------------------------------------------------
c Recover C from (X^-1)*C
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
c
c Density Matrix
c
        kk=0
c
        do k=1,NCO
          do i=1,M
            kk=kk+1
            MO_coef_at(kk)=X(i,M2+k)
          enddo
        enddo
c
        kk=0
        do j=1,M
          do i=j,M
            kk=kk+1
            RMM(kk)=0.D0
c
c one factor of 2 for alpha+beta
            if(i.eq.j) then
             ff=2.D0
c another factor of 2 for direct triangular sum (j>i) w/ real basis
            else
             ff=4.D0
            endif
c
            do k=1,NCO
              RMM(kk)=RMM(kk)+ff*X(i,M2+k)*X(j,M2+k)
            enddo
          enddo
        enddo
c
        call g2g_timer_sum_stop('initial guess')
      endif

c End of Starting guess (No MO , AO known)-------------------------------
c
      if ((timedep.eq.1).and.(tdrestart)) then
        call g2g_timer_sum_start('TD')
        call TD()
        call g2g_timer_sum_stop('TD')
        return
      endif
c
c Precalculate two-index (density basis) "G" matrix used in density fitting
c here (S_ij in Dunlap, et al JCP 71(8) 1979) into RMM(M7)
c Also, pre-calculate G^-1 if G is not ill-conditioned into RMM(M9)
c
      call g2g_timer_sum_start('Coulomb G matrix')
      call int2()
      call g2g_timer_sum_stop('Coulomb G matrix')
c
**
c
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
         call g2g_timer_stop('int3mem')
         call g2g_timer_sum_stop('Coulomb precalc')
      endif
****
c---------------------------------------------------------------------
c Now, damping is performed on the density matrix
c The first 4 iterations ( it may be changed, if necessary)
c when the density is evaluated on the grid, the density
c matrix is used ( slower), after that it is calculated
c using the vectors . Since the vectors are not damped,
c only at the end of the SCF, the density matrix and the
c vectors are 'coherent'
c---------------------------------------------------------------
c LEVEL SHIFT CASE, contruction of initial vectors ------------------
c
      if (SHFT) then
c
        write(*,*) 'Level SHIFT is not suported'
        stop
      endif
c
      if (DIIS.and.alloqueo) then
        alloqueo=.false.
c       write(*,*) 'eme=', M
       allocate(rho1(M,M),rho(M,M),fockm(MM,ndiis),
     >  FP_PFm(MM,ndiis),EMAT(ndiis+1,ndiis+1),bcoef(ndiis+1)
     >  ,suma(MM))
      endif
      call g2g_timer_sum_stop('Initialize SCF')
      do i=1,MM
         if(RMM(i).ne.RMM(i)) stop 'NAN en RHO'
         if(RMM(M5+i).ne.RMM(M5+i)) stop 'NAN en fock'
      enddo
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c      write(*,*) 'empiezo el loop',NMAX
c-------------------------------------------------------------------
c-------------------------------------------------------------------
      do 999 while (good.ge.told.and.niter.le.NMAX)

        call g2g_timer_start('Total iter')
        call g2g_timer_sum_start('Iteration')
        call g2g_timer_sum_start('Fock integrals')
        niter=niter+1
        if(niter.le.ndiis) then
          ndiist=niter
        else
          ndiist=ndiis
        endif

c      if (MEMO) then
c
c Fit density basis to current MO coeff and calculate Coulomb F elements
c
            call g2g_timer_sum_start('Coulomb fit + Fock')
            call int3lu(E2)
            call g2g_timer_sum_pause('Coulomb fit + Fock')
c
c XC integration / Fock elements
c
            call g2g_timer_sum_start('Exchange-correlation Fock')
            call g2g_solve_groups(0,Ex,0)
            call g2g_timer_sum_pause('Exchange-correlation Fock')

c-------------------------------------------------------
        E1=0.0D0
c
c REACTION FIELD CASE --------------------------------------------
c
        call g2g_timer_start('actualiza rmm')
c----------------------------------------------------------------
c E1 includes solvent 1 electron contributions
        do k=1,MM
          E1=E1+RMM(k)*RMM(M11+k-1)
        enddo
        call g2g_timer_sum_pause('Fock integrals')
        call g2g_timer_sum_start('SCF acceleration')
c
c
c now, we know S matrix, and F matrix, and E for a given P
c 1) diagonalize S, get X=U s^(-1/2)
c 2) get U F Ut
!c 3) diagonalize F
c 4) get vec ( coeff) ---->  P new
c 5) iterate
c call diagonalization routine for S , get after U s^(-1/2)
c where U matrix with eigenvectors of S , and s is vector with
c eigenvalues
c
c here in RMM(M5) it is stored the new Fock matrix
c test damping on Fock matrix
c
c      DAMP=gold
        if (niter.eq.1) then
          DAMP=0.0D0
        endif

c-----------------------------------------------------------------------------------------
c If DIIS is turned on, update fockm with the current transformed F' (into ON
c basis) and update FP_PFm with the current transformed [F',P']
c-----------------------------------------------------------------------------------------
        if (DIIS) then
          call g2g_timer_sum_start('DIIS')
          call g2g_timer_sum_start('DIIS prep')
c-----------------------------------------------------------------------------------------
c Expand F into square form
c (for better memory access patterns in coming multiplications)
c-----------------------------------------------------------------------------------------
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
         Fmtx=fock ! TODO: Solo importa quedarse con la ultima

c-----------------------------------------------------------------------------------------
c Expand density matrix into full square form (before, density matrix was set up for triangular sums
c (sum j>=i) so off-diagonal elements need to be divided by 2 to get square-form numbers)
c (for better memory access patterns in coming multiplications)
c-----------------------------------------------------------------------------------------
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
#ifdef CUBLAS
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
c-----------------------------------------------------------------------------------------
c now, scratch = A = F' * P'; scratch1 = A^T
c [F',P'] = A - A^T
c-----------------------------------------------------------------------------------------
          ! update FP_PFm with [F',P']
          do j=ndiis-(ndiist-1),ndiis-1
            do i=1,MM
              FP_PFm(i,j)=FP_PFm(i,j+1)
            enddo
          enddo
#ifdef CUBLAS
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
c
c-------------Decidiendo cual critero de convergencia usar-----------
c-----------------------------------------------------------------------------------------
c iF DIIS=T
c Do simple damping 2nd iteration; DIIS afterwards
c-----------------------------------------------------------------------------------------
c IF DIIS=F
c Always do damping (after first iteration)
c-----------------------------------------------------------------------------------------
        if (niter.gt.2.and.(DIIS)) then
          hagodiis=.true.
        endif

c-----------------------------------------------------------------------------------------
c If we are not doing diis this iteration, apply damping to F, save this
c F in RMM(M3) for next iteration's damping and put F' = X^T * F * X in RMM(M5)
c-----------------------------------------------------------------------------------------
        if(.not.hagodiis) then
          call g2g_timer_start('Fock damping')
          if(niter.ge.2) then
            do k=1,MM
              kk=M5+k-1
              kk2=M3+k-1
              RMM(kk)=(RMM(kk)+DAMP*RMM(kk2))/(1.D0+DAMP)
            enddo
          endif
c the newly constructed damped matrix is stored, for next iteration
c in RMM(M3)
c
         do k=1,MM
            kk=M5+k-1
            kk2=M3+k-1
            RMM(kk2)=RMM(kk)
          enddo
c
! xnano=X^T
!          do i=1,M
!            ! X is upper triangular
!            do j=1,i
!              xnano(i,j)=X(j,i)
!            enddo
!          enddo

! RMM(M5) gets F' = X^T * F * X
!          do j=1,M
!            do i=1,M
!              X(i,M+j)=0.D0
!            enddo
!            do k=1,j
!              ! xnano is lower triangular
!              do i=k,M
!                X(i,M+j)=X(i,M+j)+Xnano(i,k)*RMM(M5+j+(M2-k)*(k-1)/2-1)
!              enddo
!            enddo
!c
!            do k=j+1,M
!              ! xnano is lower triangular
!              do i=k,M
!                X(i,M+j)=X(i,M+j)+Xnano(i,k)*RMM(M5+k+(M2-j)*(j-1)/2-1)
!              enddo
!            enddo
!
!          enddo
!c
!          kk=0
!          do i=1,M
!            do k=1,M
!              xnano(k,i)=X(i,M+k)
!            enddo
!          enddo
!!
!          do j=1,M
!            do i=j,M
!              kk=kk+1
!              RMM(M5+kk-1)=0.D0
!              ! X is upper triangular
!              do k=1,j
!                RMM(M5+kk-1)=RMM(M5+kk-1)+Xnano(k,i)*X(k,j)
!              enddo
!            enddo
!          enddo
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
         Fmtx=fock ! TODO: Solo importa quedarse con la ultima
!         write(601,*) Fmtx


#ifdef CUBLAS
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
c
c now F contains transformed F
c diagonalization now
c
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M13+i-1)=0.D0
        enddo

c---- LEVEL SHIFT
c
        if (SHFT) then
c       shi=shi*0.99
c adition of level shifts
c constant to diagonal (virtual) elements
          do i=NCO+1,M
            ii=i+(i-1)*(M2-i)/2
            RMM(M5+ii-1)=RMM(M5+ii-1)+shi
          enddo
        endif

        call g2g_timer_stop('actualiza rmm')

c----------Si hagodiis(ver mas arriba) es true entonces sigo-----------------------
c        write(*,*) 'good < dgtrig DIIS!!! PARA LA SIGUIENTE ITERACION'
        call g2g_timer_start('diis')
c--------Pasar columnas de FP_PFm a matrices y multiplicarlas y escribir EMAT-------

        if(DIIS) then
          call g2g_timer_sum_start('DIIS')
          call g2g_timer_sum_start('DIIS Fock update')
          deallocate(EMAT)
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
c-------Escribimos en xnano y znano dos conmutadores de distintas iteraciones------
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

            EMAT(ndiist,kk)=0.0d0
            if(kk.ne.ndiist) EMAT(kk,ndiist)=0.0d0
            do l=1,M
              EMAT(ndiist,kk)=EMAT(ndiist,kk)+rho1(l,l)
              if (kk.ne.ndiist) then
                EMAT(kk,ndiist)=EMAT(ndiist,kk)
              endif
            enddo
          enddo

          do i=1,ndiist
            EMAT(i,ndiist+1)= -1.0d0
            EMAT(ndiist+1,i)= -1.0d0
          enddo
          EMAT(ndiist+1, ndiist+1)= 0.0d0

          allocate(EMAT2(ndiist+1,ndiist+1))
          EMAT2=EMAT
c        ematalloct=.true.
c********************************************************************
c   THE MATRIX EMAT SHOULD HAVE FORM
c
c      |<E(1)*E(1)>  <E(1)*E(2)> ...   -1.0|
c      |<E(2)*E(1)>  <E(2)*E(2)> ...   -1.0|
c      |<E(3)*E(1)>  <E(3)*E(2)> ...   -1.0|
c      |<E(4)*E(1)>  <E(4)*E(2)> ...   -1.0|
c      |     .            .      ...     . |
c      |   -1.0         -1.0     ...    0. |
c
c   WHERE <E(I)*E(J)> IS THE SCALAR PRODUCT OF [F*P] FOR ITERATION I
c   TIMES [F*P] FOR ITERATION J.
c
c********************************************************************
c-----Pasamos a resolver con DGELS el problema EMAT*ci=bcoef-------
          if (hagodiis) then
            do i=1,ndiist
              bcoef(i)=0
            enddo
            bcoef(ndiist+1)=-1

c----------Cálculo de parámetro optimo para DGELS-------------------
*
*
        LWORK = -1
      CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT,
     >ndiist+1, bcoef, ndiist+1, WORK, LWORK, INFO )
        LWORK = MIN( 1000, INT( WORK( 1 ) ) )


c-----Resuelve la ecuación A*X = B. (EMAT*ci=bcoef). La solución la escribe en bcoef------

      CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT,
     > ndiist+1, bcoef, ndiist+1, WORK, LWORK, INFO )

c--------Construccion de la "nueva" matriz de fock como cl de las anteriores--------------
c--------Eventualmente se puede probar con la matriz densidad-----------------------------
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

c
c F' diagonalization now
c X(1,M) will contain (X^-1)*C
c

       call g2g_timer_start('dspev')
       call g2g_timer_sum_pause('SCF acceleration')
       call g2g_timer_sum_start('diagonalization')
c ESSL OPTION ---------------------------------------------------
#ifdef essl
       call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
c#ifdef magma
c-------nano tratando de usar magma
      if(.not.allocated(fock)) allocate (fock(M,M))
      fock=0
      do j=1,M
        do k=1,j
         i=j+(M2-k)*(k-1)/2
         fock(j,k)=RMM(M5+i-1)
        enddo
      enddo
c---------------------
       LWORK=-1
#ifdef magma
      call magmaf_dsyevd('V','L',M,fock,M,RMM(M13),WORK,LWORK
     > ,IWORK,LWORK,info)
#else
      call dsyevd('V','L',M,fock,M,RMM(M13),WORK,LWORK
     > ,IWORK,LWORK,info)
#endif

       LWORK=work(1)
      LIWORK=IWORK(1)

      if(allocated(WORK2)) deallocate (WORK2)
      if(allocated(IWORK2)) deallocate (IWORK2)

       allocate (WORK2(LWORK),IWORK2(LIWORK))


#ifdef magma
      call magmaf_dsyevd('V','L',M,fock,M,RMM(M13),WORK2,LWORK
     > ,IWORK2,LIWORK,info)
#else
      call dsyevd('V','L',M,fock,M,RMM(M13),WORK2,LWORK
     > ,IWORK2,LIWORK,info)
#endif
c#else
c       call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),
c     > M,RMM(M15),info)
c#endif
#endif
       call g2g_timer_stop('dspev')
       call g2g_timer_sum_pause('diagonalization')
c       do ik=1,M
c         do jk=1,M
c         write(45,*) X(ik,M+jk),fock(ik,jk)
c
c
c         enddo
c
c       enddo
       call g2g_timer_start('coeff')
       call g2g_timer_sum_start('MO coefficients')

c-----------------------------------------------------------
c
c diagonalization now
c-----------------------------------------------------------
c Recover C from (X^-1)*C; put into xnano
c-----------------------------------------------------------
#ifdef CUBLAS
        call cumxp_r(fock,devPtrX,xnano,M)
        do i=1,M
          do j=1,M
             X(i,M2+j)=xnano(i,j)
          enddo
       enddo
#else
c new coefficients
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
!       do i=1,M
!         do k=1,M
!           xnano(i,k)=X(k,i)
!         enddo
!       enddo
!       do i=1,M
!        do j=1,M
!            X(i,M2+j)=0.D0
!            ! xnano is lower triangular
!            do k=i,M
!              X(i,M2+j)=X(i,M2+j)+xnano(k,i)*fock(k,j)
!            enddo
!          enddo
!      enddo
#endif
      call g2g_timer_stop('coeff')
      call g2g_timer_start('otras cosas')
      call g2g_timer_sum_pause('MO coefficients')
      call g2g_timer_sum_start('new density')
c
c --- For the first iteration, damping on density matrix
c Important for the case of strating guess of AO
c
      kk=0
      do k=1,NCO
        do i=1,M
          kk=kk+1
          MO_coef_at(kk)=X(i,M2+k)
          xnano(k,i)  = X(i,M2+k)
        enddo
      enddo
c
c Construction of new density matrix and comparison with old one
      kk=0
      good=0.
      call g2g_timer_start('dens_GPU')
      call density(M,NCO,X,xnano)
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
c Level Shifting
         do i=1,M
           do j=1,M
             X(i,j)=X(i,M2+j)
           enddo
         enddo
       endif
c
       if (SHFT) then
         DAMP=0.0D0
       endif
*
c--- Damping factor update -
       DAMP=DAMP0
       IDAMP=0
       if (IDAMP.EQ.1) then
         DAMP=DAMP0
         if (abs(D1).lt.1.D-5) then
           factor=dmax1(0.90D0,abs(D1/D2))
           factor=dmin1(factor,1.1D0)
           DAMP=DAMP0*factor
         endif
c
         E=E1+E2+En
         E=E+Es
c
         D2=D1
         D1=(E-E0)
c
         E0=E
         DAMP0=DAMP
       endif
c
        E=E1+E2+En
        E=E+Es
c
c
        call g2g_timer_stop('otras cosas')
        call g2g_timer_sum_pause('new density')

        if(verbose) write(6,*) 'iter',niter,'QM Energy=',E+Ex
c
        call g2g_timer_stop('Total iter')
        call g2g_timer_sum_pause('Iteration')
 999  continue

c-------------------------------------------------------------------
c
c-------------------------------------------------------------------
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

c      old3=old2

c      old2=old1
c        write(*,*) 'good final',good

c      do i=1,MM
c        old1(i)=RMM(i)
c      enddo

      if(noconverge.gt.4) then
        write(6,*)  'stop fon not convergion 4 times'
        stop
      endif

c
!    CH - Why call intsol again here? with the .false. parameter,
!    E1s is not recalculated, which would be the only reason to do
!    this again; Ens isn't changed from before...
c -- SOLVENT CASE --------------------------------------
c      if (sol) then
c      call g2g_timer_sum_start('intsol 2')
c      if(nsol.gt.0) then
c        call intsol(E1s,Ens,.false.)
c        write(*,*) 'cosillas',E1s,Ens
c        call g2g_timer_sum_stop('intsol 2')
c      endif
c      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
      Es=Es+E1s+Ens
c     endif
c--------------------------------------------------------------
c  ????
      if (MOD(npas,energy_freq).eq.0) then
      if (GRAD) then
c         if (sol) then
c         endif
c       call g2g_timer_sum_start('exchnum')
        call g2g_timer_sum_start('Exchange-correlation energy')
#ifdef G2G
#ifdef ULTIMA_CPU
        call exchnum(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >              M18,NCO,Exc,nopt)
#else
      ! Resolve with last density to get XC energy
        call g2g_new_grid(igrid)
        call g2g_solve_groups(1, Exc, 0)
c       write(*,*) 'g2g-Exc',Exc
#endif
#else
#ifdef ULTIMA_G2G
        call g2g_new_grid(igrid)
        call g2g_solve_groups(1, Exc, 0)
#else
#endif
#endif
        call g2g_timer_sum_stop('Exchange-correlation energy')
        ! -------------------------------------------------
        ! Total SCF energy =
        ! E1 - kinetic+nuclear attraction+QM/MM interaction
        ! E2 - Coulomb
        ! En - nuclear-nuclear repulsion
        ! Ens - MM point charge-nuclear interaction
        ! Exc - exchange-correlation
        ! -------------------------------------------------
        E=E1+E2+En+Ens+Exc
        if (npas.eq.1) npasw = 0
        if (npas.gt.npasw) then
          write(6,*)
          write(6,600)
          write(6,610)
          write(6,620) E1,E2-Ex,En
c         if (sol) then
c          write(6,615)
c          write(6,625) Es
c          write(6,*) 'E SCF = ', E , Exc, Ens
          npasw=npas+10
        endif
c--------------------------------------------------------------
      else
        E=E-Ex
      endif
      endif
      write(888,*) 'SCF energy = ',E

c calculation of energy weighted density matrix
c
      call g2g_timer_sum_start('energy-weighted density')
      kk=0
      do j=1,M
        do i=j,M
          kk=kk+1
          RMM(M15+kk-1)=0.D0
c
          if(i.eq.j) then
            ff=2.D0
          else
            ff=4.D0
          endif
          do k=1,NCO
            RMM(M15+kk-1)=
     >      RMM(M15+kk-1)-RMM(M13+k-1)*ff*X(i,M2+k)*X(j,M2+k)
          enddo
        enddo
      enddo
      call g2g_timer_sum_stop('energy-weighted density')

      if (MOD(npas,energy_freq).eq.0) then
c
c      if (nopt.eq.1) then
c
c PROPERTIES CALCULATION
c calculates dipole moment
c
      if (idip.eq.1) then
       call g2g_timer_sum_start('dipole')
        call dip(ux,uy,uz)
        u=sqrt(ux**2+uy**2+uz**2)
        dipxyz(1)=ux
        dipxyz(2)=uy
        dipxyz(3)=uz
c
c      write(*,*)
c      write(*,*) 'DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
c       write(69,900) ux,uy,uz,u
c      write(*,*)
c u in Debyes
       call g2g_timer_sum_stop('dipole')
      endif
c



       call g2g_timer_sum_start('Mulliken')
! MULLIKEN POPULATION ANALYSIS (FFR - Simplified)
!--------------------------------------------------------------------!
       call int1(En, RMM, Smat, Nuc, a, c, d, r, Iz, ncont, NORM, natom, M, Md, nshell,ntatom)
       call spunpack('L',M,RMM(M5),Smat)
       call spunpack('L',M,RMM(M1),RealRho)
       call fixrho(M,RealRho)
       call mulliken_calc(natom,M,RealRho,Smat,Nuc,Iz,q)
       call mulliken_write(85,natom,Iz,q)
       RhoCero=DCMPLX(RealRho)
       if (first_step) then
         RhoSave=DCMPLX(RealRho)
       endif



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! TESTIN THE FORCE
       if (allocated(kkind))  deallocate(kkind)
       if (allocated(kkinds)) deallocate(kkinds)
       if (allocated(cool))   deallocate(cool)
       if (allocated(cools))  deallocate(cools)


!       print*,'-------------------------------------IGNORE FROM HERE'
       call aint_query_gpu_level(igpu)
       if (igpu.le.1) then
         call intsol(Energy1,Energy2,.true.)
       else
         call aint_qmmm_fock(Energy1,Energy2)
       endif
       call int2()
       if (igpu.gt.2) then
         call aint_coulomb_init()
       endif
       if (igpu.eq.5) MEMO = .false.
       if (MEMO) then
          call g2g_timer_start('int3mem')
          call g2g_timer_sum_start('Coulomb precalc')
          call int3mem()
          call g2g_timer_stop('int3mem')
          call g2g_timer_sum_stop('Coulomb precalc')
       endif
       call int3lu(E2)
       call g2g_solve_groups(0,Ex,0)
!       print*,'----------------------------------------------TO HERE'
       call spunpack('L',M,RMM(M5),Fmtx)   !TODO: sacar si fock esta bien
!       write(602,*) Fmtx
       Sinv=matmul(Xmat,Xtrans)

       do iii=1,M
       do jjj=1,M
         Pmtx(iii,jjj)=CMPLX(RealRho(iii,jjj),0.0d0)
       enddo
       enddo
! NECESITO TESTFORCE PORQUE SETEA EL BASIS_DATA!
!       call testforce(Sinv,Fmtx,Pmtx)

       if (.not.allocated(qm_forces_ds)) then
          allocate(qm_forces_ds(3,natom))
          qm_forces_ds=0.0d0
       endif


       print*,natom,3
       print*,size(nucvel,2),size(nucvel,1)
       print*,size(qm_forces_ds,2),size(qm_forces_ds,1)
       print*,'ALMOST';STOP!FFR-DELETE
       do ii=1,natom
       do kk=1,3
          nucvel(kk,ii)=nucvel(kk,ii)-
     >    dt*qm_forces_ds(kk,ii)/atom_mass(ii)
       enddo
       enddo
       print*,'ALMOST';STOP!FFR-DELETE

       call calc_forceDS(natom,M,nucpos,nucvel,Pmtx,Fmtx,Sinv,Bmat
     >                  ,qm_forces_ds)


!       print*,'----------------------------------------------DONE'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


! NOTE: If 'mulliken_calc' is renamed as 'mulliken', the code will
! malfunction. I DON'T KNOW WHY.
!--------------------------------------------------------------------!
       call g2g_timer_sum_stop('Mulliken')
!       do kk=1,natom
!         q(kk)=real(Iz(kk))
!       enddo
!       call lowdinpop(M,natom,RealRho,sqsm,Nuc,q)
!       call mulliken_write(85,natom,Iz,q)
       endif

c
c        endif
c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
c
c        if (icharge.eq.1) then
c          Q1=-(2*NCO+Nunp)
c         do n=1,natom
c          Q1=Q1+Iz(n)
c         enddo
c          stop
c         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c     >            c,a,RMM,map,Q1)
c        endif
c
c--------------------------------------------------------------
c outputs final  MO ---------------------

      if (MOD(npas,restart_freq).eq.0) then
      call g2g_timer_sum_start('restart write')
      rewind 88
      do l=1,M
        do n=1,M
          X(indexii(l),M+n)=X(l,M2+n)
        enddo
      enddo
c
      do l=1,M
        write(88,400) (X(l,M+n),n=1,NCO)
      enddo
      call g2g_timer_sum_stop('restart write')
      endif
c-------------------------------------------------
c writes down MO coefficients and orbital energies
      if(1.gt.2) then
        write(29,*) 'ORBITAL COEFFICIENTS eh AND ENERGIES, CLOSED SHELL'
        do n=1,NCO
          write(29,850) n,RMM(M13+n-1)
          write(29,400) (X(l,M+n),l=1,M)
        enddo
        do n=NCO+1,M
          write(29,851) n,RMM(M13+n-1)
          write(29,400) (X(l,M+n),l=1,M)
        enddo
        close(29)
      endif

      if (cube_dens.or.cube_orb.or.cube_elec) then
        call g2g_timer_sum_start('cube gen')
        call cubegen(M15,Xnano)
        call g2g_timer_sum_stop('cube gen')
      endif

c
c-------------------------------------------------
c      endif
      if(DIIS) then
        deallocate (Y,Ytrans,Xtrans,fock,fockm,rho,FP_PFm,
     >  znano,EMAT, bcoef, suma,rho1, scratch, scratch1)
      endif
      deallocate (xnano,rmm5,rmm15)

      if (MEMO) then
        deallocate (kkind,kkinds)
        deallocate(cool,cools)
      endif
      if(allocated(WORK2)) deallocate (WORK2)

c       E=E*627.509391D0

      if(timedep.eq.1) then
        call g2g_timer_sum_start('TD')
        call TD()
        call g2g_timer_sum_stop('TD')
      endif
#ifdef CUBLAS
      call CUBLAS_FREE(devPtrX)
      call CUBLAS_FREE(devPtrY)
      call CUBLAS_SHUTDOWN
#endif
!
!--------------------------------------------------------------------!
      call g2g_timer_stop('SCF')
      call g2g_timer_sum_stop('Finalize SCF')
      call g2g_timer_sum_stop('SCF')
 500  format('SCF TIME ',I6,' sec')
 450  format ('SCF ENERGY = ',F19.12)
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
 900  format(3(F15.9,2x),2x,F15.9)
 777  format(4(F8.4,2x))
 778  format('C',2x,3(F8.4,2x))
 776  format (3(F8.4,2x))
 756  format(2x,I3,2x,f8.4,2x,f8.4,2x,f8.4)
 556  format ('evar',2x,(7(f10.5,2x)))
 345  format(2x,I2,2x,3(f10.6,2x))
 346  format(2x,4(f10.6,2x))
 682  format(2x,f15.10)
  88  format(5(2x,f8.5))
  45  format(E15.6E4)
  91  format(F14.7,4x,F14.7)
c
      !call g2g_timer_sum_stop('SCF');
      return
      end
C  -------------------------
