c SCF subroutine ----------------------------------
c DIRECT VERSION
c calls all integrals generator subroutines : 1 el integrals,
c 2 el integrals, exchange fitting , so it gets S matrix, F matrix
c and P matrix in lower storage mode ( symmetric matrices)
c
c Dario Estrin, 1992
c Modified by
c Nano and Will 2014 (nanolebrero@gmail.com , wagudelos@gmail.com)
c---------------------------------------------------
      subroutine SCFOP(E,dipxyz)
      use garcha_mod
      use mathsubs
      use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC, 
     >     FOCK_ECP_read,FOCK_ECP_write,IzECP
      use faint_cpu77, only: int1, int2, intsol, int3mem, int3lu

      REAL*8:: E2,En,E,Es,Ex,Exc,E1s,Ens
      dimension work(1000)
      real*8, dimension (:,:), ALLOCATABLE ::xnano,znano
      real*8, dimension (:), ALLOCATABLE :: rmm5,rmm15,rmm13,
     >   bcoef_a, bcoef_b, suma, FP_PFa_v, FP_PFb_v
      real*8, dimension (:,:), allocatable :: fock_a,fock_b,fock_am,
     >   fock_bm,rho_a,rho_b,rho1_a,rho1_b,FP_PF_a,FP_PF_b,
     >   FP_PFa_m,FP_PFb_m,EMATa,EMATb,EMAT2_a,EMAT2_b,
     >   Y,Ytrans,Xtrans
c
      integer ndiist
c       dimension d(natom,natom)
      logical  hagodiis,alloqueo, ematalloc
c       REAL*8 , intent(in)  :: qmcoords(3,natom)
c       REAL*8 , intent(in)  :: clcoords(4,nsolin)
      INTEGER :: ErrID,iii,jjj
      LOGICAL :: docholesky
      REAL*8,ALLOCATABLE :: MatrixVec(:),TestMatrix(:)
      INTEGER            :: LWORK2
      REAL*8,ALLOCATABLE :: WORK2(:)
        logical :: just_int3n,ematalloct
      INTEGER :: igpu

      call g2g_timer_start('SCF')
      write(*,*) '======>>>> INGRESO A SCFop <<<<=========='
c------------------------------------------------------------------
c
c Pointers
c
c first P
      MM=M*(M+1)/2 
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M
      M2a=3*M
c
c Number of OM up
      NCOa=NCO
c Number of OM down
      NCOb=NCO+Nunp
c------------------------------------------------
      M1=1
c now F alpha
      M3=M1+MM
c now S, F beta also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues of alpha spin  in open shell case)
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
c Least squares
      M17=M15+MM
c vectors of MO alpha
      M18=M17+MMd
c vectors of MO beta
      M18b=M18+M*NCOa
c weights (in case of using option )
      M19=M18b+M*NCOb
c new Fock matrix alpha

c RAM storage of two-electron integrals (if MEMO=T)
      M20=M19+natom*50*Nang

c new Fock matrix beta
      M21=M20+MM
c eigenvalues (beta spin in open shell case)
      M22=M21+MM
c
      M23 = M22 +  M  

c------------------------------------------------
c Initializations/Defaults
c------------------------------------------------
      allocate (znano(M,M),xnano(M,M))
c      allocate(rmm5(MM),rmm13(m),rmm15(mm))

c      just_int3n = .false.
      alloqueo = .true.
      ematalloc=.false.
      hagodiis=.false.

      E=0.0D0
      E1=0.0D0
      En=0.0D0
      E2=0.0D0
      Es=0.0D0

c      Ndens=1
c      npas=npas+1
c  QUE ES ngeo ?????      
      ngeo=ngeo+1
c Number of electrons
      Nel=2*NCO+Nunp

      good=1.00D0
      niter=0
      D1=1.D0
      D2=1.D0
      DAMP0=GOLD
      DAMP=DAMP0
c
c      Qc=0.0D0
c      do i=1,natom
c        Qc=Qc+Iz(i)
c      enddo
c      Qc=Qc-Nel
c      Qc2=Qc**2
C----------------------------------------
c Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano
c        write(*,*) 'que pasa?'
  
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
 
       call g2g_reload_atom_positions(igrid2)


      call aint_query_gpu_level(igpu)
      if (igpu.gt.1) call aint_new_step()
c
c
c-------------------------------------------------------
c H CORE - 1 electron matrix elements and solvent 
c
      call int1(En)

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

      E1=0.D0
      do k=1,MM
        E1=E1+RMM(k)*RMM(M11+k-1)
      enddo
c
c-------------------------------------------------------

c Diagonalization of S matrix, after this is not needed anymore
c
c      docholesky=.true.

c      IF (docholesky) THEN
c        call g2g_timer_start('cholesky')
c        PRINT*,'DOING CHOLESKY'
c        ALLOCATE(MatrixVec(MM))
c        DO iii=1,MM
c          MatrixVec(iii)=RMM(M5+iii-1)
c        ENDDO

c        CALL DPPTRF('L',M,MatrixVec,ErrID)
c        PRINT*,ErrID
c        ALLOCATE(Y(M,M),Ytrans(M,M))
c        DO iii=1,M;DO jjj=1,M
c          Y(iii,jjj)=0
c          IF (jjj.LE.iii) THEN
c            Y(iii,jjj)=MatrixVec(iii+(2*M-jjj)*(jjj-1)/2)
c          ENDIF
c          Ytrans(jjj,iii)=Y(iii,jjj)
c        ENDDO;ENDDO

c        CALL DTPTRI('L','N',M,MatrixVec,ErrID)
c        PRINT*,ErrID
c        ALLOCATE(Xtrans(M,M))
c        DO iii=1,M;DO jjj=1,M
c          Xtrans(iii,jjj)=0
c          IF (jjj.LE.iii) THEN
c            Xtrans(iii,jjj)=MatrixVec(iii+(2*M-jjj)*(jjj-1)/2)
c          ENDIF
c            X(jjj,iii)=Xtrans(iii,jjj)
c        ENDDO;ENDDO

c        DEALLOCATE(MatrixVec)
c        PRINT*,'CHOLESKY DONE'
c        call g2g_timer_stop('cholesky')
c      ELSE

c ESSL OPTION ------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X,M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
        do ii=1,M; do jj=1,M
          X(ii,jj)=Smat(ii,jj)
        enddo; enddo
        if (allocated(WORK2)) deallocate(WORK2); allocate(WORK2(1))
        call dsyev('V','L',M,X,M,RMM(M13),WORK2,-1,info)
        LWORK2=int(WORK2(1)); deallocate(WORK2); allocate(WORK2(LWORK2))
        call dsyev('V','L',M,X,M,RMM(M13),WORK2,LWORK2,info)
#endif
c-----------------------------------------------------------
c 
c X transformation matrix , canonical orthogonalization
c LINEAR DEPENDENCY ELIMINATION
        allocate (Y(M,M),Ytrans(M,M),Xtrans(M,M))
c
        do j=1,M
          if (RMM(M13+j-1).lt.1.0D-06) then
            write(*,*) 'LINEAR DEPENDENCY DETECTED ACA !!!!'
            do i=1,M
              X(i,j)=0.0D0
              Y(i,j)=0.0D0
            enddo
          else
            do i=1,M
              X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))
              Y(i,j)=X(i,j)*(RMM(M13+j-1))
            enddo
          endif
        enddo
c QUE ES ESTO ????
        do i=1,M
          do j=1,M
            X(i,M2a+j)=X(i,j)
            Ytrans(i,j)=Y(j,i)
            Xtrans(i,j)=X(j,i)
          enddo
        enddo
c      ENDIF  
c
c ======>>>>>> CASE OF NO STARTING GUESS PROVIDED,  <<<<<=========
c   1 E FOCK MATRIX USED
c
      if((.not.VCINP).and.primera) then
        primera=.false.
c QUE HACE ACA ?????
        do i=1,M
          do j=1,M
            X(i,M+j)=0.D0
            do k=1,j
              X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+j+(M2-k)*(k-1)/2-1)
            enddo
c
            do k=j+1,M
              X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+k+(M2-j)*(j-1)/2-1)
            enddo
c
          enddo
        enddo

c Start use of RMM(M5) as Fock matrix (F')

        kk=0
        do j=1,M
          do i=j,M
            kk=kk+1
            RMM(M5+kk-1)=0.D0
            do k=1,M
              RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
            enddo
          enddo
        enddo
        
c
c diagonalization now
c
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M13+i-1)=0.D0
        enddo
c
c ESSL OPTION
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c LAPACK OPTION -----------------------------------------
#ifdef pack
c
c QUE VaRIABleS ENTRAN Y CUALES SALEN ?????
c
        call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
c-----------------------------------------------------------
        do i=1,M
          do j=1,M
            X(i,M2+j)=0.D0
            do k=1,M
c calculando  C=XC'   
              X(i,M2+j)=X(i,M2+j) + X(i,k)*X(k,M+j)
            enddo
          enddo
        enddo
c
c Density Matrix
c alpha and beta coefficients set equal
c

        kk=0
        do k=1,NCOa
          do i=1,M
             kk=kk+1
             MO_coef_at(kk) = X(i,M2+k)
          enddo
        enddo
c
        kk=0
        do k=1,NCOb
          do i=1,M
             kk=kk+1
             MO_coef_at_b(kk) = X(i,M2+k)
          enddo
        enddo
c
        kk=0
        do j=1,M
          do i=j,M
c
            kk=kk+1
            RMM(kk)=0.D0
c
            do k=1,NCOa
              RMM(kk)=RMM(kk)+X(i,M2+k)*X(j,M2+k)
              rhoalpha(kk)=rhoalpha(kk)+X(i,M2+k)*X(j,M2+k)
            enddo
c
            do k=1,NCOb
              RMM(kk)=RMM(kk)+X(i,M2+k)*X(j,M2+k)
              rhobeta(kk)=rhobeta(kk)+X(i,M2+k)*X(j,M2+k)
            enddo
c
            if (i.ne.j) then
              RMM(kk)=2.0D0*RMM(kk)
              rhoalpha(kk)=2.0D0*rhoalpha(kk)
              rhobeta(kk)=2.0D0*rhobeta(kk)
            endif
c
          enddo
        enddo
      endif
c
c
#ifdef PRINT_MATRICES
c------ IMPRIMIENDO DENSIDADES ---------------------------------
      kk=0
      do i=1,M
        do j=i,M
          kk=kk+1 
          write(*,'(I4,X,I4,X,I4,X,F8.5,X,F8.5,X,F8.5)') 
     <          kk,i,j,RMM(kk),rhoalpha(kk),rhobeta(kk)
        enddo
      enddo
#endif
c
c End of Starting guess (No MO , AO known)-------------------------------
!
!
!------------------------------------------------------------------------------!
! Precalculate two-index (density basis) "G" matrix used in density fitting
! here (S_ij in Dunlap, et al JCP 71(8) 1979) into RMM(M7)
! Also, pre-calculate G^-1 if G is not ill-conditioned into RMM(M9)
!
      call int2()
!
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
      if (MEMO) then
         call g2g_timer_start('int3mem')
         call int3mem()   !*1
!         call int3mems()  !*2
         call g2g_timer_stop('int3mem')
      endif
!
! *1) Large elements of t_i put into double-precision cool here
! Size criteria based on size of pre-factor in Gaussian Product Theorem
! (applied to MO basis indices)
!
! *2) Small elements of t_i put into single-precision cools here
!
!------------------------------------------------------------------------------!
c
c-------------------------------------------------------------------      
      if (DIIS.and.alloqueo) then
        alloqueo=.false.

        allocate(rho1_a(M,M),rho_a(M,M),fock_a(M,M),fock_am(MM,ndiis),
     >           rho1_b(M,M),rho_b(M,M),fock_b(M,M),fock_bm(MM,ndiis),
     >           FP_PF_a(M,M), FP_PFa_v(MM),FP_PFa_m(MM,ndiis),
     >           FP_PF_b(M,M), FP_PFb_v(MM),FP_PFb_m(MM,ndiis),
     >           EMATa(ndiis+1,ndiis+1),EMATb(ndiis+1,ndiis+1),
     >           bcoef_a(ndiis+1),bcoef_b(ndiis+1),
     >           suma(MM))
      endif
c-------------------------------------------------------------------
c
c      write(*,*) 'empiezo el loop=========>>>>>>',NMAX

      do 999 while (good.ge.told.and.niter.le.NMAX)
        call g2g_timer_start('Total iter')
        niter=niter+1
       
        if(niter.le.ndiis) then
          ndiist=niter
        else
          ndiist=ndiis
        endif

c        write(*,*)'niter,ndiist,ndiis',niter,ndiist,ndiis

c-------------------------------------------------------
        call int3lu(E2)
        call g2g_solve_groups(0,Ex,0)
c-------------------------------------------------------
        E1=0.0D0
c
c REACTION FIELD CASE --------------------------------------------
c
        call g2g_timer_start('actualiza rmm')
c----------------------------------------------------------------
c E1 includes the solvent electrostatic matrix elements (sol T)
        do k=1,MM
          E1=E1+RMM(k)*RMM(M11+k-1)
        enddo
c
c      DAMP=gold
        if (niter.eq.1) then
          DAMP=0.0D0
        endif

c        write(*,*)"ENTRANDO A DIIS: ",DIIS 
        if (DIIS) then
c
c----  Pasamos matriz de densidad a base ON, antes copio la matriz densidad en la matriz rho  ------
c aca vamos a tener que dividir por dos los terminos no diagonales
c 
          do j=1,M
            do k=1,j
              if (j.eq.k) then
                rho_a(j,k)=rhoalpha(j+(M2-k)*(k-1)/2)
                rho_b(j,k)=rhobeta(j+(M2-k)*(k-1)/2)
              else
                rho_a(j,k)=(rhoalpha(j+(M2-k)*(k-1)/2))/2
                rho_b(j,k)=(rhobeta(j+(M2-k)*(k-1)/2))/2
              endif
            enddo
c
            do k=j+1,M
              rho_a(j,k)=rhoalpha(k+(M2-j)*(j-1)/2)/2
              rho_b(j,k)=rhobeta(k+(M2-j)*(j-1)/2)/2
            enddo
c
          enddo

          rho_a=basechange(M,Ytrans,rho_a,Y)
          rho_b=basechange(M,Ytrans,rho_b,Y)

c
c------------Ahora tenemos rho transformado en la base ON y en forma cuadrada-----------------------------
c
c-------------------------Escritura de fock cuadrada--------------------------------------
c-----------Parte de abajo a la izquierda(incluyendo terminos diagonales)-----------------
c
          do j=1,M
            do k=1,j
              fock_a(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
              fock_b(j,k)=RMM(M3+j+(M2-k)*(k-1)/2-1)
            enddo
c-----------Parte de arriba a la derecha de la matriz (sin incluir terminos diagonales)---
            do k=j+1,M
              fock_a(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
              fock_b(j,k)=RMM(M3+k+(M2-j)*(j-1)/2-1)
            enddo
          enddo


          fock_a=basechange(M,Xtrans,fock_a,X)
          fock_b=basechange(M,Xtrans,fock_b,X)

c--------------En este punto ya tenemos F transformada en base de ON y en su forma cuadrada-----
c
c--------------Acumulamos las matrices de fock (ver si está bien)-------------------------------
c--------fockm es la matriz que va a acumular en sus col. sucesivas matrices de fockes----------
c
          do j=1,ndiis-1
            do i=1,MM
              fock_am(i,j)=fock_am(i,j+1)
              fock_bm(i,j)=fock_bm(i,j+1)
            enddo
          enddo

          do j=1,M
            do k=1,j
              i=j+(M2-k)*(k-1)/2
              fock_am(i,ndiis)=fock_a(j,k)
              fock_bm(i,ndiis)=fock_b(j,k)
            enddo
          enddo
c---------------------------------------------------------------------------------
c--rho(j,k) y fock(j,k) son las matrices densidad y de fock respect (forma cuadrada)--
c---------Calculo de conmutadores [F,P]-------------------------------------------

          FP_PF_a=commutator(fock_a,rho_a)
          FP_PF_b=commutator(fock_b,rho_b)


c---------Pasar Conmutador a vector (guardamos la media matriz de abajo)------------------------------------------------
c#######OJO, SAQUE EL -1########
          do j=1,M
            do k=1,j
              FP_PFa_v(j+(M2-k)*(k-1)/2)=FP_PF_a(j,k)
              FP_PFb_v(j+(M2-k)*(k-1)/2)=FP_PF_b(j,k)
            enddo
          enddo
c----------Acumulamos en las columnas de FP_PFm los sucesivos conmutadores----------

          do j=1,ndiis-1
            do i=1,MM
              FP_PFa_m(i,j)=FP_PFa_m(i,j+1)
              FP_PFb_m(i,j)=FP_PFb_m(i,j+1)
            enddo
          enddo
c
          do i=1,MM
            FP_PFa_m(i,ndiis)=FP_PFa_v(i)
            FP_PFb_m(i,ndiis)=FP_PFb_v(i)
          enddo
        endif
c
c-------------Decidiendo cual critero de convergencia usar-----------
        if (niter.gt.2.and.(DIIS)) then
          hagodiis=.true.
        endif

c<============================
c<======= Alpha case =======
c<============================

        if(.not.hagodiis) then
          if(niter.ge.2) then

c Damping Alpha Matrix
c
            do k=1,MM
              RMM(M5+k-1)=(RMM(M5+k-1) + DAMP*RMM(M20+k-1))/(1.D0+DAMP)
            enddo
          endif
c
c the newly constructed damped matrix is stored, for next iteration
c in RMM(M20) and RMM(M21) for alpha and beta spin respectively.
c
          do k=1,MM
            RMM(M20+k-1)=RMM(M5+k-1)
          enddo

c calculates F'=Xt F X

          do i=1,M
            do j=1,M
              X(i,M+j)=0.D0
              do k=1,j
                X(i,M+j)=X(i,M+j) + X(k,i)*RMM(M5+j+(M2-k)*(k-1)/2-1)
              enddo
              do k=j+1,M
                X(i,M+j)=X(i,M+j) + X(k,i)*RMM(M5+k+(M2-j)*(j-1)/2-1)
              enddo
            enddo
          enddo
c
          kk=0
          do j=1,M
            do i=j,M
              kk=kk+1
              RMM(M5+kk-1)=0.D0
              do k=1,M
                RMM(M5+kk-1)=RMM(M5+kk-1) + X(i,M+k)*X(k,j)
              enddo
            enddo
          enddo
c now F contains transformed F
        endif
 
        if(DIIS) then
          call g2g_timer_start('diis')
          deallocate(EMATa)
          allocate(EMATa(ndiist+1,ndiist+1))

          if(niter.gt.1.and.niter.le.ndiis) then

            EMATa=0
            do k=1,ndiist-1
              do kk=1,ndiist-1
                EMATa(K,KK)=EMAT2_a(K,KK)
              enddo
            enddo
            deallocate (EMAT2_a)

          elseif(niter.gt.ndiis) then

            do k=1,ndiist-1
              do kk=1,ndiist-1
                EMATa(K,KK)=EMAT2_a(K+1,KK+1)
              enddo
            enddo
            deallocate (EMAT2_a)

          endif

          k=ndiist
          do kk=1,ndiist
            kknueva=kk+(ndiis-ndiist)
c-------Escribimos en xnano y znano dos conmutadores de distintas iteraciones------
            do i=1,M
              do j=1,i
                xnano(i,j)=FP_PFa_m(i+(M2-j)*(j-1)/2,ndiis)
                znano(i,j)=FP_PFa_m(i+(M2-j)*(j-1)/2,kknueva)
              enddo
              do j=i+1,M
                xnano(i,j)=FP_PFa_m(j+(M2-i)*(i-1)/2,ndiis)
                znano(i,j)=FP_PFa_m(j+(M2-i)*(i-1)/2,kknueva)
              enddo
            enddo

            call matmuldiag(xnano,znano,rho1_a,M)
c              xnano=matmul(xnano,znano)

            EMATa(ndiist,kk)=0.
            if(kk.ne.ndiist) EMATa(kk,ndiist)=0.
            do l=1,M
              EMATa(ndiist,kk)=EMATa(ndiist,kk)+rho1_a(l,l)
              if (kk.ne.ndiist) then
                EMATa(kk,ndiist)=EMATa(ndiist,kk)
              endif
            enddo
          enddo

          do i=1,ndiist
            EMATa(i,ndiist+1)= -1.0
            EMATa(ndiist+1,i)= -1.0
          enddo
          EMATa(ndiist+1, ndiist+1)= 0.0

          allocate(EMAT2_a(ndiist+1,ndiist+1))
          EMAT2_a=EMATa
c          ematalloct=.true.
c
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
              bcoef_a(i)=0
            enddo
            bcoef_a(ndiist+1)=-1

c----------Cálculo de parámetro optimo para DGELS-------------------
*
*
            LWORK = -1
c            write(*,*) 'work,lwork,ndiist,antes',work(1),lwork,ndiist 
            CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMATa,
     >                 ndiist+1, bcoef_a, ndiist+1, WORK, LWORK, INFO )
            LWORK = MIN( 1000, INT( WORK( 1 ) ) )
c         write(*,*)'work,lwork,ndiist,EMAT',work(1),lwork,ndiist,EMATa
c-----Resuelve la ecuación A*X = B. (EMAT*ci=bcoef). La solución la escribe en bcoef------
            LWORK = 1000
            CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMATa,
     >                 ndiist+1, bcoef_a, ndiist+1, WORK, LWORK, INFO )

c--------Construccion de la "nueva" matriz de fock como cl de las anteriores--------------
c--------Eventualmente se puede probar con la matriz densidad-----------------------------
            suma=0
c       
            do j=1,ndiist
              jnuevo=j+(ndiis-ndiist)
              do i=1,MM
                suma(i)=suma(i)+bcoef_a(j)*fock_am(i,jnuevo)
              enddo
            enddo

            do i=1,MM
              RMM(M5+i-1)=suma(i)
            enddo

              call g2g_timer_stop('diis')

          endif
        endif
c
c now F contains transformed F
c diagonalization now
c
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M13+i-1)=0.D0
        enddo
c
c         if (SHFT) then
c           if (niter.ge.2) then
c adition of level shifts
c constant to diagonal (virtual) elements
c             do i=NCOa+1,M
c               ii=i+(i-1)*(M2-i)/2
c               RMM(M5+ii-1)=RMM(M5+ii-1)+shi
c             enddo
c           endif
c         endif
                 
c ESSL OPTION ------------------------------------------------
        call g2g_timer_start('dspev')
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
        call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
c  
        call g2g_timer_stop('dspev')
c-----------------------------------------------------------
c
c diagonalization now
c
c new coefficients
c
c        call g2g_timer_start('coeff')
        do i=1,M
          do j=1,M
c
            X(i,M2+j)=0.D0
            do k=1,M
              X(i,M2+j)=X(i,M2+j) + X(i,k)*X(k,M+j)
            enddo
          enddo
        enddo
c
c
        kk=0
        do k=1,NCOa
          do i=1,M
            kk=kk+1
            MO_coef_at(kk) = X(i,M2+k)
          enddo
        enddo
c
c xxxxxxx aca poner que escriba ------------------------------
c        if ((good.le.5.0D0*told.and.nopt.eq.0).or.niter.eq.nmax) then
c          do l=1,M
c            do n=1,M
c              X(indexii(l),M+n)=X(l,M2+n)
c            enddo
c          enddo

c          write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, SPIN ALPHA'
c this option in order to look orbitals
c      do n=1,NCO+5
c this other if one is interested in fargment  populations
c          do n=1,M
c            write(29,850) n,RMM(M13+n-1)
c            write(29,400) (X(l,M+n),l=1,M)
c          enddo
c        endif
c--------------------------------------------------------------
c        if (SHFT) then
c          if (niter.ge.2) then
c Level Shifting
c            do i=1,M
c              do j=1,M
c                X(i,j)=X(i,M2+j)
c              enddo
c            enddo
c         endif
c        endif

c<============================
c<======= Beta case =======
c<============================

        if(.not.hagodiis) then
          if(niter.ge.2) then

c Damping  Beta Matrix
c
            do k=1,MM
              RMM(M3+k-1)=(RMM(M3+k-1) + DAMP*RMM(M21+k-1))/(1.D0+DAMP)
            enddo
          endif
c
c the newly constructed damped matrix is stored, for next iteration
c in RMM(M20) and RMM(M21) for alpha and beta spin respectively.
c
          do k=1,MM
            RMM(M21+k-1)=RMM(M3+k-1)
          enddo

c calculates F'=Xt F X

          do i=1,M
            do j=1,M
              X(i,M+j)=0.D0
              do k=1,j
                X(i,M+j)=X(i,M+j)+X(k,M2a+i)*RMM(M3+j+(M2-k)*(k-1)/2-1)
              enddo
c
              do k=j+1,M
                X(i,M+j)=X(i,M+j)+X(k,M2a+i)*RMM(M3+k+(M2-j)*(j-1)/2-1)
              enddo
            enddo
          enddo
c
          kk=0
          do j=1,M
            do i=j,M
              kk=kk+1
              RMM(M3+kk-1)=0.D0
              do k=1,M
                RMM(M3+kk-1)=RMM(M3+kk-1)+X(i,M+k)*X(k,M2a+j)
              enddo
            enddo
          enddo
c now F contains transformed F
        endif

        if(DIIS) then
          call g2g_timer_start('diis')
          deallocate(EMATb)
          allocate(EMATb(ndiist+1,ndiist+1))

          if(niter.gt.1.and.niter.le.ndiis) then
            EMATb=0
            do k=1,ndiist-1
              do kk=1,ndiist-1
                EMATb(K,KK)=EMAT2_b(K,KK)
              enddo
            enddo
            deallocate (EMAT2_b)

          elseif(niter.gt.ndiis) then

            do k=1,ndiist-1
              do kk=1,ndiist-1
                EMATb(K,KK)=EMAT2_b(K+1,KK+1)
              enddo
            enddo
            deallocate (EMAT2_b)

          endif

          k=ndiist
          do kk=1,ndiist
            kknueva=kk+(ndiis-ndiist)
c-------Escribimos en xnano y znano dos conmutadores de distintas iteraciones------
            do i=1,M
              do j=1,i
                xnano(i,j)=FP_PFb_m(i+(M2-j)*(j-1)/2,ndiis)
                znano(i,j)=FP_PFb_m(i+(M2-j)*(j-1)/2,kknueva)
              enddo
              do j=i+1,M
                xnano(i,j)=FP_PFb_m(j+(M2-i)*(i-1)/2,ndiis)
                znano(i,j)=FP_PFb_m(j+(M2-i)*(i-1)/2,kknueva)
              enddo
            enddo

            call matmuldiag(xnano,znano,rho1_b,M)
c              xnano=matmul(xnano,znano)

            EMATb(ndiist,kk)=0.
            if(kk.ne.ndiist) EMATb(kk,ndiist)=0.
            do l=1,M
              EMATb(ndiist,kk)=EMATb(ndiist,kk)+rho1_b(l,l)
              if (kk.ne.ndiist) then
                EMATb(kk,ndiist)=EMATb(ndiist,kk)
              endif
            enddo
          enddo

          do i=1,ndiist
            EMATb(i,ndiist+1)= -1.0
            EMATb(ndiist+1,i)= -1.0
          enddo
          EMATb(ndiist+1, ndiist+1)= 0.0

          allocate(EMAT2_b(ndiist+1,ndiist+1))
          EMAT2_b=EMATb
c          ematalloct=.true.
c
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
              bcoef_b(i)=0
            enddo
            bcoef_b(ndiist+1)=-1

c----------Cálculo de parámetro optimo para DGELS-------------------
*
*
            LWORK = -1
            CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMATb,
     >                 ndiist+1, bcoef_b, ndiist+1, WORK, LWORK, INFO )
            
            LWORK = MIN( 1000, INT( WORK( 1 ) ) )
c           write(*,*) 'work y lwork',work(1),LWORK
c-----Resuelve la ecuación A*X = B. (EMAT*ci=bcoef). La solución la escribe en bcoef------
            LWORK = 1000
            CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMATb,
     >                 ndiist+1, bcoef_b, ndiist+1, WORK, LWORK, INFO )

c--------Construccion de la "nueva" matriz de fock como cl de las anteriores--------------
c--------Eventualmente se puede probar con la matriz densidad-----------------------------
            suma=0
c       
            do j=1,ndiist
              jnuevo=j+(ndiis-ndiist)
              do i=1,MM
                suma(i)=suma(i)+bcoef_b(j)*fock_bm(i,jnuevo)
              enddo
            enddo

            do i=1,MM
              RMM(M3+i-1)=suma(i)
            enddo

          endif
          call g2g_timer_stop('diis')
        endif

c now F contains transformed F
c
c diagonalization now
c
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M22+i-1)=0.D0
        enddo
c
c        if (SHFT) then
c          if (niter.ge.2) then
c adition of level shifts
c constant to diagonal (virtual) elements
c            do i=NCOb+1,M
c              ii=i+(i-1)*(M2-i)/2
c              RMM(M3+ii-1)=RMM(M3+ii-1)+shi
c            enddo
c          endif
c        endif
c
c ESSL OPTION------------------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M3),RMM(M22),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
        call dspev('V','L',M,RMM(M3),RMM(M22),X(1,M+1),M,RMM(M15),info)
#endif
c-----------------------------------------------------------
c diagonalization now
c
c new coefficients
c
        do i=1,M
          do j=1,M
            X(i,M2+j)=0.D0
            do k=1,M
              X(i,M2+j)=X(i,M2+j)+X(i,M2a+k)*X(k,M+j)
            enddo
          enddo
        enddo

        kk=0
        do k=1,NCOb
          do i=1,M
            kk=kk+1
            MO_coef_at_b(kk) = X(i,M2+k)
          enddo
        enddo
c
c xxxxxxx  aca poner que escriba -------------------------------
c       if ((good.le.5.0D0*told.and.nopt.eq.0).or.niter.eq.nmax) then
c         do l=1,M
c           do n=1,M
c             X(indexii(l),M+n)=X(l,M2+n)
c           enddo
c         enddo
c
c         write(29,*)
c         write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, SPIN BETA'
c this option is for looking orbitals
c      do n=1,NCO+5
c this other for fragment populations
c         do n=1,M
c           write(29,850) n,RMM(M22+n-1)
c           write(29,400) (X(l,M+n),l=1,M)
c         enddo
c         rewind(29)
c       endif
c------------------------------------------------------------
c
c       if (SHFT) then
c     if (niter.ge.2) then
c Level Shifting
c          do i=1,M
c            do j=1,M
c              X(i,M2a+j)=X(i,M2+j) 
c            enddo
c          enddo
c        endif
c     endif

c      call g2g_timer_stop('coeff') 
c<<<<<========================================
c<<<<<========================================
c
c-----------------------------------------
c Construction of new density matrix and comparison with old one
c
        kk=0
        good=0.0D0
        do j=1,M
          do i=j,M
            kk  = kk + 1
            tmp = RMM(kk)
            RMM(kk)     = 0.D0
            rhoalpha(kk)= 0.D0
            rhobeta(kk) = 0.D0
c
           do k=1,NCOa
             k0 = M*(k-1)
             ki = k0 + i
             kj = k0 + j
             RMM(kk)     = RMM(kk)     +MO_coef_at(ki)*MO_coef_at(kj)
             rhoalpha(kk)= rhoalpha(kk)+MO_coef_at(ki)*MO_coef_at(kj)
           enddo
c
           do k=1,NCOb
             k0 = M*(k-1)
             ki = k0 + i
             kj = k0 + j
             RMM(kk)    = RMM(kk)    +MO_coef_at_b(ki)*MO_coef_at_b(kj)
             rhobeta(kk)= rhobeta(kk)+MO_coef_at_b(ki)*MO_coef_at_b(kj)
           enddo
c
           if (i.ne.j) then
             RMM(kk)      = 2.0D0*RMM(kk)
             rhoalpha(kk) = 2.0D0*rhoalpha(kk)
             rhobeta(kk)  = 2.0D0*rhobeta(kk)
           endif
c
           del=RMM(kk)-tmp
           if (i.ne.j) then
             del=del*sqrt(2.D0)
           endif
           good=good+del**2
         enddo
       enddo
c
       good=sqrt(good)/float(M)
c
c--- Damping factor update - 
       DAMP=DAMP0
c      IDAMP=0
c      if (IDAMP.EQ.1) then
c        DAMP=DAMP0
c        if (abs(D1).lt.1.D-5) then
c          factor=dmax1(0.90D0,abs(D1/D2))
c          factor=dmin1(factor,1.1D0)
c          DAMP=DAMP0*factor
c        endif
c
c        E=E1+E2+En
c        E=E+Es
c
c        D2=D1
c        D1=(E-E0)
c
c        E0=E
c        DAMP0=DAMP
c      endif

       E=E1+E2+En
       E=E+Es

       write(*,300) niter,DAMP,E
       write(*,"(A,X,F23.16,X,F23.16,X,F23.16,X,F23.16,X,F23.16)")
     >         'En,E1,E2,Ex,Etotal',En,E1,E2,Ex,E1+E2+En+Ex

       call g2g_timer_stop('Total iter')
 999  continue
c 995   continue
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      if (niter.ge.NMAX) then
        write(6,*) 'NO CONVERGENCE AT ',NMAX,' ITERATIONS'
        noconverge=noconverge + 1
        converge=0
c       goto 995
      else
        write(6,*) 'CONVERGED AT',niter,'ITERATIONS'
        noconverge = 0
        converge=converge+1
      endif

      if(noconverge.gt.4) then
        write(6,*)  'stop fon not convergion 4 times'
        stop
      endif


      Es=Es+E1s+Ens
      if (MOD(npas,energy_freq).eq.0) then
      if (GRAD) then

#ifdef G2G
       call g2g_new_grid(igrid)
       call g2g_solve_groups(1, Exc, 0)
#endif


!------------------------------------------------------------------------------!
! Calculation of final energy 

       Es=Ens   ! NucleusQM-CHarges MM

       call int1(En) ! One electron Kinetic (with aint >3) or Kinetic + Nuc-elec (aint >=3)

       if(nsol.gt.0.and.igpu.ge.1) then ! Computing the E1-fock without the MM atoms
          call aint_qmmm_init(0,r,pc)
          call aint_qmmm_fock(E1s,Ens)
       endif
       E1s=0.D0
       do k=1,MM
          E1s=E1s+RMM(k)*RMM(M11+k-1)  ! E1s (here) is the 1e-energy without the MM contribution 
       enddo

       Es=Es+E1-E1s ! Es is the QM/MM energy computated as total 1e - E1s + QMnuc-MMcharge
 


       E=E1+E2+En+Ens+Exc

        write(*,*) "Etot = ", E
        write(*,*) "E1 = ", E1
        write(*,*) "E2 = ", E2
        write(*,*) "Exc = ", Exc


        if (npas.eq.1) npasw = 0
        if (npas.gt.npasw) then
          if (ecpmode) then
               Eecp=0.d0
               do k=1,MM
                 Eecp=Eecp+RMM(k)*(VAAA(k)+VAAB(k)+VBAC(k))
               enddo
          Es=Es-Eecp  ! 
          end if
           call WriteEnergies(E1,E2,En,Eecp,Exc,Es,ecpmode,E_restrain)
          npasw=npas+10
        endif



       else
          E=E-Ex
       endif
       endif
!------------------------------------------------------------------------------!
! Calculation of energy weighted density matrix
!
       kk=M15-1
       do j=1,M
          do i=j,M
             kk=kk+1
             RMM(kk) = 0.D0
             do k=1,NCOa
                k0 = M*(k-1)
                ki = k0 + i
                kj = k0 + j
                RMM(kk) = RMM(kk) - 
     >                    MO_coef_at(ki)*MO_coef_at(kj)*RMM(M13+k-1)
             end do
c
             do k=1,NCOb
                k0 = M*(k-1)
                ki = k0 + i
                kj = k0 + j
                RMM(kk)=RMM(kk) - 
     >                  MO_coef_at_b(ki)*MO_coef_at_b(kj)*RMM(M22+k-1)
             end do
c
             if (i.ne.j) then
                RMM(kk)=2.0D0*RMM(kk)
             endif
c
          end do
       end do
!------------------------------------------------------------------------------!

!
!
!------------------------------------------------------------------------------!
! TODO - PROPERTIES CALCULATION
! calculates dipole moment
! calculates Mulliken poputations
! UNPAIRED SPIN POPULATION
! ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
!------------------------------------------------------------------------------!
!
      if(DIIS) then
        deallocate (Y,Ytrans,Xtrans,rho1_a,rho_a,fock_a,fock_am,
     >           rho1_b,rho_b,fock_b,fock_bm,
     >           FP_PF_a, FP_PFa_v,FP_PFa_m,
     >           FP_PF_b, FP_PFb_v,FP_PFb_m,
     >           EMATa,EMATb,EMAT2_a,EMAT2_b,bcoef_a,bcoef_b,suma)
      endif

      if (MEMO) then
         deallocate(kkind,kkinds)
         deallocate(cool,cools)
      endif

 500  format('SCF TIME ',I6,' sec')
 450  format ('SCF ENERGY = ',F14.7)
 400  format(4(E14.7E2,2x))
 600  format('  ENERGY CONTRIBUTIONS IN A.U.')
 610  format(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
 615  format(2x,'SOLVENT')
 620  format(F14.7,4x,F14.7,4x,F14.7)
 625  format(F14.7)
 760  format(I3,9x,I3,6x,F10.4)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
 300  format(I5,E14.6,2x,F14.7)
 850  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
 900  format(3(F10.4,2x),2x,F10.4)
c---- DEBUGGINGS
c      write(*,*) 'Exc, integrated and calculated',Exc,Ex
c      write(*,*) 'Coulomb energy',E2-Ex
c
      call g2g_timer_stop('SCF')
       return
       end
C  -------------------------                                            
