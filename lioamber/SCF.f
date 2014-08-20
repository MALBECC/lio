c SCF subroutine ----------------------------------
c DIRECT VERSION
c Calls all integrals generator subroutines : 1 el integrals,
c 2 el integrals, exchange fitting , so it gets S matrix, F matrix
c and P matrix in lower storage mode ( symmetric matrices)
c 
c Dario Estrin, 1992
c---------------------------------------------------
       subroutine SCF(E,dipxyz)
       use garcha_mod
c      use qmmm_module, only : qmmm_struct, qmmm_nml
c
      implicit real*8 (a-h,o-z)
      integer:: l
       dimension q(natom),work(1000)
       REAL*8 , intent(inout)  :: dipxyz(3) 
       real*8, dimension (:,:), ALLOCATABLE ::xnano,znano
       real*8, dimension (:), ALLOCATABLE :: rmm5,rmm15,rmm13,
     >   bcoef, suma, FP_PFv
       real*8, dimension (:,:), allocatable :: fock,fockm,rho,FP_PF,
     >   FP_PFm,EMAT,Y,Ytrans,Xtrans,rho1,EMAT2,xmm
c
       integer ndiist
c       dimension d(natom,natom)
       logical  hagodiis,alloqueo, ematalloc
c       REAL*8 , intent(in)  :: qmcoords(3,natom)
c       REAL*8 , intent(in)  :: clcoords(4,nsolin)
        INTEGER :: ErrID,iii,jjj
        LOGICAL :: docholesky
        REAL*8,ALLOCATABLE :: MatrixVec(:),TestMatrix(:)
!! CUBLAS
!#ifdef cublas
        integer sizeof_real
        parameter(sizeof_real=8)
        integer stat
        integer*8 devPtrX, devPtrY
        external CUBLAS_INIT, CUBLAS_SET_MATRIX
        external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
        integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX
!#endif
!!
      call g2g_timer_start('SCF')
      just_int3n = .false.
       alloqueo = .true.
       ematalloc=.false.
       hagodiis=.false.
c      if(verbose)  write(6,*) 'ntatom',ntatom,nsol,natom

c------------------------------------------------------------------
c
c Pointers
c
c chequeo -----------
c
      Ndens=1
c---------------------
c       write(*,*) 'M=',M
       allocate (xnano(M,M))
       
       npas=npas+1       
      E=0.0D0
      E1=0.0D0
      En=0.0D0
      E2=0.0D0
      Es=0.0D0
      
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
c
      Nel=2*NCO+Nunp
c

      allocate(rmm5(MM),rmm13(m),rmm15(mm))
      allocate(fock(M,M),xmm(M,M))
cc
           
      good=1.00D0
      niter=0
      D1=1.D0
      D2=1.D0
      DAMP0=GOLD
      DAMP=DAMP0

c
       Qc=0.0D0
       do 99 i=1,natom
       Qc=Qc+Iz(i)
 99    continue
       Qc=Qc-Nel
      Qc2=Qc**2

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

      if (predcoef.and.npas.gt.3) then
      if (.not.OPEN) then
c         rewind 84
        if(verbose) write(*,*) 'prediciendo densidad'
         do i=1,MM
        RMM(i)=(3*old1(i))-(3*old2(i))+(old3(i))

c        write(84,346) old1(i),old2(i),old3(i),RMM(i)
        enddo
       endif
       endif
        
      call int1(En)
       if(nsol.gt.0) then
        call g2g_timer_start('intsol')
      call intsol(E1s,Ens,.true.)
        call g2g_timer_stop('intsol')
         endif
c
c test ---------------------------------------------------------
         E1=0.D0
        do 302 k=1,MM
 302     E1=E1+RMM(k)*RMM(M11+k-1)


c
c Diagonalization of S matrix, after this is not needed anymore
c
       docholesky=.true.
       call g2g_timer_start('cholesky')
       IF (docholesky) THEN
         PRINT*,'DOING CHOLESKY'
         ALLOCATE(MatrixVec(MM))
         DO iii=1,MM
           MatrixVec(iii)=RMM(M5+iii-1)
         ENDDO

         CALL DPPTRF('L',M,MatrixVec,ErrID)
         PRINT*,ErrID
         ALLOCATE(Y(M,M),Ytrans(M,M))
         DO iii=1,M;DO jjj=1,M
           Y(iii,jjj)=0
           IF (jjj.LE.iii) THEN
             Y(iii,jjj)=MatrixVec(iii+(2*M-jjj)*(jjj-1)/2)
           ENDIF
           Ytrans(jjj,iii)=Y(iii,jjj)
         ENDDO;ENDDO

         CALL DTPTRI('L','N',M,MatrixVec,ErrID)
         PRINT*,ErrID
         ALLOCATE(Xtrans(M,M))
         DO iii=1,M;DO jjj=1,M
           Xtrans(iii,jjj)=0
           IF (jjj.LE.iii) THEN
             Xtrans(iii,jjj)=MatrixVec(iii+(2*M-jjj)*(jjj-1)/2)
           ENDIF
           X(jjj,iii)=Xtrans(iii,jjj)
         ENDDO;ENDDO

         DEALLOCATE(MatrixVec)
         PRINT*,'CHOLESKY DONE'
       ELSE


c ESSL OPTION ------------------------------------------
        do i=1,MM
         rmm5(i)=RMM(M5+i-1)
c        write(56,*) RMM(M15+1)
        enddo
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X,M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
       call dspev('V','L',M,RMM5,RMM13,X,M,RMM15,info)
#endif
c-----------------------------------------------------------
c 
c LINEAR DEPENDENCY ELIMINATION
         allocate (Y(M,M),Ytrans(M,M),Xtrans(M,M))
c
         do i=1,MM
         RMM(M5+i-1)=rmm5(i)
c        write(56,*) RMM(M15+1)
         enddo

        do 10 j=1,M
         RMM(M13+j-1)=rmm13(j)
          if (RMM(M13+j-1).lt.1.0D-06) then
          write(*,*) 'LINEAR DEPENDENCY DETECTED'
         do 11 i=1,M
         X(i,j)=0.0D0
 11      Y(i,j)=0.0D0
          else
         do 12 i=1,M
         X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))
 12      Y(i,j)=X(i,j)*(RMM(M13+j-1))
         endif
 10      continue
         do i=1,M
            do j=1,M
            Ytrans(i,j)=Y(j,i)
            Xtrans(i,j)=X(j,i)   
            enddo
         enddo

        ENDIF                   

       call g2g_timer_stop('cholesky')
!! CUBLAS ---------------------------------------------------------------------!
!#ifdef cublas
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
!#endif
!------------------------------------------------------------------------------!

c
c CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
c
       if((.not.ATRHO).and.(.not.VCINP).and.primera) then
        primera=.false.
        do 220 i=1,M
        do 220 j=1,M
         X(i,M+j)=0.D0
         do 222 k=1,j
         X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+j+(M2-k)*(k-1)/2-1)
  222     continue
c
         do 223 k=j+1,M
  223     X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+k+(M2-j)*(j-1)/2-1)
c
  220     continue
c
         kk=0
        do 230 j=1,M
        do 230 i=j,M
         kk=kk+1
         RMM(M5+kk-1)=0.D0
        do 232 k=1,M
 232      RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
 230     continue
c
c diagonalization now
c
        do 249 i=1,M
         RMM(M15+i-1)=0.D0
 249      RMM(M13+i-1)=0.D0

c
c ESSL OPTION
        do i=1,MM
         rmm5(i)=RMM(M5+i-1)
        enddo
           rmm15=0
           rmm13=0       
           xnano=0
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c LAPACK OPTION -----------------------------------------
#ifdef pack
c
       call dspev('V','L',M,RMM5,RMM13,Xnano,M,RMM15,info)
#endif
       do i =1,M
         do j=1,M
        X(i,M+j)=xnano(i,j)
         enddo
          enddo
c-----------------------------------------------------------
         do i=1,MM
         RMM(M5+i-1)=rmm5(i)
        enddo

       do 250 i=1,M
       do 250 j=1,M
c
        X(i,M2+j)=0.D0
       do 252 k=1,M
  252    X(i,M2+j)=X(i,M2+j)+X(i,k)*X(k,M+j)
  250    continue
c
c Density Matrix
c
       kk=0
c
       do 261 k=1,NCO
       do 261 i=1,M
       kk=kk+1
 261     RMM(M18+kk-1)=X(i,M2+k)
c
      kk=0
      do 330 j=1,M
      do 330 i=j,M
      kk=kk+1
      RMM(kk)=0.D0
c
      if(i.eq.j) then
       ff=2.D0
      else
       ff=4.D0
      endif
c
      do 330 k=1,NCO
       RMM(kk)=RMM(kk)+ff*X(i,M2+k)*X(j,M2+k)
 330  continue
c
c
      endif


c End of Starting guess (No MO , AO known)-------------------------------
c
        if ((timedep.eq.1).and.(tdrestart)) then
        call TD()
        return
        endif


      call int22()
c
**
      if (MEMO) then
         call g2g_timer_start('int3mem')
         call int3mem() 
         call int3mems()
         call g2g_timer_stop('int3mem')
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
       allocate(rho1(M,M),rho(M,M),fockm(MM,ndiis),FP_PF(M,M),
     >  FP_PFv(MM),FP_PFm(MM,ndiis),EMAT(ndiis+1,ndiis+1),bcoef(ndiis+1)
     >  ,suma(MM),znano(M,M))
      endif
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c
c      write(*,*) 'empiezo el loop',NMAX
      do 999 while (good.ge.told.and.niter.le.NMAX)
       call g2g_timer_start('Total iter')
      niter=niter+1
       if(niter.le.ndiis) then
         ndiist=niter
        else
         ndiist=ndiis
       endif
c

c
c      if (MEMO) then
            call int3lu(E2)
            call g2g_solve_groups(0,Ex,0)
c-------------------------------------------------------
      E1=0.0D0
c
c REACTION FIELD CASE --------------------------------------------
c
        call g2g_timer_start('actualiza rmm')
c----------------------------------------------------------------
c E1 includes solvent 1 electron contributions
        do 303 k=1,MM
 303     E1=E1+RMM(k)*RMM(M11+k-1)

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

      if (DIIS) then
c--------------Pasamos matriz de densidad a base ON, antes copio la matriz densidad en la matriz rho------
c aca vamos a tener que dividir por dos los terminos no diagonales
c ACA TAMOS write(*,*) 'TAMOS ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        do j=1,M
           do k=1,j
            if (j.eq.k) then
              rho(j,k)=RMM(j+(M2-k)*(k-1)/2)
            else
              rho(j,k)=(RMM(j+(M2-k)*(k-1)/2))/2
             endif
c
           enddo
        do k=j+1,M
              rho(j,k)=RMM(k+(M2-j)*(j-1)/2)/2
c
             enddo
        enddo
!#ifdef cublas
           call g2g_timer_start('cumatmul')
           call cumxtf(rho,devPtrY,rho,M)
           call cumfx(rho,devPtrY,rho,M)
           call g2g_timer_stop('cumatmul')
!#else
! with matmul:
!       rho=matmul(ytrans,rho)
!       rho=matmul(rho,y)
! with matmulnanoc
!            call matmulnanoc(rho,Y,rho,M)
!            rho=rho1
!--------------------------------------!
!#endif
c 
c------------Ahora tenemos rho transformado en la base ON y en forma cuadrada-----------------------------
c-------------------------Escritura de fock cuadrada--------------------------------------
c-----------Parte de abajo a la izquierda(incluyendo terminos diagonales)-----------------
        do j=1,M
        do k=1,j
c
        fock(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
c
        enddo
c-----------Parte de arriba a la derecha de la matriz (sin incluir terminos diagonales)---
        do k=j+1,M
c
        fock(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
c
        enddo
        enddo


!#ifdef cublas
            call cumxtf(fock,devPtrX,fock,M)
            call cumfx(fock,DevPtrX,fock,M)
!#else
!         call matmulnano(fock,X,rho1,M)
!          fock=rho1                     ! RHO1 lo uso como scratch
!#endif


c--------------En este punto ya tenemos F transformada en base de ON y en su forma cuadrada-----
c--------------Acumulamos las matrices de fock (ver si está bien)-------------------------------
c--------fockm es la matriz que va a acumular en sus col. sucesivas matrices de fockes----------

      do j=1,ndiis-1
       do i=1,MM
        fockm(i,j)=fockm(i,j+1)
       enddo

      enddo

       do j=1,M
        do k=1,j
c
         i=j+(M2-k)*(k-1)/2
       fockm(i,ndiis)=fock(j,k)
c
        enddo
        enddo

c---------------------------------------------------------------------------------

c--rho(j,k) y fock(j,k) son las matrices densidad y de fock respect (forma cuadrada)--

c---------Calculo de conmutadores [F,P]-------------------------------------------
 
!#ifdef cublas
         call cuconmut_r(fock,rho,FP_PF,M)
!#else
!         call conmut(fock,rho,FP_PF,M)
!#endif

c---------Pasar Conmutador a vector (guardamos la media matriz de abajo)------------------------------------------------
c#######OJO, SAQUE EL -1########
         do j=1,M
           do k=1,j
         FP_PFv(j+(M2-k)*(k-1)/2)=FP_PF(j,k)
           enddo
         enddo
c----------Acumulamos en las columnas de FP_PFm los sucesivos conmutadores----------
           do j=1,ndiis-1
              do i=1,MM
c
             FP_PFm(i,j)=FP_PFm(i,j+1)
              enddo
           enddo
         do i=1,MM
c            
            FP_PFm(i,ndiis)=FP_PFv(i)
          enddo
      endif
c
c-------------Decidiendo cual critero de convergencia usar-----------
        if (niter.gt.2.and.(DIIS)) then
          hagodiis=.true.
        endif

        if(.not.hagodiis) then
        if(niter.ge.2) then
       do 135 k=1,MM
       kk=M5+k-1
       kk2=M3+k-1
 135   RMM(kk)=(RMM(kk)+DAMP*RMM(kk2))/(1.D0+DAMP)

       endif
        
c the newly constructed damped matrix is stored, for next iteration
c in RMM(M3)
c
      do 137 k=1,MM
       kk=M5+k-1
       kk2=M3+k-1
 137   RMM(kk2)=RMM(kk)
c
c
!----Fijarse si esto esta bien hecho!!
!        do i=1,M
!        do j=1,M
!        X(i,M+j)=0.D0
!        xnano(i,j)=X(j,i)
!        enddo
!        enddo     
!
!
!        do 20 j=1,M
!           do i=1,M
!             X(i,M+j)=0.D0
!           enddo
!         do 22 k=1,j
!         do 22 i=1,M
!        X(i,M+j)=X(i,M+j)+Xnano(i,k)*RMM(M5+j+(M2-k)*(k-1)/2-1)
!  22     continue
!c
!         do 23 k=j+1,M
!         do 23 i=1,M
!  23     X(i,M+j)=X(i,M+j)+Xnano(i,k)*RMM( M5+k+(M2-j)*(j-1)/2-1)
c
!  20     continue
!c
!         kk=0
!        do i=1,M
!        do k=1,M
!        xnano(k,i)=X(i,M+k)
!        enddo
!        enddo
!        do 30 j=1,M
!        do 30 i=j,M
!         kk=kk+1
!         RMM(M5+kk-1)=0.D0
!        do 32 k=1,M
! 32      RMM(M5+kk-1)=RMM(M5+kk-1)+Xnano(k,i)*X(k,j)
! 30     continue
!------------------------------------------------------------------------------!
! Here we obtain the fock matrix in the molecular orbital (MO) basis.
! where U matrix with eigenvectors of S , and s is vector with
! eigenvalues
            fock=0
            call g2g_timer_start('fock')
            do j=1,M
              do k=1,j
                 fock(k,j)=RMM(M5+j+(M2-k)*(k-1)/2-1)
              enddo
              do k=j+1,M
                 fock(k,j)=RMM(M5+k+(M2-j)*(j-1)/2-1)
              enddo
            enddo
!#ifdef cublas
            call cumxtf(fock,devPtrX,fock,M)
            call cumfx(fock,DevPtrX,fock,M)
!#else
!            fock=matmul(xtrans,fock)
!            fock=matmul(fock,xmm)
!#endif
            call g2g_timer_stop('fock')
!
c Fock triangular matrix contained in RMM(M5,M5+1,M5+2,...,M5+MM) is copied to square matrix fock.
            do j=1,M
               do k=1,j
                  RMM(M5+j+(M2-k)*(k-1)/2-1)=fock(j,k)
               enddo
               do k=j+1,M
                  RMM(M5+k+(M2-j)*(j-1)/2-1)=fock(j,k)
               enddo
            enddo
c Now fock is stored in molecular orbital basis.
!-----------------------------------------------!
       endif
c
c now F contains transformed F
c diagonalization now
c
        do 49 i=1,M
         RMM(M15+i-1)=0.D0
 49      RMM(M13+i-1)=0.D0
c
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
             deallocate(EMAT)
             allocate(EMAT(ndiist+1,ndiist+1))
         if(niter.gt.1.and.niter.le.ndiis) then     

         EMAT=0
         do k=1,ndiist-1
         do kk=1,ndiist-1
          EMAT(K,KK)=EMAT2(K,KK)

         enddo
         enddo
         deallocate (EMAT2)
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
!#ifdef cublas
                   call cumatmul_r(xnano,znano,rho1,M)
!#else
!                   call matmuldiag(xnano,znano,rho1,M)
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


             allocate(EMAT2(ndiist+1,ndiist+1))
              EMAT2=EMAT
            ematalloct=.true.
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
      CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, ndiist+1,
     > bcoef, ndiist+1, WORK, LWORK, INFO )
      LWORK = MIN( 1000, INT( WORK( 1 ) ) )


c-----Resuelve la ecuación A*X = B. (EMAT*ci=bcoef). La solución la escribe en bcoef------

      CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, ndiist+1,
     > bcoef, ndiist+1, WORK, LWORK, INFO )

c--------Construccion de la "nueva" matriz de fock como cl de las anteriores--------------
c--------Eventualmente se puede probar con la matriz densidad-----------------------------
         suma=0
c       
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
      endif

  
       call g2g_timer_start('dspev')
c ESSL OPTION ---------------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
       call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
       call g2g_timer_stop('dspev')
       call g2g_timer_start('coeff')

c-----------------------------------------------------------
c
c diagonalization now
c
c
c new coefficients
       do i=1,M
       do k=1,M
           xnano(i,k)=X(k,i)
       enddo
       enddo
c
       do 50 i=1,M
       do 50 j=1,M
c
        X(i,M2+j)=0.D0
       do 52 k=1,M
  52    X(i,M2+j)=X(i,M2+j)+xnano(k,i)*X(k,M+j)
  50    continue
 
       call g2g_timer_stop('coeff') 
       call g2g_timer_start('otras cosas')
c
c --- For the first iteration, damping on density matrix 
c Important for the case of strating guess of AO
c

c
       kk=0
       do 61 k=1,NCO
       do 61 i=1,M
        kk=kk+1
 61     RMM(M18+kk-1)=X(i,M2+k)
c
c Construction of new density matrix and comparison with old one
      kk=0
      good=0.
c
      do 130 j=1,M
      do 130 i=j,M
      kk=kk+1
      tmp=RMM(kk)
      RMM(kk)=0.
c
      if(i.eq.j) then
       ff=2.D0
      else
       ff=4.D0
      endif
      do 139 k=1,NCO
       RMM(kk)=RMM(kk)+ff*X(i,M2+k)*X(j,M2+k)
 139  continue
      del=RMM(kk)-tmp
      if (i.ne.j) then
       del=del*sq2
      endif
      good=good+del**2
 130  continue
c
      good=sqrt(good)/float(M)
      
      if (SHFT) then
c Level Shifting
      do 141 i=1,M
      do 141 j=1,M
       X(i,j)=X(i,M2+j)
 141  continue
      endif
c
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
         fac=dmax1(0.90D0,abs(D1/D2))
         fac=dmin1(fac,1.1D0)
        DAMP=DAMP0*fac
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
       write(*,*) 'good =', good
       call g2g_timer_stop('otras cosas')

       if(verbose) write(6,*) 'iter',niter,'QM Energy=',E+Ex
c
       call g2g_timer_stop('Total iter')
 999   continue
c 995   continue

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
!#ifdef cublas
         call CUBLAS_FREE ( devPtrX )
         call CUBLAS_FREE ( devPtrY )
!#endif
        old3=old2

        old2=old1
c        write(*,*) 'good final',good

       do i=1,MM

         old1(i)=RMM(i)
       enddo


        if(noconverge.gt.4) then 
        write(6,*)  'stop fon not convergion 4 times'
       stop
       endif


c
c -- SOLVENT CASE --------------------------------------
c      if (sol) then
        call g2g_timer_start('intsol 2')
       if(nsol.gt.0) then
      call intsol(E1s,Ens,.false.)
c        write(*,*) 'cosillas',E1s,Ens
        call g2g_timer_stop('intsol 2')
         endif
c      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
      Es=Es+E1s+Ens
c      endif
c--------------------------------------------------------------
       if (GRAD) then
c         if (sol) then
c         endif
c       call g2g_timer_start('exchnum')
#ifdef G2G
#ifdef ULTIMA_CPU
       call exchnum(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >              M18,NCO,Exc,nopt)
#else
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
c calculation of energy weighted density matrix
c
      kk=0
      do 307 j=1,M
      do 307 i=j,M
      kk=kk+1
      RMM(M15+kk-1)=0.D0
c
      if(i.eq.j) then
       ff=2.D0
      else
       ff=4.D0
      endif
c
      do 309 k=1,NCO
       RMM(M15+kk-1)=RMM(M15+kk-1)-RMM(M13+k-1)*ff*X(i,M2+k)*X(j,M2+k)
 309  continue
c
 307   continue
c
c
c
c      if (nopt.eq.1) then
c
c PROPERTIES CALCULATION
c calculates dipole moment
c
       if (idip.eq.1) then
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
       endif
c
c calculates Mulliken poputations
c       if (ipop.eq.1) then
        call int1(En)
c
c--------------------------------------------------------------
        do n=1,natom
         q(n)=Iz(n)
        enddo
c
        do i=1,M
c
        do j=1,i-1
         kk=i+(M2-j)*(j-1)/2
         t0=RMM(kk)*RMM(M5+kk-1)/2.D0
         q(Nuc(i))=q(Nuc(i))-t0
        enddo
c
         kk=i+(M2-i)*(i-1)/2
         t0=RMM(kk)*RMM(M5+kk-1)
         q(Nuc(i))=q(Nuc(i))-t0
c
         do j=i+1,M
         kk=j+(M2-i)*(i-1)/2
         t0=RMM(kk)*RMM(M5+kk-1)/2.D0
         q(Nuc(i))=q(Nuc(i))-t0
         enddo
        enddo
c
         write(85,*) 'MULLIKEN POPULATION ANALYSIS'
         write(85,770)

        do n=1,natom
         write(85,760) n,Iz(n),q(n)
        enddo
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

      rewind 88
      do 420 l=1,M
      do 420 n=1,M
            
 420   X(indexii(l),M+n)=X(l,M2+n)
c
      do 225 l=1,M
 225   write(88,400) (X(l,M+n),n=1,NCO)
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
c
c-------------------------------------------------
c      endif
       if(DIIS) then
       deallocate (fockm,rho,FP_PF,FP_PFv,FP_PFm,
     > znano,EMAT, bcoef,suma,rho1)
       endif
       deallocate(Y,Ytrans,Xtrans,fock,xmm)
       deallocate (xnano,rmm5,rmm13,rmm15)

         deallocate (kkind,kkinds)
         deallocate(cool,cools)
      
c       E=E*627.509391D0 

       if(timedep.eq.1) then
       call TD()
       endif
c
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
  45  format(E14.6E4)
  91  format(F14.7,4x,F14.7)
c
      call g2g_timer_stop('SCF')
       return
       end
C  -------------------------                                            
