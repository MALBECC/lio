c TD subroutine ----------------------------------
c DIRECT VERSION
c 
c Dario Estrin, 1992
c Nano, Dario, Uriel, Damian 2012
c La informacion de entrada que recibe esta subrutina es una matriz dencidad inicial obtenida del estado fundamental del sistema junto con sus coordenadas (que tienen que
c ser las mismas que con las que se obtuvo la dencidad de entrada). En el input se especifica tambien el numero de pasos totales de propagacion temporal (nstep) que 
c queremos que haga nuestro programa y el tiempo de integracion de cada paso (tdstep).
c La subrutina esta dividida en dos: La primera parte lee la matriz dencidad inicial, calcula las integrales de overlap y de un electron y escribe las matrices para transformar
c los operadores de fock y dencidad de la base ortonormal a la base de orbitales atomicos y viceversa. Ademas, en el caso de usar excitaciones de Yabana-Bersrch tambien
c se perturba la matriz dencidad inicial en forma de delta de dirac.
c En la segunda parte del programa se hace la propagacion temporal. El esquema es el siguiente:
c
c
c                 rhold(t=i-1) --------------------- 
c                                                  |
c                                                  |
c                rho(t=i)--(int3lu)-->Fock(i)---(Propagador)--->rhonew(t=i+1)
c                     |                                               |
c                     -------------------------------------------------
c
c En cada paso de la propagacion guardamos las componentes del momento dipolar sobre cada uno de los ejes cartecianos.
c-------------------------------------------------------------------------------------------------

       subroutine TD()

c       use latom
       use garcha_mod
      implicit real*8 (a-h,o-z)
c
      integer istep
      complex *8 :: Im,Ix
      real*8 :: t, E2
       real*8, dimension (:,:), ALLOCATABLE ::xnano2,
     > xmm,xtrans,ytrans,Y,fock
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: rho,rhonew,rhold,
     > xnano,rho1
       real*8, dimension (:,:), ALLOCATABLE :: elmu,F1a,F1b
       COMPLEX*8, dimension (:), ALLOCATABLE :: kick, kickc
       real*8, dimension (:), ALLOCATABLE :: ernuc
       dimension q(natom)

c
c auxiliars
c X scratch space in matrices
c
c test
c
c      dimension d(ntq,ntq)

      call g2g_timer_start('TD')
      call g2g_timer_start('inicio')

      just_int3n = false

c------------------------------------------------------------------
c
c Pointers
c
c chequeo -----------
c
      Ndens=1
c---------------------
      E=0.0D0
      E1=0.0D0
      En=0.0D0
      E2=0.0D0
      idip=1
      ngeo=ngeo+1
      Im=(0,2.00000000000)
      sq2=sqrt(2.D0)
      MM=M*(M+1)/2 
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M

c      allocate(kkind(MM))
      allocate (xnano(M,M),xnano2(M,M),fock(M,M),rhonew(M,M),
     >  rhold(M,M),rho(M,M),xmm(M,M),xtrans(M,M),Y(M,M),ytrans(M,M),
     >  rho1(M,M))
c      allocate(d(natom,natom))
      if(propagator.eq.2) allocate (F1a(M,M),F1b(M,M))
c-----------------------------------------------
c      open(unit=49,file='P',status='old')
c------Leemos la matriz dencidad (P) y la guardamos en el RMM(1-->MM)--- Ojo que guardamos la matriz triangular sin multiplicar por 2 a los terminos no diagonales!!!------
c        do j=1,M
c         do k=1,j
c
c          read(49,*) RMM(j+(M2-k)*(k-1)/2)
c
c          enddo
c         do k=j+1,M
c
c          read(49,*) RMM(k+(M2-j)*(j-1)/2)
c
c         enddo
c        enddo
c----------------Copio la matriz dencidad en la matriz rho------------------------
c        open(unit=33,file='rho')
        do j=1,M
         do k=1,j-1
c
          rho(j,k)=RMM(j+(M2-k)*(k-1)/2)/2
c         write(33,*) rho(j,k)
c
         enddo
          rho(j,j)=RMM(j+(M2-k)*(j-1)/2)
         do k=j+1,M
c
          rho(j,k)=RMM(k+(M2-j)*(j-1)/2)/2
c         write(33,*) rho(j,k)
c
         enddo
        enddo
c---------------------------------------------------------------------------------

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
c
* RAM storage of two-electron integrals (if MEMO=T)
      M20 = M19 + natom*50*Nang   
c
      Nel=2*NCO+Nunp
c Initializations/Defaults
c xmm es la primer matriz de (M,M) en el vector X
c
      write(*,*) ' TD CALCULATION  '
      write(*,*) ' Number of functions in the basis set =', M
c----------------------------------
      good=1.00D0
      niter=0
      D1=1.D0
      D2=1.D0
      DAMP0=GOLD
      DAMP=DAMP0
c-----------------------------------
c
       Qc=0.0D0
      do 99 i=1,natom
       Qc=Qc+Iz(i)
 99   continue
      Qc=Qc-Nel
      Qc2=Qc**2
C----------------------------------------
c Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano
        do i=1,natom
         natomc(i)=0
c        write(*,*) natomc(i)
         do j=1,natom

          d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
c         write(*,*) d(i,j),atmin(i),atmin(j)
          zij=atmin(i)+atmin(j)
          ti=atmin(i)/zij
          tj=atmin(j)/zij
          alf=atmin(i)*tj
          rexp=alf*d(i,j)
c         write(*,*) 'vecinos',rexp,rmax
          if (rexp.lt.rmax) then
           natomc(i)=natomc(i)+1
           jatc(natomc(i),i)=j
          endif 
c        
         enddo
c       write(*,*) i,'natomc=',natomc(i)
        enddo
       do ii=nshell(0),1,-1
        nnps(nuc(ii))=ii
       enddo
       do ii=nshell(0)+nshell(1),nshell(0)+1,-1
c
        nnpp(nuc(ii))=ii
       enddo
c
       do ii=M,nshell(0)+nshell(1)+1,-1
        nnpd(nuc(ii))=ii
       enddo
c -------------------------------------------------------------
c H H core, 1 electron matrix elements
c
      call int1(En)
c
c -- SOLVENT CASE --------------------------------------
c      if (sol) then
      call intsol(E1s,Ens,.true.)
c end prueba------
c      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
c      endif
c test ---------------------------------------------------------
         E1=0.D0
        do 302 k=1,MM
 302     E1=E1+RMM(k)*RMM(M11+k-1)
c
c Diagonalization of S matrix, after this is not needed anymore
c s es M13!
c dspev diagonaliza matrices
c ESSL OPTION ------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X,M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
       call dspev('V','L',M,RMM(M5),RMM(M13),X,M,RMM(M15),info)
#endif
c-----------------------------------------------------------
c Paso 1: obtenemos X e Y para las transformaciones.
c LINEAR DEPENDENCY ELIMINATION. Este primer tramo de codigo tiene que ir dentro del loop.
c S es la matriz de overlap
c s es la matriz diagonal de autovalores de S
c U es la matriz con autovectores de S
c donde calcula X=U s^(-1/2)?
       call g2g_timer_start('inicio1')
        do 10 j=1,M
          if (RMM(M13+j-1).lt.1.0D-06) then
           write(*,*) 'LINEAR DEPENDENCY DETECTED'
           do 11 i=1,M
            X(i,j)=0.0D0
 11         Y(i,j)=0.0D0
           else
            do 12 i=1,M
             Y(i,j)=X(i,j)*sqrt(RMM(M13+j-1))
 12          X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))       
          endif
 10      continue
c-------------------xmm va a ser una matriz cuadrada con todos los elementos de matriz de la primer parte de X-------------------
           do i=1,M
            do j=1,M
             xmm(i,j)=X(i,j)
            enddo
           enddo
       do i=1,M
        do j=1,M
         xtrans(j,i)=X(i,j)
         ytrans(j,i)=Y(i,j)
        enddo
       enddo
c-----------Escribimos las componentes del campo externo-------------------------------------
       write(*,*) 'fx =', fx
       write(*,*) 'fy =', fy
       write(*,*) 'fz =', fz
c-------------------------------------------------------------------------------------------
c-REACTION FIELD CASE: Type=Dirac-delta (Yabana-Bertsch PRB 1996) (OFF)-----(en la matriz dencidad en base de o.a.)--------------------------
c-1.Generamos el vector 'ernuc' cuya i-esima posicion indica el producto escalar de las coordenadas del nucleo sobre el que esta centrada la i-esima funcion con el campo externo aplicado.
c-Variables agregadas: ernuc, kick, kickc, ix.
c       allocate (ernuc(M),kick(M),kickc(M))
c        if (field) then
c         if (exter) then
c          ix=(0,1)
c          do i=1,M
c           ernuc(i)=R(nuc(i),1)*fx
c           ernuc(i)=ernuc(i)+(R(nuc(i),2)*fy)
c           ernuc(i)=ernuc(i)+(R(nuc(i),3)*fz)
c-2. Ahora generamos el vector 'kick' y su conjugado 'kickc' para aplicar el campo en forma de delta de dirac: rho(perturbado)=kick*rho*kickc
c           kick(i)=exp(ix*ernuc(i))
c           kickc(i)=exp((-1)*(ix*ernuc(i)))
c          enddo
c-3. Ahora generamos rho(perturbado).
c         do i=1,M
c          do j=1,M
c           rho(i,j)=kick(i)*rho(i,j)
c           rho(i,j)=rho(i,j)*kickc(j)
c          enddo
c         enddo
c        endif
c       endif
c         deallocate(ernuc,kick,kickc)
c-----------Ahora multiplicamos los terminos no diagonales por 2 para tener en el rmm la matriz de rho que necesitamos-----
        if(1.eq.2) then
        do j=1,M
         do k=j,M
          if(j.eq.k) then
           RMM(k+(M2-j)*(j-1)/2)=REAL(rho(j,k))
          else
           RMM(k+(M2-j)*(j-1)/2)=2*(REAL(rho(j,k)))
          endif
         enddo
        enddo
        endif
c---------------Prueba para ver si el RMM(1--->MM) esta bien escrito-------------------------------------------------------
c        open(unit=10,file='RMM1')
c        i=0
c        do j=1,M
c         do k=1,j
c          write(10,*) RMM(j+(M2-k)*(k-1)/2)
c         enddo
c         do k=j+1,M
c          write(10,*) RMM(k+(M2-j)*(j-1)/2)
c         enddo
c        enddo
c---------------transformacion de rho a la base ortonormal---------------- 
c         rho=matmul(ytrans,rho)
c         rho=matmul(rho,y)
c-----------------------transformacion de rho a la base ortonormal con matmulnano-------------------------------

         call matmulnanoc(rho,Y,rho1,M)
         rho=rho1
c-------------------------------------------------------------------------
c-REACTION FIELD CASE: Type=Dirac-delta (Yabana-Bertsch PRB 1996) (OFF)-----(en la matriz dencidad en base de o.m.)--------------------------
c-1.Generamos el vector 'ernuc' cuya i-esima posicion indica el producto escalar de las coordenadas del nucleo sobre el que esta centrada la i-esima funcion con el campo externo aplicado.
c        allocate (ernuc(M),kick(M),kickc(M))
c-Variables agregadas: ernuc, kick, kickc, ix.
c        if (field) then
c         if (exter) then
c          ix=(0,1)
c         do i=1,M
c          write(69,*) nuc(i)
c         enddo
c         do i=1,M
c          ernuc(i)=R(nuc(i),1)*fx
c          ernuc(i)=ernuc(i)+(R(nuc(i),2)*fy)
c          ernuc(i)=ernuc(i)+(R(nuc(i),3)*fz)
c         write(71,*) ernuc(i)
c-2. Ahora generamos el vector 'kick' y su conjugado 'kickc' para aplicar el campo en forma de delta de dirac: rho(perturbado)=kick*rho*kickc
c          kick(i)=exp(ix*ernuc(i))
c          kickc(i)=exp((-1)*(ix*ernuc(i)))
c         write(72,*) kickc(i)
c         enddo
c-3. Ahora generamos rho(perturbado).
c         do i=1,M
c          do j=1,M
c           write(73,*) rho(i,j)
c           rho(i,j)=kick(i)*rho(i,j)
c           rho(i,j)=rho(i,j)*kickc(j)
c          write(74,*) rho(i,j)
c          enddo
c         enddo
c        endif
c       endif
c       deallocate(ernuc,kick,kickc)
c-----------Ahora escribimos en xnano la rho(perturbada) transformada a la base de o.a.-----
c     rho1=matmul(xmm,rho)
c      rho1=matmul(rho1,xtrans)
c      rho1=REAL(rho1)
c---------Ahora escribimos en xnano la rho(perturbada) transformada a la base de o.a. con matmulnano-----
c         call matmulnano(rho,xtrans,rho1,M)
c          rho1=REAL(rho1)
c--------------------------------------------------------------------------------------------------------
c      do j=1,M
c       do k=j,M
c         if(j.eq.k) then
c          RMM(k+(M2-j)*(j-1)/2)=rho1(j,k)
c         else
c          RMM(k+(M2-j)*(j-1)/2)=2*rho1(j,k)
c         endif
c        enddo
c       enddo
c-----Ahora tenemos rho(perturbada) en la base de o.m.-----------------

c---------Prueba----------------------------------------------------------
c         do i=1,M
c          do j=1,M
c           write(31,*) rho(i,j)
c          enddo
c         enddo
c---------------copiamos rho en rhold--------------------------------------------------------
c         do i=1,M
c          do j=1,M
c           rhold(i,j)=rho(i,j)
c          enddo
c         enddo
c-----------------------------------------prueba---------------------------------------------------
c         do j=1,M
c          do k=1,j
c           write(17,*) rhold(j,k)
c          enddo
c          do k=j+1,M
c           write(17,*) rhold(j,k)
c          enddo
c         enddo
c--------------------------------------------------------------------------------------------------
          call g2g_timer_start('int22')
      call int22()
          call g2g_timer_stop('int22')
c--------------------------------------------------------------------------------------------------
c
**
c      if (MEMO) then

      call g2g_timer_start('int3mmem')
         call int3mem()
         call int3mems()
c         call int3lu(E2)
c          write(*,*) 'EEEE2',E2
c      endif
      call g2g_timer_stop('int3mmem')
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
c CASE  SAVING BASIS FUNCTIONS ON DISK
c ONLY IN LEAST SQUARES
c      if (.not.integ) then
c       if (dens) then
c        write(*,*) 'in write'
c      call nwrite(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,NCO,
c     >           NCO,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)
c
c       else
c       write(*,*) 'in write'
c      call write(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,
c     >           NCO,NCO,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)
c
c       endif
c      endif
c LEVEL SHIFT CASE, contruction of initial vectors ------------------
c
      if (SHFT) then
c
       stop 'Level SHFT no CORRESPONDE'
      endif
      call g2g_timer_stop('inicio')

c#####################ACA EMPIEZA LA PROPAGACION TEMPORAL##############################################
      write(*,*) 'PROPAGATION'
c--------------------------------------------------------------------
      do 999 istep=1, ntdstep
      call g2g_timer_start('iteracion')
c---Escribimos el numero de paso------------------------------------
      write (*,*) 'istep =', istep
c----Escribimos el tiempo total de evolucion-----------------------
      t=(istep-1)*tdstep
      t=t*0.02419
      write(*,*) 'evolution time (fs)  =', t
c----------------------------------------------------------------
c
c prueba----------
c      xi(1)=0.0
c      xi(2)=0.0
c      xi(3)=0.0
c      call efield(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c    >            c,a,RMM,xi,V,Ef)
c      write(*,*) Ef(1),Ef(2),Ef(3)
c end prueba------
c
         if(propagator.eq.1.or.(propagator.eq.2.and.istep.lt.200)) then 
         call int3lu(E2)
         call g2g_solve_groups(0,Ex,0)
         endif
c-----------hasta aca-----------------------------------------------
      E1=0.0D0
c#################Tipos de perturbacion###################################
c REACTION FIELD CASE - Type=dirac-delta (OFF)--------------------------------------------------
c
c        if(istep.eq.4.or.istep.eq.5) then
c         if (field) then
c         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,Nel,
c     >       ux,uy,uz)
c
c          if (exter) then
c           g=1.0D0
c           fac=2.54D0
c          else
c           g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
c           Fx=ux/2.54D0
c           Fy=uy/2.54D0
c           Fz=uz/2.54D0
c           fac=(2.54D0*2.00D0)
c          endif
c
c        call dip2(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,g,Fx,Fy,Fz,
c     >     nopt,OPEN)
c
c        E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/fac -
c     >    0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
c         endif
c        endif
c
c        call g2g_timer_start('actualiza rmm')
c
c REACTION FIELD CASE - Type=gaussian (OFF)--------------------------------------------------
c
c        do i=1,MM
c        write(39,*) RMM(M5+i-1)
c        enddo
        if(istep.lt.100) then
         if (field) then
         call dip(ux,uy,uz)
c
          if (exter) then
           g=1.0D0
           fac=2.54D0
           fxx=fx*exp(-0.2*(real(istep-50))**2)
           fyy=fy*exp(-0.2*(real(istep-50))**2)
           fzz=fz*exp(-0.2*(real(istep-50))**2)
c
           write(*,*) fxx,fyy,fzz
          else
           g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
           Fx=ux/2.54D0
           Fy=uy/2.54D0
           Fz=uz/2.54D0
           fac=(2.54D0*2.00D0)
          endif
c
         call dip2(g,Fxx,Fyy,Fzz)
c
        E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/fac -
     >    0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
         endif
        endif
c        do i=1,MM
c        write(38,*) RMM(M5+i-1)
c        enddo
c        stop
c-----------------------------------------------------------------------
c E1 includes solvent 1 electron contributions
        do 303 k=1,MM
 303     E1=E1+RMM(k)*RMM(M11+k-1)
        write(*,*) '1 electron contribution',E1       
c-----------------------------------------------------------------------
c Paso 2: obtengo F en base canonica.
c where U matrix with eigenvectors of S , and s is vector with
c eigenvalues
c ojo: la matriz de fock se construye a partir de la matriz dencidad escrita en la base de orbitales atomicos
        do j=1,M
         do k=1,j
c
          xnano(k,j)=RMM(M5+j+(M2-k)*(k-1)/2-1)
c
          enddo
         do k=j+1,M
c
          xnano(k,j)=RMM(M5+k+(M2-j)*(j-1)/2-1)
c
         enddo
        enddo
        do i=1,M
         do j=1,M
          X(i,M+j)=0.D0
          xnano2(i,j)=X(j,i)
         enddo
        enddo
c
       call g2g_timer_start('actualiza rmm1')
        do 20 j=1,M
         do 20 k=1,M
          do 22 i=1,M
           X(i,M+j)=X(i,M+j)+X(k,i)*xnano(k,j)
 22       continue
c
  20     continue
c
        call g2g_timer_stop('actualiza rmm1')
c
        kk=0
        do i=1,M
         do k=1,M
          xnano(k,i)=X(i,M+k)
         enddo
        enddo
        do 30 j=1,M
         do 30 i=j,M
          kk=kk+1
          RMM(M5+kk-1)=0.D0
         do 32 k=1,M
 32       RMM(M5+kk-1)=RMM(M5+kk-1)+Xnano(k,i)*X(k,j)
 30      continue

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
c--------------Ahora tenemos F transformada en base de ON y en su forma cuadrada---------------
c
         !guardamos cosas para hacer Magnus!!!

         if(propagator.eq.2) then
         if(istep.eq.50) then
 
         F1a=fock         
         endif

         if(istep.eq.150) then

         F1b=fock         

         endif         
         endif
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         E=E1+E2+En
         if (sol) then
          E=E+Es
         endif

         if(propagator.eq.1.or.(propagator.eq.2.and.istep.lt.200)) then 

         write(*,*) 'Leap frog'
       call g2g_timer_start('propagando')
        
c-------En el primer paso hacemos un Euler para atras para generar rhold---->SOLO EN EL PRIMER PASO!!!-----
c         if(istep.eq.1) then
c          rhold=rho+(tdstep*Im*(matmul(fock,rho)))
c          rhold=rhold-(tdstep*Im*(matmul(rho,fock)))
c         endif
          if(istep.eq.1) then
          call conmutc(fock,rho,rhold,M)
          rhold=rho+tdstep*(Im*rhold)
          endif
c----------------------------------------------------------------------------
c######PROPAGACION DE LA MATRIZ DENCIDAD--->PROPAGADOR "LEAP FROG"#######################################
        Im=(0.0d0,2.0d0)
c
c         rhonew=rhold-(tdstep*Im*(matmul(fock,rho)))
c         rhonew=rhonew+(tdstep*Im*(matmul(rho,fock)))
c-----------------------------------------------------------------------------------------------------------
           
          call g2g_timer_start('conmutc 1')

          call conmutc(fock,rho,rhonew,M)
          rhonew=rhold-tdstep*(Im*rhonew)

          call g2g_timer_stop('conmutc 1')
c------------------Actualizamos las rho (rhold-->rho, rho-->rhonew)------------------------------------------
        do i=1,M
         do j=1,M
          rhold(i,j)=rho(i,j)
          rho(i,j)=rhonew(i,j)
         enddo
       enddo
C########################################FIN DEL PROPAGADOR############################################
c--------------------------Para que escriba la matriz dencidad cada 100 pasos (OFF)-------------------------------
c       open(unit=123,file='RHO')
c       if(mod (istep,100) == 0) then
c        do i=1,M
c         do k=1,M
c          write(123,*) rho(i,k)
c         enddo
c        enddo
c       endif
c----------------------------------------------------------------------------------------------------------
c-----------Escribimos la matriz dencidad completa en el ultimo paso para probar idempotencia y traza-------
        if (istep.eq.ntdstep) then
         open(unit=44,file='rholast')
         do j=1,M
          do k=1,j
           write(44,*) rho(j,k)
          enddo
          do k=j+1,M
           write(44,*) rho(j,k)
          enddo
         enddo
        endif
c------------------------------------------------------------------------------------------------------------
c
c------ Aca en xnano va a guardar la matriz dencidad en la base de orbitales atomicos----------------------
c       rho1=matmul(xmm,rho)
c       rho1=matmul(rho1,xtrans)
c       rho1=REAL(rho1)
ct-------Aca en rho1 va a guardar la matriz dencidad en la base de orbitales atomicos con matmulnano---------
       call g2g_timer_start('matmulnanoc')
         call matmulnanoc(rho,xtrans,rho1,M)
         rho1 = REAL(rho1)
       call g2g_timer_stop('matmulnanoc')
       call g2g_timer_stop('propagando')


c------------------------------------------------------------------------------------------------------------
c aca escribe la matriz densidad en el RMM(M1). Esto da lugar al paso 5 en el que la matriz de fock se calcula a partir de la matriz dencidad obtenida en la propagacion.
c------------------Nos quedamos con la matriz triangular de arriva------------------------------------------
        do j=1,M
         do k=j,M
          if(j.eq.k) then
           RMM(k+(M2-j)*(j-1)/2)=rho1(j,k)
          else
           RMM(k+(M2-j)*(j-1)/2)=rho1(j,k)*2
          endif
         enddo
        enddo
c-------------------------------TEST--------------------------------------------------------------------------------------------
c        do j=1,M
c         do k=1,j
c          RMM(j+(M2-k)*(k-1)/2)=xnano(j,k)
c         enddo
c        do k=j,M
c         if(k.eq.j) then
c          RMM(k+(M2-j)*(j-1)/2)=xnano(j,k)
c         else
c          RMM(k+(M2-j)*(j-1)/2)=xnano(j,k)*2
c         endif
c         enddo
c        enddo
c------------------------------Para que escriba la matriz dencidad cada 100 pasos---------------------------       
c        open(unit=124,file='rho')
c        if(mod (istep,100) == 0) then
c         do i=1,MM
c          write(124,*) RMM(i)
c         enddo
c        endif
c-----------Escribimos la matriz dencidad completa en el ultimo paso para probar idempotencia y traza-------
c        if (istep.eq.nstep) then
c         open(unit=44,file='rholast')
c         do j=1,M
c          do k=1,j
c
c           write(44,*) RMM(j+(M2-k)*(k-1)/2)
c
c          enddo
c          do k=j+1,M
c
c           write(44,*) RMM(k+(M2-j)*(j-1)/2)
c
c          enddo
c         enddo
c        endif
c---------------------------FIN TEST-------------------------------------------------------------------------
         else
         write(*,*) 'Magnus'
         write(71,*) 'separo'
         write(72,*) 'separo'
         write(71,*) F1a 
         write(72,*) F1b 
            call predictor(F1a,F1b,fock,rho,Xtrans,X)

            call magnus(fock,rho,rhonew,M,NBCH)

         F1a=F1b
         F1b=fock
         rho=rhonew

         endif

c --------------------CALCULO DE MOMENTO DIPOLAR(NUEVO)-----------------------------------
c calculates dipole moment
c OJO QUE CON SISTEMAS CARGADOS HAY QUE CAMBIAR COSAS!!!! (nano)
       call g2g_timer_start('calculandomu')
       if(istep.eq.1) then
c       open(unit=133,file='mdip')
        open(unit=134,file='x.dip')
        open(unit=135,file='y.dip')
        open(unit=136,file='z.dip')
c
c       write(133,*) '#DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
        write(134,*) '#DIPOLE MOMENT, X COMPONENT (DEBYES) vs t(fs)'
        write(135,*) '#DIPOLE MOMENT, Y COMPONENT (DEBYES) vs t(fs)'
        write(136,*) '#DIPOLE MOMENT, Z COMPONENT (DEBYES) vs t(fs)'

      allocate(elmu(3,MM))

      call dipmem(ux,uy,uz,elmu)
       endif
c       if (mod(istep,10)==0) then
c        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
c     >       Nel,ux,uy,uz)
c      u=sqrt(ux**2+uy**2+uz**2)


c      t=istep*tdstep
c      t=t*0.024
c
        ux=0
        uy=0
        uz=0

       do jjj=1,MM

       ux=ux - elmu(1,jjj)*RMM(jjj)
       uy=uy - elmu(2,jjj)*RMM(jjj)
       uz=uz - elmu(3,jjj)*RMM(jjj)

       enddo

      do 299 i=1,natom
c       Qc=Qc+Iz(i)
       ux=ux+Iz(i)*r(i,1)
       uy=uy+Iz(i)*r(i,2)
       uz=uz+Iz(i)*r(i,3)
 299   continue

        call dip(ux,uy,uz)

        ux=ux*2.54D0
        uy=uy*2.54D0
        uz=uz*2.54D0

        u=sqrt(ux**2+uy**2+uz**2)
c----Para que escriba los componentes del momento dipolar en archivos distintos-------
       write(134,901) t,ux
       write(135,901) t,uy
       write(136,901) t,uz
       call g2g_timer_stop('calculandomu')
c-------------------------------------------------------------------------------------
c u in Debyes
c------------------------------------------------------------------------------------
c calculates Mulliken poputations
       if (ipop.eq.1) then
        call int1(En)
c
c--------------------------------------------------------------
c###########MULLIKEN POPULATION ANALYSIS#######################
c        if(istep.eq.1) then
c        open(unit=125,file='Mulliken')
c        write(125,*) 'MULLIKEN POPULATION ANALYSIS'
c        endif
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
c         write(125,770)
c
        do n=1,natom
c         write(125,760) n,Iz(n),q(n)
        enddo
c
        write(*,*)
        endif
c###########END OF MULLIKEN ANALYSIS###########################
c----------------------------------------------------------------------------------
       if (nopt.ne.3) then
       write(*,300) niter,DAMP,E
       endif
c
c      write(*,*) 'Coulomb E',E2-Ex,Ex
       if (write1) then
c
        open(unit=3,file='restart')
c outputs final  MO ---------------------
c
       do 320 l=1,M
        do 320 n=1,NCO
 320     X(indexii(l),M+n)=X(l,M2+n)
c
       write(3,*) niter,E
c     
       do 325 l=1,M
 325    write(3,400) (X(l,M+n),n=1,NCO)
c
       close(3)
      endif
c
       call g2g_timer_stop('iteracion')
       write(*,*)
 999   continue
c##########ACA TERMINO LA PROPAGACION#############################
 995   continue
c
c
         if (memo) then
         deallocate (kkind,kkinds)
         deallocate(cool,cools)
         deallocate(elmu)        
c
         endif
       if(propagator.eq.2) deallocate (F1a,F1b)
c
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
			 write(*,*) 'ultimo paso SCF'
c       call g2g_timer_start('exchnum')
c
c--------------------------------------------------------------
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
      if (nopt.eq.0) then
c
c PROPERTIES CALCULATION
c calculates dipole moment
c
c
c calculates Mulliken poputations
       if (ipop.eq.1) then
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
         write(*,*) 'MULLIKEN POPULATION ANALYSIS'
         write(*,770)

        do n=1,natom
         write(*,760) n,Iz(n),q(n)
        enddo
c
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
c--------------------------------------------------------------
c outputs final  MO ---------------------
      do 420 l=1,M
c     do 420 n=1,NCO+3
      do 420 n=1,M
 420   X(indexii(l),M+n)=X(l,M2+n)
c
      do 225 l=1,M
 225   write(2,400) (X(l,M+n),n=1,NCO)
c-------------------------------------------------
c writes down MO coefficients and orbital energies
        write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
       do n=1,NCO
        write(29,850) n,RMM(M13+n-1)
        write(29,400) (X(l,M+n),l=1,M)
       enddo
       do n=NCO+1,M
        write(29,851) n,RMM(M13+n-1)
        write(29,400) (X(l,M+n),l=1,M)
       enddo
       close(29)
c
c-------------------------------------------------
      endif
c
c

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
c---- DEBUGGINGS
c      write(*,*) 'Exc, integrated and calculated',Exc,Ex
c      write(*,*) 'Coulomb energy',E2-Ex
c
      call g2g_timer_stop('TD')

        deallocate(xnano,fock,rho,d)
       return
       end
C  -------------------------                                            
