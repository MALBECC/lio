!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% intECP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This routines calculate Fock matrix elements for effective core              !
! potential (ECP).                                                             !
!                                                                              !
! EXTERNAL INPUT: type of procedure.                                           !
!   · tipodecalculo:                                                           !
!        · tipodecalculo=1 allocate common variables and make 1 center         !
!          terms integrals (AAA)                                               !
!        · tipodecalculo=2 make 2 center terms integrals (ABB & BBA)           !
!        · tipodecalculo=3 make 3 center terms integrals (ABC & ABA)           !
!                                                                              !
! INTERNAL INPUT: system, basis set and ECP information.                       !
!   · system information:                                                      !
!        · natom: number of QM atoms.                                          !
!   · basis set information:                                                   !
!        · ncont(i): number of contractions of function i-th.                  !
!        · a(i,ii): basis function exponents.                                  !
!        · c(i,ii): basis function coefficients.                               !
!        · Cnorm(i,ii): basis function coefficients, renormalized.             !
!        · nshell(0:3): number of basis functions per shell (s,p,d).           !
!        · Nuc(i): atomic index corresponding to function i.                   !
!        · Lxyz(i,1:3) exponents of the angular part for i function            !
!             · 1 is (x/r), 2 is (y/r) and 3 is (z/r)                          !
!   · ECP information:                                                         !
!        · ecptypes: number of atoms with an ECP                               !
!        · IzECP(i): original nuclear charge of the i-th atom                  ! 
!        · ZlistECP(k): nuclear charge of k-th atom defined in input with ECP  !
!        · Lmax(Z): max angular moment of ECP for the atom with nuclear        !
!          charge Z                                                            !
!        · expnumbersECP(Z,l): contraciones of ECP for nuclear charge Z,       !
!          angular moment l                                                    !
!        · nECP(Z,l): ECP polinomial part exponent                             !
!        · bECP(Z,l): ECP Gaussian part exponents.                             !
!        · aECP(Z,l): ECP coefficients.                                        !
!                                                                              !
!   · ECP calculation variables:                                               !
!        · cutECP: activates cutoff criteria in ECP integrals                  !
!        · cutecp2: turns on cutoff for 2 center integrals                     !
!        · cut2_0: cutoff value for 2 center integrals                         !
!        · cutecp3: turns on cutoff for 3 center integrals                     !
!        · Fulltimer_ECP: turns on timers                                      !
!        · cut3_0: cutoff value for 3 center integrals                         !
!                                                                              !
!   · debug:                                                                   !
!        · ecp_debug: turnd on debug mode                                      !
!        · local_nonlocal: Select local/semilocal contributions only           !
!             · =1 only local term is computed,=2 only semilocal term is       !
!               computed                                                       !
!                                                                              !
! OUTPUTS:                                                                     !
!   · VAAA: 1 center terms integrals of ECP (<A|A|A>)                          !
!   · VAAB: 2 center terms integrals of ECP (<A|A|B>, <B|A|A>)                 !
!   · VBAC: 3 center terms integrals of ECP (<B|A|C> & <B|A|B>)                !
!                                                                              !
! V 2.00 November 2019 Refactor to new 'Pedronic' Lio structure                !
! V 1.00 February 2016 Final version for ECP energy calculations               !
! V 0.95 december 2015 cutoff for ECP interacions, linear scaling              !  
! V 0.91 october 2015 optimized version for energy calculations                !
! V 0.9 september 2015 1st functional version for energy calculations          !
!                                                                              !
! Nicolas Foglia                                                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! References:                                                                  !
! J. Chem. Phys. 65, 3826 (1976); http://dx.doi.org/10.1063/1.432900           !
! J. Chem. Phys. 111, 8778 (1999); http://dx.doi.org/10.1063/1.480225          !
! Foglia N. O. <<Simulacion computacional de dinamica electronica y            !
! reactividad quımica>>. PhD Tesis. University of Buenos Aires, 2019           ! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


module subm_intECP
contains

SUBROUTINE intECP(tipodecalculo)
! Rutina principal, llama a las otras rutinas
! tipodecalculo=1 alocatea variables y calcula terminos de un centro (AAA)
! tipodecalculo=2 calcula terminos de dos centros (ABB)
! tipodecalculo=3 calcula terminos de tres centros (ABC)

   USE ECP_mod, ONLY : defineparams
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: tipodecalculo

   IF (tipodecalculo .EQ. 1) THEN
      CALL intECPAAA() !calcula terminos de un centro (AAA)
   ELSEIF (tipodecalculo .EQ. 2) THEN
      CALL intECPAAB() !calcula terminos de dos centros (AAB)
   ELSEIF (tipodecalculo .EQ. 3) THEN
      CALL intECPABC() !calcula terminos de tres centros (ABC)
   ELSE
      CALL WRITE_POST(4)   
   ENDIF
END SUBROUTINE intECP


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    1 Center terms <Xa|Va|Xa>    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE intECPAAA() !calcula los terminos de Fock para bases y pseudopotenciales en el mismo atomo
   USE basis_data, ONLY : ncont, nshell, nuc
!a(i,ni) exponente de la funcion de base i, contraccion ni
!c(i,ni) coeficiente de la funcion de base i, contraccion ni
!ncont(i) cantidad de contracciones de la funcion de base i
!nshell(i) cantidad de funciones i=1 s, i=2, p, i=3, d
!nuc(i) atomo al que corresponde la base i
   USE ECP_mod, ONLY : ecptypes, IzECP, Lxyz, VAAA,ecp_debug, local_nonlocal, Cnorm,ZlistECP
!nECP, bECP, aECP valores del pseudo potencial
! V = Σ aECP * r^b * exp(-bECP r^2)
! estan guardados como: xECP(Z,l,i) Z carga del nucleo, l momento angular del ecp, i numero de funcion del ecp con Z,l

! ecptypes cantidad de atomos con ECP
! IzECP(i) carga nuclear sin corregir por el Zcore para el atomo i
! Lmax(Z) L maximo del ECP
! Lxyz(i,j) contiene los exponentes de la parte angular de la funcion de base i

! <r|x> = Σ c x^lx y^ly z^lz *e^-ar^2, j=1 lx, j=2, ly, j=3 lz 

!VAAA contiene los terminos <A|A|A> del pseudo potencial

!ecp_debug activa el modo de debug
!local_nonlocal=1 solo calcula parte local,=2 solo calcula parte semilocal

!Cnorm(i,ni) coeficiente de la funcion de base i, contraccion ni correctamente normalizados
!ZlistECP(k) carga del atomo k con ECP

   IMPLICIT NONE
   INTEGER :: i,j,k,ii,ji,pos
   DOUBLE PRECISION :: AAA, acum !variables auxiliares

   INTEGER :: lxi,lxj,lyi,lyj,lzi,lzj !l?$  pootencia de la parte angular de la base
   INTEGER :: M !M cantidad de funciones de base

   AAA=0.d0
   acum=0.d0
   M=nshell(0)+nshell(1)+nshell(2)
   DO i=1, M !barre funciones de la base 
      DO j=1, i !barre el otro coef de la base j<=i ya que la matriz tiene que ser simetrica
         IF (nuc(i) .EQ. nuc(j)) THEN !solo calcula si los terminos corresponden al mismo atomo
            DO k=1, ecptypes !barre atomos con ecp
               IF (IzECP(nuc(i)) .EQ. ZlistECP(k)) THEN !solo calcula si el atomo tiene ECP
                  DO ii=1, ncont(i)
                     DO ji=1, ncont(j) !ji y ii barren contracciones de las funcion de base
                        lxi=Lxyz(i,1)
                        lxj=Lxyz(j,1)
                        lyi=Lxyz(i,2)
                        lyj=Lxyz(j,2)
                        lzi=Lxyz(i,3)
                        lzj=Lxyz(j,3) !exponentes de la parte angular
                        AAA=AAA_LOCAL(i,j,k,ii,ji,lxi+lxj,lyi+lyj,lzi+lzj) !contribucion local al ECP
                        AAA= AAA + AAA_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj) !contribucion semilocal al ECP

                        IF (ecp_debug .AND. local_nonlocal .EQ. 1) THEN ! solo para debugueo
                           AAA=AAA_LOCAL(i,j,k,ii,ji,lxi+lxj,lyi+lyj,lzi+lzj)
                        ELSE IF (ecp_debug .AND. local_nonlocal .EQ. 2) THEN
                           AAA= AAA_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj)
                        END IF

                        acum=acum+AAA*Cnorm(j,ji) !multiplica por el coeficiente de la base
                     END DO
                     pos=i+(1-j)*(j-2*M)/2 !posicion en el vector que guarda la matriz de FOCK triangular
                     VAAA(pos) = VAAA(pos) + acum*Cnorm(i,ii)
                     acum=0.d0
                  END DO
               END IF
            END DO
         END IF
      END DO
   END DO
END SUBROUTINE intECPAAA


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION AAA_LOCAL(i,j,k,ii,ji,lx,ly,lz)!       ͚
!Calcula el termino local del pseudopotencial  [<xi|V(LM)|xj>] = a(LM) ʃr^n exp(-alpha *r^2) dr *
!                                                                      ̊ 
! * ʃ((x/r)^ni (y/r)^li (z/r)^mi dΩ

!los coef de la base se multiplican en la rutina que llama a esta
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, angularint
!nECP, bECP, aECP valores del pseudo potencial
! V = Σ aECP * r^b * exp(-bECP r^2)
! estan escritos como: ?ECP(Z,l,i) Z carga del nucleo, l momento angular del ecp, i numero de funcion del ecp con Z,l
!ZlistECP(k) carga del atomo k con ECP
!Lmax(Z) maximo momento angular del pseudopotencial para el atomo con carga nuclear Z
!expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y momento angular l 
   IMPLICIT NONE
   INTEGER :: n,Z,l
   INTEGER, INTENT(IN) :: i,j,k,ii,ji,lx,ly,lz
!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
!Z carga del nucleo

   DOUBLE PRECISION :: Ccoef !exponente de la gaussiana para la integral angular
   INTEGER :: w !variable auxiliar

   Z=ZlistECP(k)
   L=Lmax(Z)
   n=lx+ly+lz
   AAA_LOCAL=0.d0
   DO w =1, expnumbersECP(z,l) !barre todos los terminos del Lmaximo
      Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
      AAA_LOCAL=AAA_LOCAL+aECP(z,L,w)*angularint(lx,ly,lz)*Q0(n+nECP(z,l,w),Ccoef)
   END DO
   RETURN
END FUNCTION AAA_LOCAL

DOUBLE PRECISION FUNCTION AAA_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj)
!calcula el termino semilocal del pseudopotencial 
!  l                                    l      ͚
!  Σ[<xi|lm> V(l-LM) <lm|xj>] = C(l-LM) Σ    ʃr^n exp(-alpha *r^2) dr *
! m=-l                                 m=-l  ̊ 

! * ʃ((x/r)^ni (y/r)^li (z/r)^mi Ylm dΩ * ʃ((x/r)^nj (y/r)^lj (z/r)^mj Ylm dΩ

!los coef de la base se multiplican en la rutina que llama a esta

   USE basis_data, ONLY : a
   USE ECP_mod, ONLY :ZlistECP,Lmax,aECP,nECP,bECP, expnumbersECP

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
!lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
   INTEGER :: Z !carga del nucleo
   DOUBLE PRECISION :: A2, Acoef
   INTEGER :: l,m,term,n !variables auxiliares

   AAA_SEMILOCAL=0.d0
   A2=0.d0
   Z=ZlistECP(k)
   n=lxi+lxj+lyi+lyj+lzi+lzj

   DO l = 0 , Lmax(z)-1 !barre todos los l de la parte no local
      DO m=-l,l !barre m
         A2=A2+Aintegral(l,m,lxi,lyi,lzi)*Aintegral(l,m,lxj,lyj,lzj)!A2 contiene la parte angular de la integral
      END DO

      IF ( A2 .NE. 0.d0) THEN !solo calcula cuando la parte angular no es cero
         DO term=1, expnumbersECP(z,l) !barre contracciones del ECP para el atomo con carga z y l del ecp
            Acoef=bECP(z,L,term)+a(i,ii)+a(j,ji) !Acoef es el exponente de la integral radial
            AAA_SEMILOCAL=AAA_SEMILOCAL+A2*aECP(z,L,term)*Q0(n+nECP(z,l,term),Acoef)
!                ͚
! Q0(n,alpha)= ʃ r^n exp(-alpha *r^2)dr
!              ̊ 
!aECP coeficiente del ECP
         END DO
         A2=0.d0
      END IF
   END DO
   RETURN
END FUNCTION AAA_SEMILOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    2 Center terms <Xa|Va|Xb>    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!calcula los terminos de Fock en el caso en que una base este centrada en un atomo con ECP
!<Xa|Va|Xb>,<Xb|Va|Xa> NO incluye el caso <Xb|Va|Xb>


SUBROUTINE intECPAAB()
   USE basis_data, ONLY : nshell,nuc,a,ncont
!nshell(i) cantidad de funciones i=1 s, i=2, p, i=3, d
!nuc(i) atomo al que corresponde la base i
!a(i,ni) exponente de la funcion de base i, contraccion ni
!ncont(i) cantidad de contracciones de la funcion de base i
   USE ECP_mod, ONLY : ecptypes,cutECP,IzECP,cutecp2,distx, disty, distz,Lxyz, VAAB,ecp_debug, local_nonlocal, Cnorm, &
   Fulltimer_ECP,tsemilocal,tlocal,Tiempo,tQ1,cut2_0,ZlistECP,pi
! ecptypes cantidad de atomos con ECP
! cutECP activa cuts en las integrales de ECP
! IzECP(i) carga nuclear sin corregir por el Zcore para el atomo i
! cutecp2 cut limite para evitar 0 * NAN
! dist? distancia (Bohr) en la direccion ? entre los atomos i y j dist?(i,j)=?i-?j
! Lxyz(i,j) exponentes de la parte angular de la funcion de base i
! VAAB terminos <A|A|B> y <B|A|A> del pseudo potencial
! ecp_debug activa el modo de debug
! local_nonlocal=1 solo calcula parte local,=2 solo calcula parte semilocal
! Cnorm coeficientes de la base renormalizados
! Fulltimer_ECP activa los timers para int. radiales
! tsemilocal,tlocal,Tiempo,tQ1 auxiliares para el calculo de tiempos
! cut2_0 valor de corte para las integrales de 2 centros, solo alcula si [(xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2]*a(j,ji) <= cut2_0 
! ZlistECP(k) carga del atomo k con ECP

   IMPLICIT NONE
   INTEGER :: i,j,k,M,ii,ji,lxi,lxj,lyi,lyj,lzi,lzj,pos
!M cantidad total de funciones de base
   DOUBLE PRECISION :: Distcoef, AAB, acum, dx,dy,dz
   IF (Fulltimer_ECP) THEN
      tsemilocal=0.d0
      tlocal=0.d0
      tQ1=0.d0
   END IF

   M=nshell(0)+nshell(1)+nshell(2)
   Distcoef=0.d0
   AAB=0.d0
   acum=0.d0
   DO i=1, M !barre funciones de la base 
      DO j=1, i !barre el otro coef de la base j<=i ya que la matriz tiene que ser simetrica
         IF (nuc(i) .NE. nuc(j)) THEN !solo calcula si las bases estan en atomos distintos
            DO k=1, ecptypes ! barre atomos con ecp
               IF (IzECP(nuc(i)) .EQ. ZlistECP(k) .OR. IzECP(nuc(j)) .EQ. ZlistECP(k)) THEN !solo calcula si el atomo tiene ECP
                  dx=distx(nuc(i),nuc(j))
                  dy=disty(nuc(i),nuc(j))
                  dz=distz(nuc(i),nuc(j))
                  Distcoef=(dx**2.d0 + dy**2.d0 + dz**2.d0)
                  lxi=Lxyz(i,1)
                  lxj=Lxyz(j,1)
                  lyi=Lxyz(i,2)
                  lyj=Lxyz(j,2)
                  lzi=Lxyz(i,3)
                  lzj=Lxyz(j,3) !exponentes de la parte angular

   IF (IzECP(nuc(i)) .EQ. ZlistECP(k)) THEN !calculo para ECP en i 
      DO ji=1, ncont(j) !barre contracciones de la base j
         IF ( .NOT. cutECP .OR. ((Distcoef*a(j,ji) .LT. cutecp2) .AND. Distcoef*a(j,ji) .LE. cut2_0)) THEN ! cut por contraccion de la base
! solo calcula los terminos que luego se multipliquen por un factor que no sea demasiado pequeño ya que si los atomos quedan muy
! separados el termino VAAB se obtiene como el producto de un numero muy grande por uno que es casi 0
            acum=0.d0
            DO ii=1, ncont(i) !ii barre contracciones de las funcion de base i
               AAB=AAB_LOCAL(i,j,k,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
               AAB= AAB + AAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
               IF (ecp_debug .AND. local_nonlocal .EQ. 1) THEN ! solo para debugueo
                  AAB=AAB_LOCAL(i,j,k,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
               ELSE IF (ecp_debug .AND. local_nonlocal .EQ. 2) THEN
                  AAB=AAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
               END IF
               acum=acum+AAB*Cnorm(i,ii) !multiplica por el coeficiente de la base
               AAB=0.d0
            END DO
            pos=i+(1-j)*(j-2*M)/2 !posicion en el vector que guarda la matriz de FOCK triangular
            VAAB(pos) = VAAB(pos) + acum*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji)) !multiplica por el otro coef de la base
         END IF
      END DO
   END IF

   IF (IzECP(nuc(j)) .EQ. ZlistECP(k)) THEN !calculo para ECP en j
      DO ii=1, ncont(i) ! barre contracciones de las funcion de base i
         IF ( .NOT. cutECP .OR. ((Distcoef*a(i,ii) .LT. cutecp2) .AND.Distcoef*a(i,ii) .le. cut2_0)) THEN !cut por contraccion de la base
! solo calcula los terminos que luego se multipliquen por un factor que no sea demasiado pequeño ya que si los atomos quedan muy
! separados el termino VAAB se obtiene como el producto de un numero muy grande por uno que es casi 0
            acum=0.d0
            DO ji=1, ncont(j) !barre contracciones de las funcion de base j
               AAB=AAB_LOCAL(j,i,k,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)
               AAB=AAB + AAB_SEMILOCAL(j,i,ji,ii,k,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz)

               IF (ecp_debug .AND. local_nonlocal .EQ. 1 .AND. ecp_debug) THEN !solo para debugueo
                  AAB=AAB_LOCAL(j,i,k,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)
               ELSE IF (ecp_debug .AND. local_nonlocal .EQ. 2) THEN
                  AAB=AAB_SEMILOCAL(j,i,ji,ii,k,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz)
               END IF
               acum=acum+AAB*Cnorm(j,ji) !multiplica por el coeficiente de la base j
               AAB=0.d0
            END DO

            pos=i+(1-j)*(j-2*M)/2 !posicion en el vector que guarda la matriz de FOCK triangular
            VAAB(pos) = VAAB(pos) + acum*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii)) !multiplica por el coeficiente de la base i
         END IF
      END DO
   END IF


               END IF
            END DO
         END IF
      END DO
   END DO

   IF (Fulltimer_ECP) THEN
      Tiempo = tQ1
      CALL WRITE_POST(11)
      Tiempo = tlocal
      CALL WRITE_POST(9)
      Tiempo = tsemilocal
      CALL WRITE_POST(10)
   END IF
END SUBROUTINE intECPAAB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


DOUBLE PRECISION FUNCTION AAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
! calcula el termino semi-local del pseudopotencial centrado en i

!  l                                   
!  Σ[<xi|lm> Vi(l-LM) <lm|xj>] 
! m=-l                               

! los coef de la base se multiplican en la rutina que llama a esta

   USE basis_data, ONLY : a !a(i,ni) exponente de la funcion de base i, contrccion ni
   USE ECP_mod, ONLY :ZlistECP,Lmax,aECP,nECP,bECP, expnumbersECP,Qnl,Fulltimer_ECP,tsemilocal,tQ1
! ZlistECP(k) carga del atomo k con ECP
! Lmax(Z) maximo momento angular del pseudopotencial para el atomo con carga nuclear Z
! Vl= Σ aECP * r^nECP * exp(-bECP r^2)
! expnumbersECP(Z,l) terminos del ECP para el atomo con carga nuclear Z y momento angular l del ECP
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!            ̊ 
! Fulltimer_ECP activa los timers para int. radiales
! tsemilocal,tQ1 auxiliares para el calculo de tiempos

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj
! i,j funciones de la base
! ii,ji numero de contraccion de la funcion
! k atomo con ecp
! lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz ! dx, dy, dz distancia entre nucleos
   DOUBLE PRECISION, DIMENSION(3) :: Kvector

   INTEGER :: l,m, term, lx,ly,lz, lambda,lmaxbase !auxiliades para ciclos
   INTEGER :: Z,n !Z= carga del nucleo
   DOUBLE PRECISION :: acumang, acumint, AABx, AABy, AABz, Kmod,Ccoef, auxdistx,auxdisty,auxdistz
!auxiliares
   INTEGER :: lambmin !minimo valor de lambda para integral angular no nula

   DOUBLE PRECISION :: t1,t2,t1q, t2q !auxiliares para timers

   IF (Fulltimer_ECP) CALL cpu_time ( t1 )
   AAB_SEMILOCAL=0.d0
   Z=ZlistECP(k)
   n=lxi+lxj+lyi+lyj+lzi+lzj
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod=2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)
   lmaxbase=lxj+lyj+lzj
   AABx=0.d0
   AABy=0.d0
   AABz=0.d0
   acumint=0.d0
   acumang=0.d0

   DO l = 0 , Lmax(z)-1 !barre todos los l de la parte no local
      DO term=1, expnumbersECP(z,l) !barre contracciones del ECP para el atomo con carga z y l del ecp
         Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
         Qnl=0.d0
         IF (Fulltimer_ECP) CALL cpu_time ( t1q )
         CALL Qtype1N(Kmod,Ccoef,lmaxbase+l,necp(Z,l,term)+n,necp(Z,l,term)) !calcula integrales radiales
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!  
         IF (Fulltimer_ECP) THEN
            CALL cpu_time ( t2q )
            tQ1=tQ1 +t2q-t1q
         END IF

         DO lx=0,lxj !barre potencias por expansion del binomio de Newton (x - dx)^lxj
            auxdistx=dx**(lxj-lx)
            IF (auxdistx .NE. 0.d0) THEN !si el factor de distancia es 0 deja de calcular
               DO ly=0,lyj
                  auxdisty=dy**(lyj-ly)
                  IF (auxdisty .NE. 0.d0) THEN 
                     DO lz=0,lzj
                        auxdistz=dz**(lzj-lz)
                        IF (auxdistz .NE. 0.d0) THEN 
                           acumint=0.d0
                           lambmin=0

   IF (l-lxj-lyj-lzj .GT. 0) lambmin=l-lxj-lyj-lzj !minimo valor de lambda para integral angular no nula
   DO lambda=lxj+lyj+lzj+l,lambmin,-1
      acumang=0.d0
      DO m=-l,l
         acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
      END DO

      IF (Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda) .EQ. 0.d0) THEN
         WRITE(*,*) necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda,10,lmaxbase+l
         STOP " q = 0 in aab semiloc"
      END IF
      acumint=acumint+acumang*Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda)*aECP(z,L,term)
   END DO
   AABz=AABz+acumint * auxdistz* comb(lzj,lz)
   acumint=0.d0

                        END IF
                     END DO
                     AABy=AABy+AABz * auxdisty * comb(lyj,ly)
                     AABz=0.d0
                  END IF
               END DO
               AABx=AABx+AABy * auxdistx * comb(lxj,lx)
               AABy=0.d0
            END IF
         END DO
         AAB_SEMILOCAL=AAB_SEMILOCAL+AABx
         AABx=0.d0
      END DO
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tsemilocal=tsemilocal+t2-t1
   END IF
   RETURN
END FUNCTION AAB_SEMILOCAL

DOUBLE PRECISION FUNCTION AAB_LOCAL(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz)
!Calcula el termino local del pseudopotencial centrado en i [<xi|Vi(LM)|xj>]  
!los coef de la base se multiplican en la rutina que llama a esta
!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
USE basis_data, ONLY : a
USE ECP_mod, ONLY :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, Qnl, Fulltimer_ECP,tlocal,tQ1
! Vl= Σ aECP * r^nECP * exp(-bECP r^2)
! ZlistECP(k) carga del atomo k con ECP
! Lmax(Z) maximo momento angular del pseudopotencial para el atomo con carga nuclear Z
! expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y l del ECP
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!            ̊ 

! Fulltimer_ECP activa los timers para int. radiales
! semilocal,tQ1 auxiliares para el calculo de tiempos

! angularint(i,j,k)= ʃ(x/r)^i (y/r)^j (z/r)^k dΩ
! tlocal,tQ1 auxiliares para timers
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz
   INTEGER :: z,l,Lmaxbase
! Z carga nuclear
! l maximo valor del momento angular del pseudopotencial
   INTEGER :: w,lxi,lyi,lzi,lambda !Auxiliares
   DOUBLE PRECISION :: Ccoef, acum,integral, Kmod,distcoefx, distcoefy,distcoefz
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: t1,t2,t1q,t2q !auxiliares para timers
   IF (Fulltimer_ECP) CALL cpu_time ( t1 )

   Z=ZlistECP(k)
   L=Lmax(Z)
   Lmaxbase=lx+ly+lz+kxi+kyi+kzi
   AAB_LOCAL=0.d0
   acum=0.d0
   integral=0.d0
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod= 2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)

   DO w =1, expnumbersECP(z,l) !barre todos los terminos del Lmaximo
      Qnl=0.d0
      Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
      IF (Fulltimer_ECP) CALL cpu_time ( t1q )
      CALL Qtype1N(Kmod,Ccoef,Lmaxbase,necp(Z,l,w)+Lmaxbase,necp(Z,l,w)+kxi+kyi+kzi) !calcula integrales radiales

      IF (Fulltimer_ECP) THEN
         CALL cpu_time ( t2q )
         tQ1=tQ1 +t2q-t1q
      END IF

      DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
         distcoefx=dx**(lx-lxi)
         IF ( distcoefx .NE. 0.d0) THEN !si el factor de distancia es 0 deja de calcular
            DO lyi=0,ly
               distcoefy=dy**(ly-lyi)
               IF ( distcoefy .NE. 0.d0) THEN
                  DO lzi=0,lz
                     distcoefz=dz**(lz-lzi)
                     IF ( distcoefz .NE. 0.d0) THEN

   DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
      integral=integral + OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) * Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
      IF (Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda) .EQ. 0.d0)  STOP " q = 0 in aab loc"
   END DO
   acum= acum + integral*distcoefx * distcoefy * distcoefz *comb(lx,lxi) *comb(ly,lyi) * comb(lz,lzi)
   integral=0.d0

                     END IF
                  END DO
               END IF
            END DO
         END IF
      END DO
      AAB_LOCAL=AAB_LOCAL+aECP(z,L,w)*acum 
      acum=0.d0
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tlocal=tlocal+t2-t1
   END IF

   RETURN
END FUNCTION AAB_LOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    3 Center terms <Xb|Va|Xc>    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!calcula los terminos de Fock en el caso en que ninguna base este centrada en el atomo con ECP <Xb|Va|Xc>, <Xb|Va|Xb>


SUBROUTINE intECPABC()
   USE basis_data, ONLY : nshell,nuc,ncont,a
   USE garcha_mod, ONLY : natom
! nshell(i) cantidad de funciones i=1 s, i=2, p, i=3, d
! nuc(i) atomo al que corresponde la base i
! ncont(i) cantidad de contracciones de la funcion de base i
! natom cantidad de atomos en el sistema
! a(i,ni) exponente de la funcion de base i, contrccion ni
   USE ECP_mod, ONLY : pi,ecptypes,cutECP,cutecp3,Lxyz,IzECP,VBAC,local_nonlocal,ecp_debug,distx,disty,distz,Cnorm, Fulltimer_ECP, & 
        tsemilocal,tlocal,Tiempo,tQ1,tQ2,Taux,cut3_0,ZlistECP
! ecptypes cantidad de atomos con ECP
! cutECP activa cuts en las integrales de ECP
! cutecp3 cut limite para evitar 0 * NAN
! Lxyz(i,j) exponentes de la parte angular de la funcion de base i
! IzECP(i) carga nuclear sin corregir por el Zcore para el atomo i
! VBAC terminos <B|A|C> y <B|A|B> del pseudo potencial
! local_nonlocal=1 solo calcula parte local,=2 solo calcula parte semilocal
! ecp_debug activa el modo de debug
! dist? distancia (Bohr) en la direccion ? entre los atomos i y j dist?(i,j)=?i-?j
! Cnorm coeficientes de la base renormalizados
! Fulltimer_ECP activa los timers para int. radiales
! tsemilocal,tlocal,Tiempo,tQ1,tQ2,Taux  auxiliares para el calculo de tiempos
! cut3_0 valor de corte para las integrales de 3 centros
! ZlistECP(k) carga del atomo k con ECP

   IMPLICIT NONE
   INTEGER :: i,j,ii,ji,M,k,ki,pos
   DOUBLE PRECISION :: Distcoef,dxi,dxj,dyi,dyj,dzi,dzj,ABC,acum
   INTEGER :: lxi,lxj,lyi,lyj,lzi,lzj

   IF (Fulltimer_ECP) THEN
      tsemilocal=0.d0
      tlocal=0.d0
      tQ1=0.d0
      tQ2=0.d0
      Taux=0.d0
   END IF

   M=nshell(0)+nshell(1)+nshell(2)
   acum=0.d0
   ABC=0.d0
   DO i=1,M !barre funciones de la base
      DO j=1,i !barre funciones de la base
         DO k=1, natom !barre todos los nucleoas del sistema
            IF (nuc(i) .NE. k .AND. nuc(j) .NE. k) THEN !solo calcula si las 2 funciones de base NO corresponden al atomo con el ECP
               DO ki=1, ecptypes !barre atomos con ecp
                  IF ( IzECP(k) .EQ. ZlistECP(ki)) THEN !solo calcula si el nucleo tiene ecp
                     dxi=-distx(nuc(i),k)
                     dyi=-disty(nuc(i),k)
                     dzi=-distz(nuc(i),k)
                     dxj=-distx(nuc(j),k)
                     dyj=-disty(nuc(j),k)
                     dzj=-distz(nuc(j),k) ! distancias al nucleo con ecp
                     lxi=Lxyz(i,1)
                     lxj=Lxyz(j,1)
                     lyi=Lxyz(i,2)
                     lyj=Lxyz(j,2)
                     lzi=Lxyz(i,3)
                     lzj=Lxyz(j,3)


   DO ii=1, ncont(i) !barre contracciones de la base i
      IF (.NOT.cutECP .OR. (a(i,ii)*(dxi**2+dyi**2+dzi**2) .LE. cut3_0)) THEN
         DO ji=1, ncont(j) !barre contracciones de la base j
            IF (.NOT.cutECP .OR. (a(j,ji)*(dxj**2+dyj**2+dzj**2) .LE. cut3_0)) THEN
               Distcoef=a(i,ii)*(dxi**2.d0 + dyi**2.d0 + dzi**2.d0) + a(j,ji)*(dxj**2.d0 + dyj**2.d0 + dzj**2.d0)

               IF (.NOT.cutECP .OR. ((Distcoef .LT. cutecp3) )) THEN
!cut por exponente de la base
!solo calcula los terminos que luego se multipliquen por un factor que no sea demasiado pequeño para evitar NAN*0
                  ABC=ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
                  ABC= ABC + 4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)

                  IF (local_nonlocal .EQ. 1 .AND. ecp_debug) THEN !solo para debugueo
                     ABC=ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
                  ELSE IF (local_nonlocal .EQ. 2 .AND. ecp_debug) THEN
                     ABC=4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
                  END IF
                  acum=acum + ABC*Cnorm(j,ji)*exp(-Distcoef)
                  ABC=0.d0

               END IF
            END IF
         END DO
         pos=i+(1-j)*(j-2*M)/2
         VBAC(pos) = VBAC(pos) + acum*Cnorm(i,ii)*4.d0*pi
         acum=0.d0
      END IF
   END DO

                  END IF 
               END DO
            END IF
         END DO
      END DO
   END DO

   IF (Fulltimer_ECP) THEN
      Tiempo=tQ1
      CALL WRITE_POST(11)
      Tiempo=tQ2
      CALL WRITE_POST(12)
      Tiempo = tlocal
      CALL WRITE_POST(9)
      Tiempo = tsemilocal
      CALL WRITE_POST(10)
      Tiempo = Taux
      CALL WRITE_POST(13)
   END IF
END SUBROUTINE intECPABC


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


DOUBLE PRECISION FUNCTION ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY :expnumbersECP, Qnl,bECP,IzECP,angularint,pi,Fulltimer_ECP,tlocal,tQ1,Lmax,necp,aECP
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji !terminos de la base
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dx1,dy1,dz1,dx2,dy2,dz2 !distancias del centro con ecp a cada nucleo
   INTEGER, INTENT(IN) :: k !numero de atomo
   INTEGER :: L, Z !L maximo del ecp, Z carga nuclear sin modificar
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: Kmod,Ccoef, integral,auxcomb,auxdist,acum, auxdista, auxdistb, auxdistc, auxdistd, auxdiste, auxdistf
   INTEGER :: lmaxbase

   INTEGER :: ac,bc,cc,dc,ec,fc,w,lambda ! variables auxiliares
   DOUBLE PRECISION :: t1,t2,t1q,t2q ! auxiliares para timers
   IF (Fulltimer_ECP) CALL cpu_time ( t1 )


   ABC_LOCAL=0.d0
   z=IzECP(k)
   L=Lmax(Z)
   lmaxbase=lxi+lyi+lzi+lxj+lyj+lzj
   Kvector=(/a(i,ii)*dx1+a(j,ji)*dx2,a(i,ii)*dy1+a(j,ji)*dy2,a(i,ii)*dz1+a(j,ji)*dz2/)
   Kvector=-2.d0*Kvector
   Kmod=sqrt(Kvector(1)**2.d0+Kvector(2)**2.d0+Kvector(3)**2.d0)
   integral=0.d0
   acum=0.d0

   DO w =1, expnumbersECP(z,l) !barre terminos del ECP para el atomo con carga nuclear Z y l del ECP
      Qnl=0.d0
      Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
      IF (Fulltimer_ECP) CALL cpu_time ( t1q )
      CALL Qtype1N(Kmod,Ccoef,Lmaxbase,necp(Z,l,w)+Lmaxbase+2,necp(Z,l,w)) !calcula integral radial
      IF (Fulltimer_ECP) THEN
         CALL cpu_time ( t2q )
         tQ1=tQ1 +t2q-t1q
      END IF

      DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
         auxdista=dx1**(lxi-ac)
         IF (auxdista .NE. 0.d0) THEN !si el factor de distancia es 0 deja de calcular
            DO bc=0,lyi
               auxdistb=dy1**(lyi-bc)
               IF (auxdistb .NE. 0.d0) THEN
                  DO cc=0,lzi
                     auxdistc=dz1**(lzi-cc)
                        IF (auxdistc .NE. 0.d0) THEN
                           DO dc=0,lxj
                           auxdistd=dx2**(lxj-dc)
                           IF (auxdistd .NE. 0.d0) THEN
                              DO ec=0,lyj
                                 auxdiste=dy2**(lyj-ec)
                                 IF (auxdiste .NE. 0.d0) THEN
                                    DO fc=0,lzj
                                       auxdistf=dz2**(lzj-fc)
                                       IF (auxdistf .NE. 0.d0) THEN

   auxdist=auxdista*auxdistb*auxdistc*auxdistd*auxdiste*auxdistf
   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         integral=integral + OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc) * Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)
         IF (Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda) .NE. Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)) &
         STOP " qnl1l2 = 0 in ABC_SEMILOCAL" !corta el calculo si integral radial es 0
      ELSE
         integral=integral + angularint(ac+dc,bc+ec,cc+fc) * Q0(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),Ccoef) *0.25d0/pi
!parche para el caso accidental en que K fuera = (0,0,0) por compensacion de a(i,ii)*dx1+a(j,ji)*dx2 (idem y,z)
!0.25d0/pi compensa el factor 4pi por el que se multiplica luego a la suma de las integrales
      END IF
   END DO
   auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
   acum=acum + auxcomb*auxdist*integral

   integral=0.d0
   auxcomb=0.d0
   auxdist=0.d0

                                       END IF
                                    END DO
                                 END IF
                              END DO
                           END IF
                        END DO
                     END IF
                  END DO
               END IF
            END DO
         END IF
      END DO

      ABC_LOCAL=ABC_LOCAL+aECP(z,L,w)*acum
      acum=0.d0
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tlocal=tlocal+t2-t1
   END IF
   RETURN
END FUNCTION ABC_LOCAL

DOUBLE PRECISION FUNCTION ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl1l2,necp,bECP,IzECP,Fulltimer_ECP,tsemilocal,tQ2,Taux,Lmax,expnumbersECP,aECP
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,k
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dxi,dyi,dzi,dxj,dyj,dzj !distancias de las nucleos con bases al nucleo con ecp
   INTEGER :: Z,l1max,l2max !Z= carga del nucleo, y momento angular de la base i y j
   DOUBLE PRECISION,DIMENSION (3) :: Kivector,Kjvector
   DOUBLE PRECISION :: Kimod, Kjmod,Ccoef
   INTEGER :: l,term,ac,bc,cc,dc,ec,fc,lambdai, lambdaj,m,lambimin, lambjmin !auxiliares ciclos
   DOUBLE PRECISION :: acumang,acumang1,acumang2,integral,auxcomb,auxdist,acum,auxdista,auxdistb,auxdistc,auxdistd,auxdiste,auxdistf
!auxiliares
   DOUBLE PRECISION :: t1,t2, t1q,t2q,t1aux,t2aux !auxiliares para timers
   IF (Fulltimer_ECP) CALL cpu_time ( t1 )
   t1aux=0.d0
   t2aux=0.d0

   ABC_SEMILOCAL=0.d0
   integral=0.d0
   acum=0.d0
   Qnl1l2=0.d0
   acumang=0.d0
   acumang1=0.d0
   acumang2=0.d0
   auxcomb=0.d0
   auxdist=0.d0

   Z=IzECP(k)

   l1max=lxi+lyi+lzi
   l2max=lxj+lyj+lzj

   Kivector=-2.d0*a(i,ii)*(/dxi,dyi,dzi/)
   Kjvector=-2.d0*a(j,ji)*(/dxj,dyj,dzj/)

   Kimod=sqrt(Kivector(1)**2.d0+Kivector(2)**2.d0+Kivector(3)**2.d0)
   Kjmod=sqrt(Kjvector(1)**2.d0+Kjvector(2)**2.d0+Kjvector(3)**2.d0)

   DO l = 0 , Lmax(z)-1 !barre todos los l de la parte no local
      DO term=1, expnumbersECP(Z,l) !barre contracciones del ECP para el atomo con carga Z y momento angular l del ecp
         Qnl1l2=0.d0
         Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
         IF (Fulltimer_ECP) CALL cpu_time ( t1q )
         call Qtype2N(Kimod,Kjmod,Ccoef,l1max+l,l2max+l,necp(Z,l,term)+l1max+l2max+2,necp(Z,l,term))
!agrega a la matriz Qnl1l2 los terminos correspondientes a un termino radiales.

         IF (Fulltimer_ECP) THEN
            CALL cpu_time ( t2q )
            tQ2=tQ2 +t2q-t1q
         END IF

         DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
            auxdista=dxi**(lxi-ac)
            IF (auxdista .NE. 0.d0) THEN !si el factor de distancia es 0 deja de calcular
               DO bc=0,lyi
                  auxdistb=dyi**(lyi-bc)
                  IF (auxdistb .NE. 0.d0) THEN
                     DO cc=0,lzi
                        auxdistc=dzi**(lzi-cc)
                        IF (auxdistc .NE. 0.d0) THEN
                           DO dc=0,lxj
                              auxdistd=dxj**(lxj-dc)
                              IF (auxdistd .NE. 0.d0) THEN
                                 DO ec=0,lyj
                                    auxdiste=dyj**(lyj-ec)
                                    IF (auxdiste .NE. 0.d0) THEN
                                       DO fc=0,lzj
                                          auxdistf=dzj**(lzj-fc)
                                          IF (auxdistf .NE. 0.d0) THEN

   IF (Fulltimer_ECP) CALL cpu_time ( t1aux )

   auxdist=auxdista*auxdistb*auxdistc*auxdistd*auxdiste*auxdistf
   auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
   lambimin=0
   lambjmin=0

   IF (l-ac-bc-cc .GT. 0) lambimin=l-ac-bc-cc !lambda minimo que no anula la integral angular
   IF (l-dc-ec-fc .GT. 0) lambjmin=l-dc-ec-fc !lambda minimo que no anula la integral angular

   DO lambdai=ac+bc+cc+l,lambimin,-2
      DO lambdaj=dc+ec+fc+l,lambjmin,-2

         acumang=0.d0
         DO m=-l,l
            acumang=acumang+OMEGA2(Kivector,lambdai,l,m,ac,bc,cc)*OMEGA2(Kjvector,lambdaj,l,m,dc,ec,fc)
         END DO
         integral=integral+acumang*Qnl1l2(ac+bc+cc+dc+ec+fc+necp(Z,l,term),lambdai,lambdaj)
         IF (Qnl1l2(ac+bc+cc+dc+ec+fc+necp(Z,l,term),lambdai,lambdaj) .EQ. 0.d0)  STOP " q = 0 in abc semiloc1"

         acumang=0.d0
      END DO
   END DO
   acum=acum + auxcomb*auxdist*integral
   integral=0.d0
   auxcomb=0.d0
   auxdist=0.d0
   IF (Fulltimer_ECP) CALL cpu_time ( t2aux )
   Taux=Taux+t2aux-t1aux

                                          END IF
                                       END DO
                                    END IF
                                 END DO
                              END IF
                           END DO
                        END IF
                     END DO
                  END IF
               END DO
            END IF
         END DO

         ABC_SEMILOCAL=ABC_SEMILOCAL+aECP(Z,L,term)*acum
         acum=0.d0
      END DO
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tsemilocal=tsemilocal+t2-t1
   END IF

   RETURN
END FUNCTION ABC_SEMILOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Functions for Radial Integrals    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION Q0(n,alpha)
!                ͚
! Q0(n,alpha)= ʃ r^n exp(-alpha *r^2)dr
!              ̊ 
   USE ECP_mod, ONLY : pi12,fac,doublefactorial
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN) :: alpha
   IF ( mod(n,2).EQ.0) THEN
      Q0=0.5d0*pi12/sqrt(alpha) * doublefactorial(n-1)/((2.d0*alpha)**(n/2))
       RETURN
   ELSE
      IF ( (n-1)/2 .lt. 0) STOP "Er factorial de un negativo en Q0"
      Q0=fac((n-1)/2)/(2.d0*alpha**((n+1)/2))
      RETURN
   END IF
END FUNCTION Q0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Functions for Angular Integrals    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION Aintegral(l,m,lx,ly,lz)
! calcula angularint(i,j,k)= ʃ(x/r)^i (y/r)^j (z/r)^k Ylm dΩ
! Ylm es el armonico esferico real ortonormal
   USE ECP_mod, ONLY :angularint
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: l,m,lx,ly,lz     
   INTEGER :: ix,iy,iz
   Aintegral=0.d0
   DO ix=0,l
      DO iy=0,l-ix
         iz=l-ix-iy
         Aintegral=Aintegral+Ucoef(l,m,ix,iy,iz)*angularint(lx+ix,ly+iy,lz+iz)
! angularint(i,j,k)= ʃ(x/r)^i (y/r)^j (z/r)^k dΩ
! Ucoef es el coeficiente de la expansion de Ylm en (x/r)^i (y/r)^j (z/r)^k 
      END DO
   END DO
   RETURN
END FUNCTION Aintegral


DOUBLE PRECISION FUNCTION OMEGA1(K,l,a,b,c)
!Esta funcion devuelve el valor de omega evaluado en el vector K
!                    l
!OMEGA1(K,l,a,b,c) = Σ Ylu(K) * ʃ (x/r)^a * (y/r)^b * (z/r)^c Ylu(Ω) dΩ
!                   u=-l

   USE ECP_mod, ONLY : angularint 
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: K
   INTEGER, INTENT(IN) :: l,a,b,c
   DOUBLE PRECISION, DIMENSION(3) :: Kun
   INTEGER :: r,s,t,u
   DOUBLE PRECISION :: SUM1, SUM2
   SUM1=0.d0
   SUM2=0.d0
   OMEGA1=0.d0
   IF ( all(K .EQ. (/0.d0,0.d0,0.d0/))) RETURN !caso especial para los terminos <A|A|A>
   IF (a.LT.0 .OR. b.LT.0 .OR. c.LT.0) RETURN
   Kun=K/sqrt(K(1)**2.d0 + K(2)**2.d0 + K(3)**2.d0)
   DO u=-l,l
      DO r=0,l
         DO s=0,l-r
            t=l-r-s
            SUM1=SUM1+Ucoef(l,u,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
            SUM2=SUM2+Ucoef(l,u,r,s,t)*angularint(a+r,b+s,c+t)
         END DO
      END DO
      OMEGA1=OMEGA1+SUM1*SUM2
      SUM1=0.d0
      SUM2=0.d0
   END DO
   RETURN
END FUNCTION OMEGA1


DOUBLE PRECISION FUNCTION OMEGA2(K,lambda,l,m,a,b,c)
!                           lambda
!OMEGA2(K,lambda,l,m,a,b,c) = Σ Y(lambda,o)(K) * ʃ(x/r)^a * (y/r)^b * (z/r)^c Y(lambda,o)(Ω) Ylm(Ω) dΩ
!                           o=-lambda

   USE ECP_mod, ONLY : angularint
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN),DIMENSION(3) :: K
   INTEGER, INTENT(IN) :: lambda,l,m,a,b,c
   DOUBLE PRECISION, DIMENSION(3) :: Kun
   INTEGER :: o,r,s,t,u,v,w
   DOUBLE PRECISION :: SUM1, SUM2
   Kun=K/sqrt(K(1)**2.d0 + K(2)**2.d0 + K(3)**2.d0)
   SUM1=0.d0
   SUM2=0.d0
   OMEGA2=0.d0

   DO o=-lambda,lambda
      DO r=0,lambda
         DO s=0,lambda-r
            t=lambda-r-s
            SUM1=SUM1+Ucoef(lambda,o,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
            DO u=0,l
               DO v = 0,l-u
                  w=l-u-v
                  SUM2=SUM2+Ucoef(lambda,o,r,s,t)*Ucoef(l,m,u,v,w)*angularint(a+r+u,b+s+v,c+t+w)
               END DO
            END DO
         END DO
      END DO
      OMEGA2=OMEGA2+SUM1*SUM2
      SUM1=0.d0
      SUM2=0.d0
   END DO
END FUNCTION OMEGA2


DOUBLE PRECISION FUNCTION Ucoef(l,m,lx,ly,lz)
   USE ECP_mod, ONLY :  l0,l1,l2,l3,l4,l5
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: l,m,lx,ly,lz

   IF (lx+ly+lz .ne. l) STOP "problem with Ucoef: lx+ly+lz not equal to l"

   IF (l .EQ. 0) THEN
      Ucoef=l0(1)
   ELSEIF (l .EQ. 1) THEN
      Ucoef=l1(int(2*lx+ly+1),m)
   ELSEIF (l .EQ. 2) THEN
      Ucoef=l2(int(0.5*lx*(7-lx)+ly+1),m)
   ELSEIF (l .EQ. 3) THEN
      Ucoef=l3(int(0.5*lx*(9-lx)+ly+1),m)
   ELSEIF (l .EQ. 4) THEN
      Ucoef=l4(int(0.5*lx*(11-lx)+ly+1),m)
   ELSEIF (l .EQ. 5) THEN
      Ucoef=l5(int(0.5*lx*(13-lx)+ly+1),m)
   ELSE
      STOP "ECP error l is greater than 5"
   END IF
   RETURN
END FUNCTION Ucoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Functions for Auxiliar Calc.    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION comb(a,b) !devuelve el combinatorio a,b; a>=b
   USE ECP_mod, ONLY :fac
   INTEGER, INTENT(IN) :: a,b
   IF (a .lt. 0 .or. b .lt. 0) STOP "comb con numeros negativos"
   IF (a .lt. b) STOP "b mayor que a en comb"
   comb=0.d0
   comb=fac(a)/((fac(b)*fac(a-b)))
   RETURN
END FUNCTION comb


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Subr. for Radial Integrals    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
SUBROUTINE Qtype1N(K,Ccoef,lmax,nmax,nmin) !this routine obtain Q(l,n,k,a) where
!              ͚
!Q(l,n,k,a)= ʃ r^n exp(-ar^2) Ml(kr) dr 
!            ̊ 
!where Ml are the modified spherical Bessel function of the first kind.

!agrega a la matriz Qnl los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl ya que hay q llamar a esta rutina por cada termino del pseudopotencial

   USE ECP_mod, ONLY :  alpha, betha, Bn, Cn, Qnl, ecp_full_range_int
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: K,Ccoef

   INTEGER :: lmax !lmax  = 0 para <s||s>, 1 para <s||p>, 2 para <s||d>, ... , 6 para <d||d> 
   INTEGER :: nmin,nmax !nmin y nmax dan el rango en que tiene q calcular Q

!variables auxiliares
   INTEGER :: n,l,i
   DOUBLE PRECISION :: acoef, gam, acum

   IF (ecp_full_range_int) THEN
      nmin=0
      nmax=10
      lmax=4
   END IF

   acoef=K/(2.d0*Ccoef)
   gam=0.5d0*exp(K**2.d0/(4.d0*Ccoef))
   Bn=0.d0
   Cn=0.d0
   CALL ByC(Acoef,Ccoef,nmin-3,nmax,Bn,Cn)

   DO n=nmin,nmax
      DO l=0,lmax
         acum=0.d0
         DO i=l,1,-2
            acum=acum+alpha(l,i)*Bn(n-i)/k**i
         END DO
         DO i=l+1,1,-2
            acum=acum+betha(l+1,i)*Cn(n-i)/k**i
         END DO

         IF (acum .EQ. 0.d0) THEN
            WRITE(*,*) "error en Qtype 1, integral radial 0. n,l",n,l
            STOP
         END IF
         Qnl(n,l)=Qnl(n,l)+acum*gam
         acum=0.d0
      END DO
   END DO
END SUBROUTINE Qtype1N


SUBROUTINE ByC(Acoef,Ccoef,nmin,nmax,Barray,Carray)
!calcula los coeficientes B y C 
   USE ECP_mod, ONLY : DAW,NEXTCOEF,pi12
   USE ESP_FUNCT, ONLY : HYB_DAW_ERR
   IMPLICIT NONE
!acoef,ccoef son los coeficientes para el calculo de B y C
! A(ccoef,acoef)= int exp(-ccoef*x^2)(x+acoef)^n dx from -acoef to inf
!nmin,nmax delimitan el rango de n
!Bn=An(ccoef,acoef) + An(ccoef,-acoef)
!Cn=An(ccoef,acoef) - An(ccoef,-acoef)
   DOUBLE PRECISION, INTENT(IN) :: Acoef,Ccoef
   INTEGER, INTENT(IN) :: nmin,nmax
   DOUBLE PRECISION, intent(INOUT),DIMENSION (-12:14) :: Barray,Carray
!variables auxiliares
   DOUBLE PRECISION :: C0sq,ncos,ca
   INTEGER :: i

   C0sq=sqrt(Ccoef)
   Barray(0)=pi12/c0sq
   Carray(0)=Barray(0)*erf(Acoef*C0sq)

   IF (nmax>0) THEN
      Barray(1)= exp(-Ccoef*Acoef**2)/Ccoef + Acoef*Carray(0)
      Carray(1)= Barray(0)*Acoef
      DO i=2,nmax
         ncos=(i-1)/(2*Ccoef)
         Barray(i)=ncos*Barray(i-2)+Acoef*Carray(i-1)
         Carray(i)=ncos*Carray(i-2)+Acoef*Barray(i-1)
      END DO
   END IF

   IF (nmin<0) THEN
      Barray(-1)=2*pi12*HYB_DAW_ERR(Acoef*C0sq)
      Carray(-1)=2*pi12*DAW(Acoef*C0sq)
      IF (nmin<-1) THEN
         ca=2*Ccoef*Acoef
         Barray(-2)=ca*Carray(-1)-2*Ccoef*Barray(0)
         Carray(-2)=2*ca*exp(-Ccoef*Acoef**2)+ca*Barray(-1)-2*Ccoef*Carray(0)
         DO i=-3,nmin,-1
            Barray(i)=NEXTCOEF(1,i,ca,exp(-Ccoef*Acoef**2),Ccoef,Carray(i+1), Barray(i+2))
            Carray(i)=NEXTCOEF(-1,i,ca,exp(-Ccoef*Acoef**2),Ccoef,Barray(i+1), Carray(i+2))
         END DO
      END IF
   END IF
END SUBROUTINE ByC



SUBROUTINE Qtype2N(Ka,Kb,Ccoef,l1max,l2max,nmax,nmin)
!this routine obtain Q(l1,l2,n,k1,k2,a) where
!                      ͚
!Q(l1,l2,n,k1,k2,a)= ʃ r^n exp(-ar^2) Ml1(k1*r) Ml2(k2*r) dr 
!                    ̊ 
!where Ml are the modified spherical Bessel function of the first kind.

!agrega a la matriz Qnl1l2 los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl1l2 ya que hay que llamar a esta rutina por cada termino del pseudopotencial
   USE ECP_mod, ONLY :  alpha, betha, rho, tau, sigma, sigmaR, Qnl1l2,ecp_full_range_int
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: Ka,Kb,Ccoef
!l1max y l2max = 0 para s, 1 para p, 2 para d, etc      
!n corresponde al exponente del paseudopotencial en r^n
   INTEGER :: l1max, l2max
   INTEGER :: nmin,nmax
!variables auxiliares
   INTEGER :: i,j,n,l1,l2
   DOUBLE PRECISION :: alfok, betok, acum1
        
   IF (ecp_full_range_int) THEN
      nmin=0
      nmax=10
      l1max=4
      l2max=4
   END IF

   CALL integrals(Ka,Kb,Ccoef,nmin-l1max-l2max-2,nmax+1)

   acum1=0.d0
   DO n=nmin,nmax
      DO l1=0,l1max
         DO l2=0,l2max
            DO i=l1,1,-2
               alfok=alpha(l1,i)/Ka**i
               DO j=l2,1,-2
                  IF (tau(n-i-j) == 0.d0) STOP "Error, no calculo tau"
                  acum1=acum1+alpha(l2,j)*tau(n-i-j)/Kb**j
               END DO
               DO j=l2+1,1,-2
                  acum1=acum1+betha(l2+1,j)*sigmaR(n-i-j)/Kb**j
                  IF (sigmaR(n-i-j) == 0.d0) STOP "Error, no calculo sigmaR"
               END DO
               Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*alfok
               acum1=0.d0
            END DO
            DO i=l1+1,1,-2
               betok=betha(l1+1,i)/Ka**i
               DO j=l2,1,-2
                  acum1=acum1+alpha(l2,j)*sigma(n-i-j)/Kb**j
                  IF (sigma(n-i-j) == 0.d0) STOP "Error, no calculo sigma"
               END DO
               DO j=l2+1,1,-2
                  acum1=acum1+betha(l2+1,j)*rho(n-i-j)/Kb**j
                  IF (rho(n-i-j) == 0.d0) STOP "Error, no calculo rho"
               END DO
               Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*betok
               acum1=0.d0
            END DO
         END DO
      END DO
   END DO
END SUBROUTINE Qtype2N



SUBROUTINE integrals(Ka,Kb,Ccoef,nmin,nmax)
!obtain integrals rho, tau, sigma.AND.sigmaR from n between nmin.AND.nmax
! rho(n) = int exp(-cr^2) * sinh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! sigma(n) = int exp(-cr^2) * sinh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf
! sigmaR(n) = int exp(-cr^2) * cosh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! tau(n) = int exp(-cr^2) * cosh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf

   USE ECP_mod, ONLY : Bn1,Bn2,Cn1,Cn2,rho, tau, sigma, sigmaR,ecp_full_range_int
   IMPLICIT NONE
   INTEGER :: nmin,nmax
   DOUBLE PRECISION, INTENT(IN) :: Ka,Kb,Ccoef
   DOUBLE PRECISION, DIMENSION(2) :: acoef,gammacoef
   acoef(1)=0.5d0*(Ka+Kb)/Ccoef
   acoef(2)=0.5d0*abs(Ka-Kb)/Ccoef
   gammacoef(1)=0.25d0*exp(Ccoef*acoef(1)**2.d0)
   gammacoef(2)=0.25d0*exp(Ccoef*acoef(2)**2.d0)
   Bn1=0.d0
   Bn2=0.d0
   Cn1=0.d0
   Cn2=0.d0
   rho=0.d0
   tau=0.d0
   sigma=0.d0
   sigmaR=0.d0

   IF (ecp_full_range_int) THEN
      nmin=-12
      nmax=14
   END IF

   CALL ByC(acoef(1),Ccoef,nmin,nmax,Bn1,Cn1)
   CALL ByC(acoef(2),Ccoef,nmin,nmax,Bn2,Cn2)
   rho=gammacoef(1)*Bn1-gammacoef(2)*Bn2
   tau=gammacoef(1)*Bn1+gammacoef(2)*Bn2
   sigma=gammacoef(1)*Cn1+sign(1.d0,Ka-Kb)*gammacoef(2)*Cn2
   sigmaR=gammacoef(1)*Cn1-sign(1.d0,Ka-Kb)*gammacoef(2)*Cn2
END SUBROUTINE integrals

end module subm_intECP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
