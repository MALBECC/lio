!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Effective Core Potential Module    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Contains commmon variables and functions for ECP subroutines                                                                     !
!                                                                                                                                  !
! V 1.00 February 2016 Final version for ECP energy calculations                                                                   !
! V 0.95 december 2015, optimized Module, and compatible with gfortran                                                             !
! V 0.9  sept 2015, first functional version for energy calculations                                                               !
!                                                                                                                                  !
! Nicolas Foglia                                                                                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"

MODULE ECP_mod
   IMPLICIT NONE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Namelist Variables    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   LOGICAL :: ecpmode !activa los pseudopotenciales
   INTEGER :: ecptypes !cantidad de typos de atomos con ECP
   INTEGER :: ECPatoms !cantidad de atomos con ECP
   CHARACTER (LEN=30) :: tipeECP !tipo de ECP usado, tiene que estar en $LIOHOME/dat/ECP/
   INTEGER, DIMENSION(128) :: ZlistECP !Z de atomos con ECP
   LOGICAL :: cutECP !activa cuts en las integrales de ECP
   LIODBLE, PARAMETER :: cutecp2=7D2, cutecp3=7D2 !cuts limite para evitar 0 * NAN
   LIODBLE :: cut2_0, cut3_0 !valores de corte para las integrales de 2 y 3 centros (AAB y BAC)
   LOGICAL :: FOCK_ECP_read, FOCK_ECP_write !activan lectura y escritura Fock
   LOGICAL :: Fulltimer_ECP !activa los timers para int. radiales
   LIODBLE :: tlocal,tsemilocal,tQ1,tQ2,Tiempo, Taux !auxiiares para timers
   INTEGER, ALLOCATABLE, DIMENSION(:) :: ECPatoms_order 
   LOGICAL :: first_steep !control for multiple lio calls

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Dbug & Verbose Variables    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   LOGICAL :: ecp_debug,ecp_full_range_int
! ecp_debug activa el modo de debugueo, ecp_full_range_int activa el calculo de integrales radiales en todo el rango disponible por
! los arrays
   INTEGER :: local_nonlocal
! =1 solo calcula terminos locales <xi|V|xj>
! =2 solo calcula terminos no locales <xi|Ylm>V<Ylm|xj>
! default =0
   INTEGER :: verbose_ECP ! controla la impresion
   INTEGER :: inf_Q, NAN_Q, inf_Q2, NAN_Q2
   LOGICAL, dimension(-12:14) :: computedBD
   LOGICAL, DIMENSION (0:10,0:4,0:4) :: computedQnl1l2, usedQnl1l2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Effective Core potential Data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   INTEGER, DIMENSION(118) :: Zcore, Lmax
! Zcore(Z) carga del core para el ECP elegido del atomo con carga nuclear Z
! Lmax(Z) L maximo del ECP elegido para el atomo con carga nuclear Z
   INTEGER, DIMENSION(118,0:5) :: expnumbersECP
! expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y momento angular l del ECP
   INTEGER, DIMENSION(118,0:5,10) :: nECP
   LIODBLE, DIMENSION(118,0:5,10) :: bECP, aECP

! Los pseudopotenciales vienen dados por:

!                     _LMAX-1         _l
! \  / = \  /     +  \    \  /       \    |l,m\/l,m|
!  \/     \/LMAX     /_    \/l-LMAX  /_   |   /\   |
!                    l=0             m=-l


!donde Vl= Σ aECP * r^nECP * exp(-bECP r^2)
!          i

! Los  xECP estan escritos como: xECP(Z,l,i) Z carga del nucleo, l = momento angular de la expansion del ecp, i numero de funcion
! del ecp con Z,l

   INTEGER, DIMENSION (:), ALLOCATABLE :: IzECP ! cargas nucleares sin corregir por la carga del core (Zcore)
   INTEGER, DIMENSION (:,:), ALLOCATABLE :: Lxyz ! exponentes de la parte angular de la funcion de base i

!|xi> =Σci x^lx y^ly z^lz *e^(-a * r^2)
!      i

!j=1 lx, j=2 ly, j=3 lz para la funcion i de la base

   LIODBLE, DIMENSION(:), ALLOCATABLE :: VAAA, VAAB, VBAC,term1e
   LIODBLE, DIMENSION(:), ALLOCATABLE :: VAAB1, VBAC1
!        LIODBLE, DIMENSION(:,:), ALLOCATABLE :: VAAAcuadrada, VAABcuadrada, VBACcuadrada
! VXXX contiene los valores de la Matriz de Fock del pseudopotencial.
! VAAA integrales de un centro (base y ecp en el mismo atomo)
! VAAB integrales de 2 centros (1 base y ecp en el mismo atomo)
! VBAC integrales de 3 centros (ninguna base en el atomo con ecp)
! term1e contiene una copia de los terminos de 1e- sin la modificacion por agregar los terminos de los pseudopotenciales
   LIODBLE, DIMENSION(:,:,:,:), ALLOCATABLE :: dVAABcuadrada, dVBACcuadrada

   LIODBLE, DIMENSION(:,:,:), ALLOCATABLE :: dHcore_AAB, dHcore_ABC
! dHcore_AAB derivadas de VAAB (position(i,j),atomo,xyz)
! dHcore_ABC derivadas de VBAC (position(i,j),atomo,xyz)



   LIODBLE, DIMENSION(:,:), ALLOCATABLE :: distx, disty, distz
!guarda la distancia en x, y, z entre los atomos i y j dist(i,j)=xi-xj
!Cuidado, esta en unidades atomicas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Normalized Basis Coeficients    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   LIODBLE, DIMENSION(:,:), ALLOCATABLE :: Cnorm
! Lio guarda en el array c los coeficientes de la base |xi> = ci x^lx y^ly z^lz *e^(-a * r^2)
! En el caso de las funciones d el coeficiente ci esta normalizado para todos los terminos como si fueran xy, xz o yz.
! Para los terminos xx, yy y zz hay que dividir por 3^0.5 en los calculos de lio este factor ya se considera, mientras que en los
! calculos con pseudo-potenciales se modifica la base copiandola a Cnorm para hacer a la rutina mas facil de adaptar a otras 
! implementaciones
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Parameters For Radial Integration    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   INTEGER, PARAMETER, DIMENSION (6,6) :: alpha = reshape((/1,0,1,0,1, &
0,0,-3,0,-10,0,-21,0,0,15,0,105,0,0,0,0,-105,0, &
-1260,0,0,0,0,945,0,0,0,0,0,0,-10395/),(/6,6/))
   INTEGER, PARAMETER, DIMENSION (7,7) :: betha = reshape((/1,0,1,0,1, &
0,1,0,-1,0,-6,0,-15,0,0,0,3,0,45,0,210,0,0,0,-15,0,-420,0,0,0,0,0,105,0,4725,0,0,0,0,0,-945,0,0,0,0,0,0,0,10395/),(/7,7/))



! alpha y betha contains coeficients for expantion of modified spherical bessel function of the first kind Mk(x) in terms of 
! sinh(x)/x^i (betha(k+1,i)) and cosh(x)/x^i (alpha(k,i)) for k between 0 and 4, it is enought for energy calculations of functions
! s to g

   LIODBLE, DIMENSION (-12:14) :: Bn1, Bn2, Cn1, Cn2, rho, tau, sigma, sigmaR
   LIODBLE, DIMENSION (-12:14) :: Bn,Cn
!                                ͚ 
! Bni(j) contain the value of  ʃ t^(n-2-j)*(exp)-c(t-ai)^2 + exp(-c(t+ai)^2) dt
!                              ̊ 

!                                ͚ 
! Cni(j) contain the value of  ʃ t^(n-2-j)*(exp)-c(t-ai)^2 - exp(-c(t+ai)^2) dt
!                              ̊ 

!            ͚  
! rho(n) = ʃ exp(-cr^2) * sinh(Ka*r)* sinh(Kb*r) r^n dr
!          ̊ 
!              ͚ 
! sigma(n) = ʃ exp(-cr^2) * sinh(Ka*r)* cosh(Kb*r) r^n dr
!            ̊ 
!               ͚ 
! sigmaR(n) = ʃ exp(-cr^2) * cosh(Ka*r)* sinh(Kb*r) r^n dr
!             ̊ 
!            ͚ 
! tau(n) = ʃ exp(-cr^2) * cosh(Ka*r)* cosh(Kb*r) r^n dr
!          ̊ 

   LIODBLE, DIMENSION (0:10,0:6,0:6) :: Qnl1l2
!                     ͚
! Qnl1l2(n,l1,l2) = ʃ Ml1(kA*r)* Ml2(kB*r)*r^n * exp(-cr^2) dr 
!                   ̊ 
! Ml(x) are modified spherical Bessel functions

   LIODBLE, DIMENSION (-1:11,-1:5) :: Qnl
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!            ̊ 

  LIODBLE, ALLOCATABLE, DIMENSION (:,:,:,:,:) :: ECP_Ang_stack !temporary array for angular integrals

! Parameters
   LIODBLE, PARAMETER :: pi=3.14159265358979312D0, pi12=1.77245385090552D0 !pi12 = pi^0.5
! factorial
   INTEGER, DIMENSION(0:15) :: fac=(/1,1,2,6,24,120,720,5040,40320, 362880,3628800,39916800,479001600,1932053504,1278945280,  &
    2004310016/)
! double factorial
   INTEGER, PARAMETER :: i16 = selected_int_kind (16)
   INTEGER (KIND=i16), PARAMETER, DIMENSION(-1:33) :: doublefactorial=(/1_i16, 1_i16,1_i16,2_i16,3_i16,8_i16,15_i16,48_i16,   &
   105_i16,384_i16,945_i16,3840_i16,10395_i16,46080_i16,135135_i16,645120_i16,2027025_i16,10321920_i16,34459425_i16,          &
   185794560_i16,654729075_i16,3715891200_i16,13749310575_i16,81749606400_i16,316234143225_i16,1961990553600_i16,             &
   7905853580625_i16,51011754393600_i16,213458046676875_i16,1428329123020800_i16,6190283353629375_i16, 42849873690624000_i16, &
   191898783962510625_i16,1371195958099968000_i16,6332659870762850625_i16/)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Parameters for Angular Integrals    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

! li contains the coeficient ul,m,lx,ly,lz for expansion of normalized real spherical harmonics Slm in terms of unitary sphere
! polinomials

! S(l,m)= Σ ul(m,lx,ly,lz) * (x/r)^lx * (y/r)^ly * (z/r)^lz
!     lx+ly+lz=l

! ul=F(f(lx,ly,lz),m)

   LIODBLE, DIMENSION (1) :: l0 = (/0.5d0/pi12/)
   LIODBLE, PARAMETER :: aux1=sqrt(3.d0)/(2d0*pi12)
   LIODBLE, DIMENSION (3,-1:1) :: l1=reshape((/0.d0,aux1,0.d0,aux1,0.d0,0.d0,0.d0,0.d0,aux1/),(/3,3/))
   LIODBLE, PARAMETER :: aux2= 0.5d0 * sqrt(15.d0/pi)
   LIODBLE, PARAMETER :: aux3= 0.25d0 * sqrt(5.d0/pi)
   LIODBLE, DIMENSION (6,-2:2) :: l2=reshape((/0.d0,0.d0,0.d0,0.d0,aux2,0.d0,0.d0,aux2,0.d0,0.d0,0.d0,0.d0,2*aux3,   &
   0.d0,-aux3,0.d0,0.d0,-aux3,0.d0,0.d0,0.d0,aux2,0.d0,0.d0,0.d0,0.d0,-0.5d0*aux2,0.d0,0.d0,0.5d0*aux2/), (/6,5/))
   LIODBLE, PARAMETER :: aux4=0.25d0 * sqrt(17.5d0/pi)
   LIODBLE, PARAMETER :: aux5=0.5d0 * sqrt(105.d0/pi)
   LIODBLE, PARAMETER :: aux6=sqrt(10.5d0/pi)
   LIODBLE, PARAMETER :: aux7=sqrt(7.0d0/pi)
   LIODBLE, DIMENSION (10,-3:3) :: l3=reshape((/0.d0,0.d0,0.d0, -aux4,0.d0,0.d0,0.d0,0.d0,3.d0*aux4,0.d0,0.d0, 0.d0, &
   0.d0,0.d0,0.d0,aux5,0.d0,0.d0,0.d0,0.d0,0.d0,aux6,0.d0,-0.25d0*aux6,0.d0,0.d0,0.d0,0.d0,-0.25d0*aux6,0.d0,0.5d0*aux7,0.d0, &
   -0.75d0*aux7,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux7,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,aux6,0.d0,-0.25*aux6,0.d0,0.d0,-0.25d0*aux6,   &
   0.d0,0.d0,-0.5d0*aux5,0.d0,0.d0,0.d0,0.d0,0.5d0*aux5,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-3.d0*aux4,0.d0,0.d0,aux4/),  &
   (/10,7/))
   LIODBLE, PARAMETER :: aux8=sqrt(35.d0/pi) 
   LIODBLE, PARAMETER :: aux9=sqrt(35.d0/(2.d0*pi))
   LIODBLE, PARAMETER :: aux10=sqrt(5.d0/pi)
   LIODBLE, PARAMETER :: aux11=sqrt(5.d0/(2.d0*pi))
   LIODBLE, PARAMETER :: aux12=1.d0/pi12
   LIODBLE, DIMENSION (15,-4:4) :: l4=reshape((/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux8,0.d0,0.d0,0.d0, &
   0.d0,0.75d0*aux8,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux9,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,2.25d0*aux9,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
   0.d0,0.d0,0.d0,0.d0,4.5d0*aux10,0.d0,-0.75d0*aux10,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux10,0.d0,0.d0,3.d0*aux11,0.d0,            &
   -2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,1.5d0*aux12,0.d0,-4.5d0*aux12,0.d0,          &
   0.5625d0*aux12,0.d0,0.d0,0.d0,0.d0,-4.5d0*aux12,0.d0,1.125d0*aux12,0.d0,0.d0,0.5625d0*aux12,0.d0,0.d0,0.d0,0.d0,0.d0,      &
   3.d0*aux11,0.d0,-2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux11,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux10,0.d0,0.375d0*aux10,0.d0, &
   0.d0,0.d0,0.d0,2.25*aux10,0.d0,0.d0,0.d0,0.d0,-0.375*aux10,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-2.25d0*aux9,0.d0,0.d0,0.d0, &
   0.d0,0.75d0*aux9,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.1875d0*aux8,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-1.125d0*aux8,0.d0,0.d0,        &
   0.1875d0*aux8/),(/15,9/))

   LIODBLE, PARAMETER :: aux13=sqrt(154.d0/pi)
   LIODBLE, PARAMETER :: aux14=sqrt(385.d0/pi)
   LIODBLE, PARAMETER :: aux15=sqrt(770.d0/pi)
   LIODBLE, PARAMETER :: aux16=sqrt(1155.d0/pi)
   LIODBLE, PARAMETER :: aux17=sqrt(165.d0/pi)
   LIODBLE, PARAMETER :: aux18=sqrt(11.d0/pi)

   LIODBLE, DIMENSION (21,-5:5) :: l5=reshape((/ &
!m=-5
   0.d0,0.d0,0.d0,0.d0,0.d0,3d0/32d0 *aux13,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-15d0/16d0 * aux13, 0.d0,0.d0,0.d0,0.d0,  &
   15.d0/32.d0 * aux13,0.d0, &
!m=-4
   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.75d0*aux14,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux14,0.d0,0.d0,0.d0,0.d0, &
!m=-3
   0.d0,0.d0,0.d0,-0.25d0*aux15, 0.d0,1d0/32d0 * aux15, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, 0.75d0*aux15, 0.d0, -1d0/16d0 * aux15, &
   0.d0,0.d0,0.d0,0.d0,-3.d0/32.d0 * aux15, 0.d0, &
!m=-2
   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-0.5d0*aux16, 0.d0, 0.25d0*aux16, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, 0.25d0*aux16, 0.d0,    &
   0.d0,0.d0,0.d0, &
!m=-1
   0.d0,0.5d0*aux17, 0.d0, -0.75d0*aux17, 0.d0, 1d0/16d0 * aux17, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, -0.75*aux17, 0.d0,           &
   1d0/8d0 * aux17, 0.d0, 0.d0, 0.d0, 0.d0, 1d0/16d0 * aux17, 0.d0, &
!m=0
   0.5d0*aux18, 0.d0, -2.5d0*aux18, 0.d0, 15d0/16d0 * aux18, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-2.5d0*aux18, 0.d0,                &
   15d0/8d0 * aux18, 0.d0, 0.d0, 0.d0, 0.d0, 15d0/16d0 * aux18, 0.d0, 0.d0, &
!m=1
   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.5d0*aux17, 0.d0,-0.75d0*aux17, 0.d0, 1d0/16d0 * aux17, 0.d0, 0.d0, 0.d0, 0.d0,             &
   -0.75d0*aux17, 0.d0, 1d0/8d0 * aux17, 0.d0, 0.d0, 1d0/16d0 *aux17, &
!m=2
   0.d0, 0.d0, -0.25d0*aux16, 0.d0, 1d0/8d0 * aux16, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0d0, 0.25d0*aux16, 0.d0, 0.d0, 0.d0, 0.d0, &
   0.d0, 0d0, -1d0/8d0 * aux16, 0.d0, 0.d0, &
!m=3
   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-0.75d0*aux15, 0.d0, 3d0/32d0*aux15, 0.d0, 0.d0, 0.d0, 0.d0, 0.25d0*aux15, 0.d0,   &
   1d0/16d0*aux15, 0.d0, 0.d0, -1d0/32d0*aux15, &
!m=4
   0.d0,0.d0,0.d0,0.d0,3d0/16d0*aux14, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-9d0/8d0*aux14, 0.d0, 0.d0,0.d0,0.d0,          &
   3d0/16d0*aux14, 0.d0, 0.d0, &
!m=5
   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,15d0/32d0*aux13, 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-15d0/16d0*aux13, 0.d0,   &
   0.d0, 3d0/32d0*aux13 /),(/21,11/))




   LIODBLE, DIMENSION (0:10,0:10,0:10) :: angularint ! angularint(i,j,k)= ʃ(x/r)^i (y/r)^j (z/r)^k dΩ

CONTAINS

SUBROUTINE defineparams() ! define parametros de variables del modulo
   IMPLICIT NONE
   inf_Q=0; NAN_Q=0; inf_Q2=0; NAN_Q2=0
! parametros auxiliares para integrales angulares
   angularint=0.d0
   angularint(0,0,0)=12.5663706143592D0
   angularint(0,0,2)=4.18879020478639D0
   angularint(0,0,4)=2.51327412287183D0
   angularint(0,0,6)=1.79519580205131D0
   angularint(0,0,8)=1.39626340159546D0
   angularint(0,0,10)=1.14239732857811D0
   angularint(0,2,0)=4.18879020478639D0
   angularint(0,2,2)=0.837758040957278D0
   angularint(0,2,4)=0.359039160410262D0
   angularint(0,2,6)=0.199466200227923D0
   angularint(0,2,8)=0.126933036508679D0
   angularint(0,2,10)=8.787671758293127D-2
   angularint(0,4,0)=2.51327412287183D0
   angularint(0,4,2)=0.359039160410262D0
   angularint(0,4,4)=0.119679720136754D0
   angularint(0,4,6)=5.439987278943365D-2
   angularint(0,4,8)=2.929223919431043D-2
   angularint(0,4,10)=1.757534351658626D-2
   angularint(0,6,0)=1.79519580205131D0
   angularint(0,6,2)=0.199466200227923D0
   angularint(0,6,4)=5.439987278943364D-2
   angularint(0,6,6)=2.092302799593602D-2
   angularint(0,6,8)=9.764079731436807D-3
   angularint(0,6,10)=5.169218681348898D-3
   angularint(0,8,0)=1.39626340159546D0
   angularint(0,8,2)=0.126933036508679D0
   angularint(0,8,4)=2.929223919431043D-2
   angularint(0,8,6)=9.764079731436809D-3
   angularint(0,8,8)=4.020503418826921D-3
   angularint(0,8,10)=1.904448987865384D-3
   angularint(0,10,0)=1.14239732857811D0
   angularint(0,10,2)=8.787671758293127D-2
   angularint(0,10,4)=1.757534351658626D-2
   angularint(0,10,6)=5.169218681348898D-3
   angularint(0,10,8)=1.904448987865384D-3
   angularint(0,10,10)=8.161924233708786D-4
   angularint(2,0,0)=4.18879020478639D0
   angularint(2,0,2)=0.837758040957278D0
   angularint(2,0,4)=0.359039160410262D0
   angularint(2,0,6)=0.199466200227923D0
   angularint(2,0,8)=0.126933036508679D0
   angularint(2,0,10)=8.787671758293127D-2
   angularint(2,2,0)=0.837758040957278D0
   angularint(2,2,2)=0.119679720136754D0
   angularint(2,2,4)=3.989324004558467D-2
   angularint(2,2,6)=1.813329092981121D-2
   angularint(2,2,8)=9.764079731436809D-3
   angularint(2,2,10)=5.858447838862085D-3
   angularint(2,4,0)=0.359039160410262D0
   angularint(2,4,2)=3.989324004558467D-2
   angularint(2,4,4)=1.087997455788673D-2
   angularint(2,4,6)=4.184605599187203D-3
   angularint(2,4,8)=1.952815946287362D-3
   angularint(2,4,10)=1.033843736269780D-3
   angularint(2,6,0)=0.199466200227923D0
   angularint(2,6,2)=1.813329092981121D-2
   angularint(2,6,4)=4.184605599187203D-3
   angularint(2,6,6)=1.394868533062401D-3
   angularint(2,6,8)=5.743576312609886D-4
   angularint(2,6,10)=2.720641411236262D-4
   angularint(2,8,0)=0.126933036508679D0
   angularint(2,8,2)=9.764079731436809D-3
   angularint(2,8,4)=1.952815946287362D-3
   angularint(2,8,6)=5.743576312609887D-4
   angularint(2,8,8)=2.116054430961537D-4
   angularint(2,8,10)=9.068804704120875D-5
   angularint(2,10,0)=8.787671758293127D-2
   angularint(2,10,2)=5.858447838862085D-3
   angularint(2,10,4)=1.033843736269780D-3
   angularint(2,10,6)=2.720641411236262D-4
   angularint(2,10,8)=9.068804704120875D-5
   angularint(2,10,10)=3.548662710308168D-5
   angularint(4,0,0)=2.51327412287183D0
   angularint(4,0,2)=0.359039160410262D0
   angularint(4,0,4)=0.119679720136754D0
   angularint(4,0,6)=5.439987278943365D-2
   angularint(4,0,8)=2.929223919431043D-2
   angularint(4,0,10)=1.757534351658626D-2
   angularint(4,2,0)=0.359039160410262D0
   angularint(4,2,2)=3.989324004558467D-2
   angularint(4,2,4)=1.087997455788673D-2
   angularint(4,2,6)=4.184605599187203D-3
   angularint(4,2,8)=1.952815946287362D-3
   angularint(4,2,10)=1.033843736269780D-3
   angularint(4,4,0)=0.119679720136754D0
   angularint(4,4,2)=1.087997455788673D-2
   angularint(4,4,4)=2.510763359512322D-3
   angularint(4,4,6)=8.369211198374407D-4
   angularint(4,4,8)=3.446145787565932D-4
   angularint(4,4,10)=1.632384846741758D-4
   angularint(4,6,0)=5.439987278943365D-2
   angularint(4,6,2)=4.184605599187203D-3
   angularint(4,6,4)=8.369211198374408D-4
   angularint(4,6,6)=2.461532705404238D-4
   angularint(4,6,8)=9.068804704120875D-5
   angularint(4,6,10)=3.886630587480375D-5
   angularint(4,8,0)=2.929223919431043D-2
   angularint(4,8,2)=1.952815946287362D-3
   angularint(4,8,4)=3.446145787565932D-4
   angularint(4,8,6)=9.068804704120875D-5
   angularint(4,8,8)=3.022934901373625D-5
   angularint(4,8,10)=1.182887570102723D-5
   angularint(4,10,0)=1.757534351658626D-2
   angularint(4,10,2)=1.033843736269780D-3
   angularint(4,10,4)=1.632384846741758D-4
   angularint(4,10,6)=3.886630587480375D-5
   angularint(4,10,8)=1.182887570102723D-5
   angularint(4,10,10)=4.258395252369803D-6
   angularint(6,0,0)=1.79519580205131D0
   angularint(6,0,2)=0.199466200227923D0
   angularint(6,0,4)=5.439987278943364D-2
   angularint(6,0,6)=2.092302799593602D-2
   angularint(6,0,8)=9.764079731436807D-3
   angularint(6,0,10)=5.169218681348898D-3
   angularint(6,2,0)=0.199466200227923D0
   angularint(6,2,2)=1.813329092981121D-2
   angularint(6,2,4)=4.184605599187203D-3
   angularint(6,2,6)=1.394868533062401D-3
   angularint(6,2,8)=5.743576312609886D-4
   angularint(6,2,10)=2.720641411236262D-4
   angularint(6,4,0)=5.439987278943364D-2
   angularint(6,4,2)=4.184605599187203D-3
   angularint(6,4,4)=8.369211198374406D-4
   angularint(6,4,6)=2.461532705404237D-4
   angularint(6,4,8)=9.068804704120873D-5
   angularint(6,4,10)=3.886630587480374D-5
   angularint(6,6,0)=2.092302799593602D-2
   angularint(6,6,2)=1.394868533062401D-3
   angularint(6,6,4)=2.461532705404238D-4
   angularint(6,6,6)=6.477717645800625D-5
   angularint(6,6,8)=2.159239215266875D-5
   angularint(6,6,10)=8.449196929305162D-6
   angularint(6,8,0)=9.764079731436807D-3
   angularint(6,8,2)=5.743576312609886D-4
   angularint(6,8,4)=9.068804704120873D-5
   angularint(6,8,6)=2.159239215266875D-5
   angularint(6,8,8)=6.571597611681793D-6
   angularint(6,8,10)=2.365775140205445D-6
   angularint(6,10,0)=5.169218681348898D-3
   angularint(6,10,2)=2.720641411236262D-4
   angularint(6,10,4)=3.886630587480374D-5
   angularint(6,10,6)=8.449196929305162D-6
   angularint(6,10,8)=2.365775140205445D-6
   angularint(6,10,10)=7.885917134018151D-7
   angularint(8,0,0)=1.39626340159546D0
   angularint(8,0,2)=0.126933036508679D0
   angularint(8,0,4)=2.929223919431043D-2
   angularint(8,0,6)=9.764079731436809D-3
   angularint(8,0,8)=4.020503418826921D-3
   angularint(8,0,10)=1.904448987865384D-3
   angularint(8,2,0)=0.126933036508679D0
   angularint(8,2,2)=9.764079731436809D-3
   angularint(8,2,4)=1.952815946287362D-3
   angularint(8,2,6)=5.743576312609887D-4
   angularint(8,2,8)=2.116054430961537D-4
   angularint(8,2,10)=9.068804704120875D-5
   angularint(8,4,0)=2.929223919431043D-2
   angularint(8,4,2)=1.952815946287362D-3
   angularint(8,4,4)=3.446145787565932D-4
   angularint(8,4,6)=9.068804704120875D-5
   angularint(8,4,8)=3.022934901373625D-5
   angularint(8,4,10)=1.182887570102723D-5
   angularint(8,6,0)=9.764079731436809D-3
   angularint(8,6,2)=5.743576312609887D-4
   angularint(8,6,4)=9.068804704120875D-5
   angularint(8,6,6)=2.159239215266875D-5
   angularint(8,6,8)=6.571597611681793D-6
   angularint(8,6,10)=2.365775140205446D-6
   angularint(8,8,0)=4.020503418826921D-3
   angularint(8,8,2)=2.116054430961537D-4
   angularint(8,8,4)=3.022934901373625D-5
   angularint(8,8,6)=6.571597611681793D-6
   angularint(8,8,8)=1.840047331270902D-6
   angularint(8,8,10)=6.133491104236341D-7
   angularint(8,10,0)=1.904448987865384D-3
   angularint(8,10,2)=9.068804704120875D-5
   angularint(8,10,4)=1.182887570102723D-5
   angularint(8,10,6)=2.365775140205446D-6
   angularint(8,10,8)=6.133491104236341D-7
   angularint(8,10,10)=1.903497239245761D-7
   angularint(10,0,0)=1.14239732857811D0
   angularint(10,0,2)=8.787671758293127D-2
   angularint(10,0,4)=1.757534351658626D-2
   angularint(10,0,6)=5.169218681348898D-3
   angularint(10,0,8)=1.904448987865384D-3
   angularint(10,0,10)=8.161924233708786D-4
   angularint(10,2,0)=8.787671758293127D-2
   angularint(10,2,2)=5.858447838862085D-3
   angularint(10,2,4)=1.033843736269780D-3
   angularint(10,2,6)=2.720641411236262D-4
   angularint(10,2,8)=9.068804704120875D-5
   angularint(10,2,10)=3.548662710308168D-5
   angularint(10,4,0)=1.757534351658626D-2
   angularint(10,4,2)=1.033843736269780D-3
   angularint(10,4,4)=1.632384846741758D-4
   angularint(10,4,6)=3.886630587480375D-5
   angularint(10,4,8)=1.182887570102723D-5
   angularint(10,4,10)=4.258395252369803D-6
   angularint(10,6,0)=5.169218681348898D-3
   angularint(10,6,2)=2.720641411236262D-4
   angularint(10,6,4)=3.886630587480374D-5
   angularint(10,6,6)=8.449196929305162D-6
   angularint(10,6,8)=2.365775140205445D-6
   angularint(10,6,10)=7.885917134018151D-7
   angularint(10,8,0)=1.904448987865384D-3
   angularint(10,8,2)=9.068804704120875D-5
   angularint(10,8,4)=1.182887570102723D-5
   angularint(10,8,6)=2.365775140205446D-6
   angularint(10,8,8)=6.133491104236341D-7
   angularint(10,8,10)=1.903497239245761D-7
   angularint(10,10,0)=8.161924233708786D-4
   angularint(10,10,2)=3.548662710308168D-5
   angularint(10,10,4)=4.258395252369801D-6
   angularint(10,10,6)=7.885917134018151D-7
   angularint(10,10,8)=1.903497239245761D-7
   angularint(10,10,10)=5.526282307487692D-8
END SUBROUTINE defineparams


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Other Functions/Subs    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE asignacion(Z,simb)
! Toma el numero atomico (Z)  y devuelve el simbolo del elemento en mayusculas (simb)
   IMPLICIT NONE
   INTEGER,INTENT(in) :: Z
   CHARACTER, INTENT(out) :: simb*3
   CHARACTER (LEN=3), DIMENSION(118) :: vec

   vec=(/'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ',  &
   'CA ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ','GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ',  &
   'ZR ','NB ','MO ','TC ','RU ','RH ','PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ',  &
   'ND ','PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ','OS ','IR ','PT ','AU ',  &
   'HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ','PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ',  &
   'FM ','MD ','NO ','LR ','RF ','DB ','SG ','BH ','HS ','MT ','DS ','UUU','UUB','UUT','UUQ','UUP','UUH','UUS','UUO'/)

   simb=vec(Z)
   RETURN
END SUBROUTINE asignacion


LIODBLE FUNCTION NEXTCOEF(sgn,n,cados,expo,c0coef,coefn1, coefn2)
! calcula el coef Bn y Cn, con n<2
! si sgn=1 calcula Bn; si sgn=-1 calcula Cn
   INTEGER, INTENT(IN) :: sgn,n
   LIODBLE, INTENT(IN) :: cados,expo,c0coef,coefn1, coefn2
   if ( -n-1 .lt. 0) stop " se pide fac(n), n<0 en NEXTCOEF"
   NEXTCOEF=(1+sgn*(-1.d0)**(-n-1))*cados**(-n-1)*expo/fac(-n-1)-2.d0*c0coef*coefn2 +cados*coefn1
   NEXTCOEF=NEXTCOEF/(-n-1)
   RETURN
END FUNCTION NEXTCOEF



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Not Made By Nicolas    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



FUNCTION DAW(XX) RESULT(fn_val)
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!----------------------------------------------------------------------

! This function program evaluates Dawson's integral,

!                       2  / x   2
!                     -x   |    t
!             F(x) = e     |   e    dt
!                          |
!                          / 0

!   for a real argument x.

!   The calling sequence for this function is

!                   Y=DAW(X)

!   The main computation uses rational Chebyshev approximations published
!   in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek and Thacher.
!   This transportable program is patterned after the machine-dependent
!   FUNPACK program DDAW(X), but cannot match that version for efficiency or
!   accuracy.  This version uses rational approximations that are
!   theoretically accurate to about 19 significant decimal digits.
!   The accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-dependent
!   constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   XINF   = largest positive machine number
!   XMIN   = the smallest positive machine number.
!   EPS    = smallest positive number such that 1+eps > 1.
!            Approximately  beta**(-p), where beta is the machine radix
!            and p is the number of significant base-beta digits in a
!            floating-point number.

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XMAX   = absolute argument beyond which DAW(X) underflows.
!            XMAX = min(0.5/xmin, xinf).
!   XSMALL = absolute argument below DAW(X)  may be represented
!            by X.  We recommend XSMALL = sqrt(eps).
!   XLARGE = argument beyond which DAW(X) may be represented by
!            1/(2x).  We recommend XLARGE = 1/sqrt(eps).

!     Approximate values for some important machines are

!                        beta  p     eps     xmin       xinf

!  CDC 7600      (S.P.)    2  48  7.11E-15  3.14E-294  1.26E+322
!  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
!  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
!  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
!                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
!   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307

!                         XSMALL     XLARGE     XMAX

!  CDC 7600      (S.P.)  5.96E-08   1.68E+07  1.59E+293
!  CRAY-1        (S.P.)  5.96E-08   1.68E+07  5.65E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  2.44E-04   4.10E+03  4.25E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  1.05E-08   9.49E+07  2.24E+307
!  IBM 3033      (D.P.)  3.73D-09   2.68E+08  7.23E+75
!  VAX 11/780    (S.P.)  2.44E-04   4.10E+03  1.70E+38
!                (D.P.)  3.73E-09   2.68E+08  1.70E+38
!   (G Format)   (D.P.)  1.05E-08   9.49E+07  8.98E+307

!*******************************************************************

! Error Returns

!  The program returns 0.0 for |X| > XMAX.

! Intrinsic functions required are:

!     ABS


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!  Latest modification: March 9, 1992

!----------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

   REAL (dp), INTENT(IN)  :: xx
   REAL (dp)              :: fn_val

! Local variables

   INTEGER    :: i
   REAL (dp)  :: frac, sump, sumq, w2, x, y
!----------------------------------------------------------------------
!  Mathematical constants.
!----------------------------------------------------------------------
   REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         six25 = 6.25_dp, one225 = 12.25_dp, two5 = 25.0_dp
!----------------------------------------------------------------------
!  Machine-dependent constants
!----------------------------------------------------------------------
   REAL (dp), PARAMETER  :: XSMALL = 1.05D-08, XLARGE = 9.49D+07,   &
                         XMAX = 2.24D+307
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!----------------------------------------------------------------------
   REAL (dp), PARAMETER  :: P1(10) = (/  &
   -2.69020398788704782410D-12, 4.18572065374337710778D-10,  &
   -1.34848304455939419963D-08, 9.28264872583444852976D-07,  &
   -1.23877783329049120592D-05, 4.07205792429155826266D-04,  &
   -2.84388121441008500446D-03, 4.70139022887204722217D-02,  &
   -1.38868086253931995101D-01, 1.00000000000000000004D+00 /)
   REAL (dp), PARAMETER  :: Q1(10) = (/  &
   1.71257170854690554214D-10, 1.19266846372297253797D-08,  &
   4.32287827678631772231D-07, 1.03867633767414421898D-05,  &
   1.78910965284246249340D-04, 2.26061077235076703171D-03,  &
   2.07422774641447644725D-02, 1.32212955897210128811D-01,  &
   5.27798580412734677256D-01, 1.00000000000000000000D+00 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [2.5, 3.5)
!----------------------------------------------------------------------
   REAL (dp), PARAMETER  :: P2(10) = (/  &
   -1.70953804700855494930D+00,-3.79258977271042880786D+01,  &
   2.61935631268825992835D+01, 1.25808703738951251885D+01,  &
  -2.27571829525075891337D+01, 4.56604250725163310122D+00,  &
  -7.33080089896402870750D+00, 4.65842087940015295573D+01,  &
  -1.73717177843672791149D+01, 5.00260183622027967838D-01 /)
   REAL (dp), PARAMETER  :: Q2(9) = (/  &
   1.82180093313514478378D+00, 1.10067081034515532891D+03,  &
  -7.08465686676573000364D+00, 4.53642111102577727153D+02,  &
   4.06209742218935689922D+01, 3.02890110610122663923D+02,  &
   1.70641269745236227356D+02, 9.51190923960381458747D+02,  &
   2.06522691539642105009D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [3.5, 5.0]
!----------------------------------------------------------------------
   REAL (dp), PARAMETER  :: P3(10) = (/  &
   -4.55169503255094815112D+00,-1.86647123338493852582D+01,  &
   -7.36315669126830526754D+00,-6.68407240337696756838D+01,  &
   4.84507265081491452130D+01, 2.69790586735467649969D+01,  &
   -3.35044149820592449072D+01, 7.50964459838919612289D+00,  &
   -1.48432341823343965307D+00, 4.99999810924858824981D-01 /)
   REAL (dp), PARAMETER  :: Q3(9) = (/  &
   4.47820908025971749852D+01, 9.98607198039452081913D+01,  &
   1.40238373126149385228D+01, 3.48817758822286353588D+03,  &
  -9.18871385293215873406D+00, 1.24018500009917163023D+03,  &
  -6.88024952504512254535D+01,-2.31251575385145143070D+00,  &
   2.50041492369922381761D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  |x| > 5.0
!----------------------------------------------------------------------
   REAL (dp), PARAMETER  :: P4(10) = (/  &
   -8.11753647558432685797D+00,-3.84043882477454453430D+01,  &
   -2.23787669028751886675D+01,-2.88301992467056105854D+01,  &
   -5.99085540418222002197D+00,-1.13867365736066102577D+01,  &
   -6.52828727526980741590D+00,-4.50002293000355585708D+00,  &
   -2.50000000088955834952D+00, 5.00000000000000488400D-01 /)
   REAL (dp), PARAMETER  :: Q4(9) = (/  &
   2.69382300417238816428D+02, 5.04198958742465752861D+01,  &
   6.11539671480115846173D+01, 2.08210246935564547889D+02,  &
   1.97325365692316183531D+01,-1.22097010558934838708D+01,  &
   -6.99732735041547247161D+00,-2.49999970104184464568D+00,  &
   7.49999999999027092188D-01 /)
!----------------------------------------------------------------------

   x = xx
   IF (ABS(x) > xlarge) THEN
      IF (ABS(x) <= xmax) THEN
         fn_val = half / x
      ELSE
         fn_val = zero
      END IF
   ELSE IF (ABS(x) < xsmall) THEN
      fn_val = x
   ELSE
      y = x * x
      IF (y < six25) THEN
!----------------------------------------------------------------------
!  ABS(X) .LT. 2.5
!----------------------------------------------------------------------
         sump = p1(1)
         sumq = q1(1)
         DO  i = 2, 10
            sump = sump * y + p1(i)
            sumq = sumq * y + q1(i)
         END DO
         fn_val = x * sump / sumq
      ELSE IF (y < one225) THEN
!----------------------------------------------------------------------
!  2.5 .LE. ABS(X) .LT. 3.5
!----------------------------------------------------------------------
         frac = zero
         DO  i = 1, 9
            frac = q2(i) / (p2(i)+y+frac)
         END DO
            fn_val = (p2(10)+frac) / x
      ELSE IF (y < two5) THEN
!----------------------------------------------------------------------
!  3.5 .LE. ABS(X) .LT. 5.0
!---------------------------------------------------------------------
         frac = zero
         DO  i = 1, 9
            frac = q3(i) / (p3(i)+y+frac)
         END DO
            fn_val = (p3(10)+frac) / x
      ELSE
!----------------------------------------------------------------------
!  5.0 .LE. ABS(X) .LE. XLARGE
!------------------------------------------------------------------
            w2 = one / x / x
            frac = zero
            DO  i = 1, 9
               frac = q4(i) / (p4(i)+y+frac)
            END DO
            frac = p4(10) + frac
            fn_val = (half + half*w2*frac) / x
      END IF
   END IF
   Return
!---------- Last line of DAW ----------
END FUNCTION DAW


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

END MODULE ECP_mod
