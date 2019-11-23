!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Effective Core Potential subroutines                                                                                             !
!                                                                                                                                  !
! This routines calculate Fock matrix elements for effective core potential (ECP).                                                 !
!                                                                                                                                  !
! V 2.00 November 2019 Refactor to new 'Pedronic' Lio structure, all integrals was moved to subm_intECP.f90                        !
! V 1.00 February 2016 Final version for ECP energy calculations                                                                   !
! V 0.95 december 2015 cutoff for ECP interacions, linear scaling                                                                  !  
! V 0.91 october 2015 optimized version for energy calculations                                                                    !
! V 0.9 september 2015 1st functional version for energy calculations                                                              !
!                                                                                                                                  !
! Nicolas Foglia                                                                                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! References                                                                                                                       !
! J. Chem. Phys. 65, 3826 (1976); http://dx.doi.org/10.1063/1.432900                                                               !
! J. Chem. Phys. 111, 8778 (1999); http://dx.doi.org/10.1063/1.480225                                                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Main ECP Routine    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      SUBROUTINE generalECP(tipodecalculo)
! Rutina principal, llama a las otras rutinas
! tipodecalculo=0 allocatea variables comunes y las lee de un restart
! tipodecalculo=1 alocatea variables y calcula terminos de un centro (AAA)
! tipodecalculo=2 calcula terminos de dos centros (ABB)
! tipodecalculo=3 calcula terminos de tres centros (ABC)
! tipodecalculo=4 desalocatea variables

       USE ECP_mod, ONLY : verbose_ECP,ecp_debug,Fulltimer_ECP,Tiempo,defineparams
       use faint_cpu, only: intECP
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: tipodecalculo
       DOUBLE PRECISION :: t1,t2

       IF (tipodecalculo .EQ. 0) THEN
        CALL READ_ECP() !lee la matriz de Fock de los pseudopotenciales de una restart
       ELSE IF (tipodecalculo .GE. 1 .AND. tipodecalculo .LE.3) THEN

        IF (tipodecalculo .EQ. 1) THEN
          CALL defineparams() !define parametros auxiliares
          CALL norm_C() !guarda los coeficientes de la base normalizados correctamente para las funciones d en un nuevo aray
          CALL obtainls() !obtiene una matriz con lx,ly y lz para cada base
          CALL WRITE_POST(3)
        ELSEIF (tipodecalculo .EQ. 2) THEN
          CALL obtaindistance() !obtiene arrays con la distancia en x,y y z entre cada par de atomos i,j
        ENDIF


        IF (Fulltimer_ECP) CALL cpu_time ( t1 )

        CALL intECP(tipodecalculo) !calcula terminos de un centro (AAA)

        IF (Fulltimer_ECP) THEN
         CALL cpu_time ( t2 )
         Tiempo = t2-t1
         CALL WRITE_POST(8)
        END IF

        CALL WRITE_POST(4+tipodecalculo)

       ELSEIF (tipodecalculo .EQ. 4) THEN
        CALL deallocateV() ! desalocatea variables de ECP

       ELSE
        CALL WRITE_POST(4)   
       ENDif

       IF (tipodecalculo.EQ.0 .OR. tipodecalculo.EQ.3) THEN
        IF ( verbose_ECP .GT. 0) THEN
         CALL WRITE_ECP_PARAMETERS()
        END IF
        IF ( verbose_ECP .GT. 1) THEN
         CALL WRITE_BASIS()
         CALL WRITE_ANG_EXP()
        END IF
        IF ( verbose_ECP .GT. 2) THEN
         CALL WRITE_DISTANCE()
         CALL WRITE_FOCK_ECP_TERMS()
         CALL WRITE_FOCK_ECP()
        END IF
       END IF

       IF ( ecp_debug .AND. (tipodecalculo.EQ.0 .OR. tipodecalculo.EQ.3)) THEN
        CALL WRITE_POST(2)
        CALL SEARCH_NAN()
       END IF
       END SUBROUTINE generalECP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Subr. for Prep.    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!rutinas para preparacion y calculo de variables previas


        SUBROUTINE allocate_ECP()
!alocatea todas las variables que van a necesitar los pseudopotenciales
        USE garcha_mod, ONLY : natom
        USE basis_data, ONLY : nshell
        USE ECP_mod, ONLY :VAAA,VAAB,VBAC,term1e,distx, disty, distz, IzECP,Lxyz,Cnorm,dVAABcuadrada, dVBACcuadrada, ECPatoms
!term1e terminos de fock de 1 electron sin agregarles VAAA
        IMPLICIT NONE
        INTEGER :: ns,np,nd,M,Mcuad
        ns=nshell(0)
        np=nshell(1)
        nd=nshell(2)
        M=ns+np+nd
        Mcuad=M*(M+1)/2
        ALLOCATE (VAAA(Mcuad),VAAB(Mcuad),VBAC(Mcuad),term1e(Mcuad))
        VAAA=0.d0
        VAAB=0.d0
        VBAC=0.d0
        term1e=0.d0
        ALLOCATE (distx(natom,natom), disty(natom,natom), distz(natom,natom))
        distx=0.d0
        disty=0.d0
        distz=0.d0
        ALLOCATE (IzECP(natom))
        ALLOCATE (Lxyz(M, 3))
        ALLOCATE (dVAABcuadrada(M,M,2,3), dVBACcuadrada(M,M,ECPatoms,3))
        END SUBROUTINE allocate_ECP


        SUBROUTINE deallocateV() !desalocatea variables de ECP
        USE ECP_mod, ONLY :VAAA,VAAB,VBAC,term1e,distx, disty, distz,IzECP,Lxyz,Cnorm, dVAABcuadrada, dVBACcuadrada
        IMPLICIT NONE
        DEALLOCATE (VAAA,VAAB,VBAC)
        DEALLOCATE (dVAABcuadrada, dVBACcuadrada)
        DEALLOCATE (term1e)
        DEALLOCATE (distx, disty, distz)
	DEALLOCATE (IzECP,Lxyz)
	DEALLOCATE (Cnorm)
        END SUBROUTINE deallocateV

	SUBROUTINE norm_C() !Escribe la matriz C corrigiendo la normalizacion de las funciones d
        USE basis_data, ONLY :nshell,c,ncont
	USE ECP_mod, ONLY :Cnorm
	IMPLICIT NONE
	INTEGER :: i,j,ns,np,nd
        DOUBLE PRECISION :: factor
        ns=nshell(0)
        np=nshell(1)
        nd=nshell(2)
	factor=sqrt(1.d0/3.d0)
	DO i=1,ns+np
	   DO j=1,ncont(i)
	      Cnorm(i,j)=c(i,j)
	   END DO
	END DO


!dx2,dxy,dyy,dzx,dzy,dzz
	DO i=ns+np+1,ns+np+nd,6
                DO j=1,ncont(i)
                        Cnorm(i,j)=c(i,j)*factor
                END DO
		DO j=1,ncont(i+1)
			Cnorm(i+1,j)=c(i+1,j)
		END DO
                DO j=1,ncont(i+2)
                        Cnorm(i+2,j)=c(i+2,j)*factor
                END DO
                DO j=1,ncont(i+3)
                        Cnorm(i+3,j)=c(i+3,j)
                END DO
                DO j=1,ncont(i+4)
                        Cnorm(i+4,j)=c(i+4,j)
                END DO
                DO j=1,ncont(i+5)
                        Cnorm(i+5,j)=c(i+5,j)*factor
                END DO
        END DO
	END SUBROUTINE norm_C

        SUBROUTINE ReasignZ()
!cambia la carga de los nucleos con pseudopotenciales sacandole la carga del core y guarda las cargas originales en IzECP
!tambien corrige la cantidad de electrones restando los que van al core
        USE garcha_mod, ONLY : Iz, natom, NCO
        USE basis_data, ONLY : nuc, nshell
        USE ECP_mod, ONLY : ZlistECP,IzECP,Zcore,ecptypes,asignacion
        IMPLICIT NONE
	CHARACTER  :: simb*3
        INTEGER :: i,j,elec_remov
	elec_remov=0
	
	WRITE(*,3900)
	WRITE(*,3901)
	WRITE(*,3902)
	WRITE(*,3903)
	WRITE(*,3904)
	WRITE(*,3905)
        DO i=1, natom !barre atomos
	 IzECP(i)=Iz(i)
         DO j=1, ecptypes !barre atomos con ecp
          IF (IzECP(i) .EQ. ZlistECP(j)) THEN
           Iz(i)=Iz(i)-Zcore(ZlistECP(j)) !cambia la carga del nucleo
           NCO=NCO-Zcore(ZlistECP(j))/2 !saca e-
	   CALL asignacion(IzECP(i),simb)
	   WRITE(*,3906) i,simb,Iz(i),Zcore(ZlistECP(j)) 
	   elec_remov=elec_remov+Zcore(ZlistECP(j))
          END IF
         END DO
        END DO

	WRITE(*,3907)
	WRITE(*,*)
	WRITE(*,3908)
	WRITE(*,3909) elec_remov
	WRITE(*,3910)


        3900 FORMAT(4x,"╔═════════════════════════════════", &
        "═════╗")
        3901 FORMAT(4x,"║          ATOM MODIFICATIONS          ║")
        3902 FORMAT(4x,"╠═══════╦══════╦═══════════╦══════", &
        "═════╣")
        3903 FORMAT(4x,"║ Atom  ║ Atom ║ Eff. Nuc. ║ Electrons ║")
        3904 FORMAT(4x,"║ Numb. ║ Type ║   Charge  ║  Removed  ║")
        3905 FORMAT(4x,"╠═══════╬══════╬═══════════╬══════", &
        "═════╣")
        3906 FORMAT(4x,"║",1x,i3,3x,"║",2x,a3,1x,"║",3x,i3,5x,"║",1x,i5,5x,"║")
        3907 FORMAT(4x,"╚═══════╩══════╩═══════════╩══════", &
        "═════╝")
        3908 FORMAT(4x,"╔═══════════════════════════╦═════", &
        "═════╗")
        3909 FORMAT(4x,"║  TOTAL ELECTRONS REMOVED  ║",1x,i5,4x,"║")
        3910 FORMAT(4x,"╚═══════════════════════════╩═════", &
        "═════╝")
        END SUBROUTINE ReasignZ


        SUBROUTINE obtainls()
!arma una matriz que contenga los exponentes de la parte angular de la base
! |x> = A x^lx y^ly z^lz *e^-ar^2
! angularL(i,j) j=1 => lx, j=2 => ly, j=3 => lz para la funcion i de la base
        USE basis_data, ONLY : nshell
        USE ECP_mod, ONLY : Lxyz
        IMPLICIT NONE
        INTEGER :: i,resto
        INTEGER :: M,lx,ly,lz
        M = nshell(0)+nshell(1)+nshell(2)
        Lxyz=0
        DO i=nshell(0)+1,M
           IF (i .LE. nshell(0)+nshell(1)) THEN !funciones p
              resto=modulo(i-nshell(0),3)
              IF (resto .EQ. 1) Lxyz(i,1)=1
              IF (resto .EQ. 2) Lxyz(i,2)=1
              IF (resto .EQ. 0) Lxyz(i,3)=1
           ELSE IF (i .le. M) THEN !funciones d
              resto=modulo(i-nshell(0)+nshell(1),6)
              IF (resto .EQ. 1) Lxyz(i,1)=2
              IF (resto .EQ. 2) THEN
                Lxyz(i,1)=1
                Lxyz(i,2)=1
              END IF
              IF (resto .EQ. 3) Lxyz(i,2)=2
              IF (resto .EQ. 4) THEN
                Lxyz(i,1)=1
                Lxyz(i,3)=1
              END IF
              IF (resto .EQ. 5) THEN
                Lxyz(i,2)=1
                Lxyz(i,3)=1
              END IF
              IF (resto .EQ. 0) Lxyz(i,3)=2
           END IF
        END DO
        END SUBROUTINE obtainls

        SUBROUTINE obtaindistance() !obtiene matrices de distancias
        USE garcha_mod, ONLY :r, natom
        USE ECP_mod, ONLY :distx, disty, distz !distx(i,j) = xi-xj
        IMPLICIT NONE
        INTEGER :: i,j
        DO i=1,natom
        DO j=1,natom
           distx(i,j)=r(i,1)-r(j,1)
           disty(i,j)=r(i,2)-r(j,2)
           distz(i,j)=r(i,3)-r(j,3)
        END DO
        END DO
        END SUBROUTINE obtaindistance



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Subr. for Write/Read Fock ECP    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	SUBROUTINE READ_ECP()
	USE basis_data, ONLY :nshell
	USE ECP_mod, ONLY : VAAA, VAAB, VBAC
	IMPLICIT NONE
	INTEGER :: ns,np,nd,k,M,MM
	logical :: hay_restart
	ns=nshell(0)
	np=nshell(1)
	nd=nshell(2)
	M=ns+np+nd
	MM=M*(M+1)/2

        INQUIRE(FILE="ECP_restart", EXIST=hay_restart)

        IF ( .NOT. hay_restart) THEN ! verifica que el archivo ECP este
         WRITE(*,*) "ECP_restart not Found"
	 STOP
        ELSE
	OPEN(UNIT=69,FILE="ECP_restart", STATUS='UNKNOWN', ACCESS='STREAM')
	DO k=1,MM
	   READ(69) VAAA(k),VAAB(k),VBAC(K)
	END DO
	CLOSE(69)
	ENDif
	END SUBROUTINE READ_ECP

	SUBROUTINE WRITE_ECP()
	USE basis_data, ONLY :nshell
	USE ECP_mod, ONLY : VAAA, VAAB, VBAC
	IMPLICIT NONE
	INTEGER :: ns,np,nd,k,M,MM
	ns=nshell(0)
	np=nshell(1)
	nd=nshell(2)
	M=ns+np+nd
	MM=M*(M+1)/2
	OPEN(UNIT=69,FILE="ECP_restart", STATUS='UNKNOWN', ACCESS='STREAM')
	DO k=1,MM
	   WRITE(69) VAAA(k),VAAB(k),VBAC(K)
	END DO
	CLOSE(69)
	END SUBROUTINE WRITE_ECP

	SUBROUTINE WRITE_POST(i)
	USE ECP_mod, ONLY :Tiempo
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: i
	WRITE(*,3911)
	IF (i.EQ.1) WRITE(*,3913)
	IF (i.EQ.2) WRITE(*,3914)
	IF (i.EQ.3) WRITE(*,3915)
	IF (i.EQ.4) WRITE(*,3916)
	IF (i.EQ.5) WRITE(*,3917)
	IF (i.EQ.6) WRITE(*,3918)
	IF (i.EQ.7) WRITE(*,3919)
	IF (i.EQ.8) WRITE(*,3920) Tiempo
	IF (i.EQ.9) WRITE(*,3921) Tiempo
	IF (i.EQ.10) WRITE(*,3922) Tiempo
	IF (i.EQ.11) WRITE(*,3923) Tiempo
	IF (i.EQ.12) WRITE(*,3924) Tiempo
        IF (i.EQ.13) WRITE(*,3925) Tiempo
	WRITE(*,3912) 


 3911   FORMAT(4x,"╔═══════════════════════════════════", &
        "═══╗")
 3912   FORMAT(4x,"╚═══════════════════════════════════", &
        "═══╝")
 3913   FORMAT(4x,"║    End of Effective Core Potential   ║")
 3914   FORMAT(4x,"║        Staring Debug Routines        ║")
 3915   FORMAT(4x,"║    Doing Eff. Core Pot. Integrals    ║")
 3916   FORMAT(4x,"║   ERROR in type of ECP calculation   ║")
 3917   FORMAT(4x,"║        1 Center Int. Finished        ║")
 3918   FORMAT(4x,"║       2 Centers Int. Finished        ║")
 3919   FORMAT(4x,"║       3 Centers Int. Finished        ║")
 3920   FORMAT(4x,"║         Time ",f12.9," s          ║")
 3921   FORMAT(4x,"║      Time Local ",f12.9," s       ║")
 3922   FORMAT(4x,"║    Time Semilocal ",f12.9," s     ║")
 3923   FORMAT(4x,"║       Radial Q1 ",f12.9," s       ║")
 3924   FORMAT(4x,"║       Radial Q2 ",f12.9," s       ║")
 3925   FORMAT(4x,"║    Ciclos Semiloc  ",f12.9," s    ║")
	END SUBROUTINE WRITE_POST

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Subr. for Verbose    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	SUBROUTINE WRITE_BASIS()
!luego agregar para q escriba el atomo, y el momento angular, y que USE los coeficiente d bien normalizados
	USE basis_data, ONLY : nshell,a,c,ncont
	IMPLICIT NONE
	INTEGER :: i,j,M
	M=nshell(0)+nshell(1)+nshell(2)
        WRITE(*,4010)
	WRITE(*,4011)
	WRITE(*,4012)
	WRITE(*,4013)
	WRITE(*,4014)
        DO i=1,M
                DO j=1,ncont(i)
                        WRITE(*,4015) i,j,a(i,j),c(i,j)
                END DO
        END DO
	WRITE(*,4016)

 4010 FORMAT(4x,"╔════════════════════════════════════", &
      "═══════╗")
 4011 FORMAT(4x,"║            NORMALIZED BASIS SET           ║")
 4012 FORMAT(4x,"╠═══════╦═══════╦════════════╦═══════", &
      "═══════╣")
 4013 FORMAT(4x,"║ Basis ║ Cont. ║  Exponent  ║ Coefficient  ║")
 4014 FORMAT(4x,"╠═══════╬═══════╬════════════╬═══════", &
      "═══════╣")
 4015 FORMAT(4x,"║",1x,i3,3x,"║",1x,i3,3x,"║",1x,f10.7,1x,"║",1x,f12.9,1x,"║")
 4016 FORMAT(4x,"╚═══════╩═══════╩════════════╩═══════", &
      "═══════╝ ")
	END SUBROUTINE WRITE_BASIS


	SUBROUTINE WRITE_FOCK_ECP_TERMS() !escribe los terminos de Fock del pseudopotencial
	USE ECP_mod, ONLY : VAAA, VAAB, VBAC
	USE basis_data, ONLY : nshell
	IMPLICIT NONE
	INTEGER :: i,j,k,M
	logical :: no_zero
	no_zero = .true.
	M=nshell(0)+nshell(1)+nshell(2)
	WRITE(*,4020)
        WRITE(*,4021)
        WRITE(*,4022)
        WRITE(*,4023)
        WRITE(*,4024)
	DO i=1,M
	DO j=1,M
	   IF (i .GE. j) THEN
	      k=i+(1-j)*(j-2*M)/2
	         IF (no_zero) THEN
	            IF ( abs(VAAA(k)+VAAB(k)+VBAC(k)) .gt. 1E-30) THEN
	               WRITE(*,4025) i,j,VAAA(k),VAAB(k),VBAC(k)
	            END IF
	         ELSE
	            WRITE(*,4025) i,j,VAAA(k),VAAB(k),VBAC(k)
	         END IF
            END IF
	END DO
	END DO
	WRITE(*,4026)

 4020 FORMAT(2x,"╔═══════════════════════════════════", &
      "════════════════════════════════════╗")
 4021 FORMAT(2x,"║                         FOCK Pseudopotencials         "&
                ,"                ║")
 4022 FORMAT(2x,"╠═════╦═════╦═══════════════════╦════", &
      "═══════════════╦═══════════════════╣")
 4023 FORMAT(2x,"║  i  ║  j  ║      <A|A|A>      ║      <A|A|B>      ║      <"&
         ,"B|A|C>      ║")
 4024 FORMAT(2x,"╠═════╬═════╬═══════════════════╬════", &
      "═══════════════╬═══════════════════╣")     
 4025 FORMAT(2x,"║",1x,i3,1x,"║",1x,i3,1x,"║",f18.15,1x,"║",f18.15,1x,"║",f18.15,1x,"║",f18.15,1x,"║")
 4026 FORMAT(2x,"╚═════╩═════╩═══════════════════╩════", &
      "═══════════════╩═══════════════════╝ ")
	END SUBROUTINE WRITE_FOCK_ECP_TERMS


	SUBROUTINE WRITE_FOCK_ECP()
	        USE ECP_mod, ONLY : VAAA, VAAB, VBAC
        USE basis_data, ONLY : nshell
        IMPLICIT NONE
        INTEGER :: i,j,k,M
        logical :: no_zero
        no_zero = .true.
        M=nshell(0)+nshell(1)+nshell(2)
        WRITE(*,4032)
        WRITE(*,4033)
        WRITE(*,4034)
        DO i=1,M
        DO j=1,M
           IF (i .GE. j) THEN
              k=i+(1-j)*(j-2*M)/2
                 IF (no_zero) THEN
                    IF ( abs(VAAA(k)+VAAB(k)+VBAC(k)) .GT. 1E-30) THEN
                       WRITE(*,4035) i,j,VAAA(k)+VAAB(k)+VBAC(k)
                    END if
                 ELSE
                    WRITE(*,4035) i,j,VAAA(k)+VAAB(k)+VBAC(k)
                 END IF
            END IF
        END DO
        END DO
        WRITE(*,4036)

	4032 FORMAT(5x,"╔═════╦═════╦═══════════════════════╗")
	4033 FORMAT(5x,"║  i  ║  j  ║ FOCK Pseudopotencials ║")
	4034 FORMAT(5x,"╠═════╬═════╬═══════════════════════╣")
	4035 FORMAT(5x,"║",1x,i3,1x,"║",1x,i3,1x,"║",2x,f18.15,3x,"║")
	4036 FORMAT(5x,"╚═════╩═════╩═══════════════════════╝ ")

	END SUBROUTINE WRITE_FOCK_ECP


	SUBROUTINE WRITE_ANG_EXP()
!escribe el exponente de la parte angular de la funcion de base:
!xi(x,y,z)= x^nx * y^ny * z^nz *f(r)
	USE ECP_mod, ONLY :lxyz
	USE basis_data, ONLY : nshell,nuc
	IMPLICIT NONE
	INTEGER :: i

	WRITE(*,4040)
	WRITE(*,4041)
	WRITE(*,4042)
        WRITE(*,4043)
        WRITE(*,4044)

	DO i=1,nshell(0)+nshell(1)+nshell(2)
	   WRITE(*,4045) i,nuc(i),lxyz(i,1),lxyz(i,2),lxyz(i,3)
	END DO

        WRITE(*,4046)

	4040 FORMAT(6x,"╔══════════════════════════╗")
	4041 FORMAT(6x,"║     Angular Exponents    ║")
	4042 FORMAT(6x,"╠═══════╦══════╦═══╦═══╦═══╣")
	4043 FORMAT(6x,"║ Basis ║ Atom ║ x ║ y ║ z ║")
	4044 FORMAT(6x,"╠═══════╬══════╬═══╬═══╬═══╣")
	4045 FORMAT(6x,"║",2x,i3,2x,"║",2x,i2,2x,"║",i2,1x,"║",i2,1x,"║",i2,1x,"║")
	4046 FORMAT(6x,"╚═══════╩══════╩═══╩═══╩═══╝ ")
	END SUBROUTINE WRITE_ANG_EXP

	SUBROUTINE WRITE_DISTANCE
	USE ECP_mod, ONLY :distx,disty,distz
	USE garcha_mod, ONLY : natom
	IMPLICIT NONE
	INTEGER :: i,j

	WRITE(*,4050)
	WRITE(*,4051)
	WRITE(*,4052)
	WRITE(*,4053)
	WRITE(*,4054)

	DO i=1,natom
	DO j=1,natom
	   WRITE(*,4055) i,j,distx(i,j),disty(i,j),distz(i,j)
	END DO
	END DO
	WRITE(*,4056)


 4050 FORMAT(2x,"╔══════════════════════════════════", &
      "════════════════════════════════════════", &
      "═══╗")
 4051 FORMAT(2x,"║                               Distances (Bohr) " &
                ,"                             ║")
 4052 FORMAT(2x,"╠════════╦════════╦════════════════", &
      "═══╦═══════════════════╦════════════════", &
      "═══╣")
 4053 FORMAT(2x,"║ atom i ║ atom j ║     distance x    ║     dista" &
                ,"nce y    ║     distance z    ║")
 4054 FORMAT(2x,"╠════════╬════════╬════════════════", &
      "═══╬═══════════════════╬════════════════", &
      "═══╣")
 4055 FORMAT(2x,"║",2x,i3,3x,"║",2x,i3,3x,"║",f18.14,1x,"║",f18.14,1x,"║",f18.14,1x,"║")
 4056 FORMAT(2x,"╚════════╩════════╩════════════════", &
      "═══╩═══════════════════╩════════════════", &
      "═══╝ ")

	END SUBROUTINE WRITE_DISTANCE

	SUBROUTINE WRITE_ECP_PARAMETERS()
	USE ECP_mod, ONLY : ecptypes,Lmax,ZlistECP,expnumbersECP,ZlistECP,nECP,bECP,aECP,asignacion,Zcore
	IMPLICIT NONE
	INTEGER :: k,l,w,z
        character :: simb*3
	WRITE(*,*)
	WRITE(*,*)
	WRITE(*,4060)
	WRITE(*,4061)
	WRITE(*,4062)

	DO k=1, ecptypes !barre todos los atomos con ECP
	z=ZlistECP(k)
	   CALL asignacion(Z,simb)
	   WRITE(*,*)
	   WRITE(*,4063)
	   WRITE(*,4064) simb,Zcore(z),Lmax(z)
	   WRITE(*,4065)

	   DO l=0,Lmax(z) !barre momento angular del ecp
	      WRITE(*,4066) l
	      WRITE(*,4067)
	      WRITE(*,4068)
	      WRITE(*,4069)

	         DO w =1, expnumbersECP(z,l) !barre funciones
	            WRITE(*,4070) necp(Z,l,w),bECP(z,L,w),aECP(z,L,w)
	         END DO
	      IF ( l .LT. Lmax(z) ) WRITE(*,4071)
	      IF ( l .EQ. Lmax(z) ) WRITE(*,4072)
	   END DO
	END DO

 4060 FORMAT(4x,"╔════════════════════════════════════", &
      "═════════╗")
 4061 FORMAT(4x,"║     EFFECTIVE CORE POTENTIAL PARAMETERS     ║")
 4062 FORMAT(4x,"╚════════════════════════════════════", &
      "═════════╝")
 4063 FORMAT(8x,"╔════════╦════════════════╦══════════", &
      "═╗")
 4064 FORMAT(8x,"║",2x,A3,3x,"║",2x,"ZCore=",i3,5x,"║",2x,"Lmax=",i1,3x,"║")
 4065 FORMAT(8x,"╠════════╩════════════════╩══════════", &
      "═╣")
 4066 FORMAT(8x,"║   L=",i1,"                               ║")
 4067 FORMAT(8x,"╠═══╦════════════════╦═══════════════", &
      "═╣")
 4068 FORMAT(8x,"║ n ║    exponent    ║  coefficient   ║")
 4069 FORMAT(8x,"╠═══╬════════════════╬═══════════════", &
      "═║")
 4070 FORMAT(8x,"║",1x,i1,1x,"║",f12.6,4x,"║",f12.6,4x,"║")
 4071 FORMAT(8x,"╠═══╩════════════════╩═══════════════", &
      "═╣")
 4072 FORMAT(8x,"╚═══╩════════════════╩═══════════════", &
      "═╝")

	END SUBROUTINE WRITE_ECP_PARAMETERS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Subr. for Debug    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	SUBROUTINE SEARCH_NAN()
	USE ECP_mod, ONLY : VAAA,VAAB,VBAC
	USE basis_data, ONLY : nshell
	IMPLICIT NONE
	INTEGER :: M,MM,i
	M=nshell(0)+nshell(1)+nshell(2)
	MM=M*(M+1)/2
	DO i=1,MM
	   IF (VAAA(i) .NE. VAAA(i)) STOP "NAN in VAAA"
	   IF (VAAB(i) .NE. VAAB(i)) STOP "NAN in VAAB"
	   IF (VBAC(i) .NE. VBAC(i)) STOP "NAN in VBAC"
	END DO
	END SUBROUTINE SEARCH_NAN

!	END SUBROUTINE generalECP
