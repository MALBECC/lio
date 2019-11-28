!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Effective Core Potential Read SUBROUTINEs                                                                                        !
!                                                                                                                                  !
! This routines READ ECP parameter located in $LIOHOME/dat/ECP/"tipeECP"                                                           !
!                                                                                                                                  ! 
! V 1.0 september 2015                                                                                                             !
!                                                                                                                                  !
! Nicolas Foglia                                                                                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ECP parameters FORMAT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                                                                                  !
!ATOM_SYMBOL(in capitals)                                                                                                          !
!Junk    LMax_of_ECP   Z_CORE                                                                                                      !
!Junk                                                                                                                              !
!Number_of_functions_for_LMax                                                                                                      !
!n1_ECP b1_ECP a1_ECP                                                                                                              !
!n2_ECP b2_ECP a2_ECP                                                                                                              !
!...                                                                                                                               !
!...                                                                                                                               ! 
!...                                                                                                                               !
!Number_of_functions_for_L=0                                                                                                       !
!n1_ECP b1_ECP a1_ECP                                                                                                              !
!n2_ECP b2_ECP a2_ECP                                                                                                              !
!...                                                                                                                               !
!...                                                                                                                               !
!...                                                                                                                               !
!Number_of_functions_for_L=1                                                                                                       !
!...                                                                                                                               !
!...                                                                                                                               !
!...                                                                                                                               !
!...                                                                                                                               !
!...                                                                                                                               !
!END                                                                                                                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                                                                                  !
!             LMax-1   l                                                                                                           !
! V = V_LMax + Σ       Σ |lm> V_l <lm|                                                                                             !
!             l=0     m=-l                                                                                                         !
!                                                                                                                                  !
!                                                                                                                                  !
! where Vl = Σ ai_ECP * r^ni_ECP * exp(-bi_ECP*r^2)                                                                                !
!            i                                                                                                                     !
!                                                                                                                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Example:                                                                                                                          !
!                                                                                                                                  !
!V                                                                                                                                 !
!V-ECP     2     10                                                                                                                !
!d-ul potential                                                                                                                    !
!  1                                                                                                                               !
!1     15.3000000             -2.5483000                                                                                           !
!s-ul potential                                                                                                                    !
!  3                                                                                                                               !
!0      1.8256900              3.5213800                                                                                           !
!2      4.9531400            225.0512300                                                                                           !
!2      2.9512300           -228.5468500                                                                                           !
!p-ul potential                                                                                                                    !
!  2                                                                                                                               !
!0     40.5262000              4.5328100                                                                                           !
!2      7.3458990             48.9660120                                                                                           !
!CR     0                                                                                                                          !
!CR-ECP     2     10                                                                                                               !
!d-ul potential                                                                                                                    !
!  1                                                                                                                               !
!1     23.8530452             -2.6648250                                                                                           !
!s-ul potential                                                                                                                    !
!  3                                                                                                                               !
!0      1.6345600              4.5679900                                                                                           !
!2      6.9382500            271.0586100                                                                                           !
!2      5.4215800           -127.6456800                                                                                           !
!p-ul potential                                                                                                                    !
!  2                                                                                                                               !
!0     41.6151600              4.2216200                                                                                           !
!2      8.2232100             55.6485300                                                                                           !
!END                                                                                                                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


      SUBROUTINE LecturaECP
       USE ECP_mod, ONLY : ecptypes, tipeECP, ZlistECP,asignacion,Zcore, Lmax, nECP, bECP, aECP,verbose_ECP
       USE garcha_mod, ONLY : Iz, natom
       IMPLICIT NONE
       INTEGER :: i,jj,ppotat,ppotati,ppotat_acum
       CHARACTER :: simb*3

       WRITE(*,*)
       WRITE(*,4160)
       WRITE(*,4161)
       WRITE(*,4162)
       WRITE(*,4163)
       WRITE(*,4164)

       ppotat=0
       ppotat_acum=0
       Lmax=0
       Zcore=0
       nECP=0
       bECP=0.d0
       aECP=0.d0

       DO i=1,ecptypes
          ppotati=0 !cantidad de atomos de un tipo con pseudopotencial
          CALL asignacion(ZlistECP(i),simb) !Toma el numero atomico (Z)  obtiene el simbolo del elemento en mayusculas (simb)
          CALL dataECPelement(ZlistECP(i),simb) !lee los valores del pseudopotencial para el "elemento" pedido
          DO jj=1,natom
             IF (Iz(jj).EQ.ZlistECP(i)) THEN
                ppotati=ppotati+1
             END IF
          END DO
          ppotat=ppotat+ppotati !cantidad de atomos toatales con pseudopotencial
       
          IF (verbose_ECP .GT. 0) THEN
             WRITE(*,4165) simb, ppotati,tipeECP
          ELSE
             IF (ppotati .GT. 0) WRITE(*,4165) simb, ppotati,tipeECP
          END if
          ppotat_acum=ppotat_acum+ppotati
          ppotati=0
       END DO
       WRITE(*,4166)
       WRITE(*,4167) ppotat_acum
       WRITE(*,4168)
       WRITE(*,*)

 4160 FORMAT(4x,"╔════════════════════════════════════", &
      "═══════════════╗")
 4161 FORMAT(4x,"║    READING EFFECTIVE CORE POTENTIAL PARAMETERS    ║")
 4162 FORMAT(4x,"╠══════════════╦═════════════════════", &
      "╦══════════════╣")
 4163 FORMAT(4x,"║   ATOM TYPE  ║   NUMBER IN SYSTEM  ║   ECP TYPE   ║")
 4164 FORMAT(4x,"╠══════════════╬═════════════════════", &
      "╬══════════════╣")
 4165 FORMAT(4x,"║",a9,5x,"║",i11,10x,"║",5x,a9 ,"║")
 4166 FORMAT(4x,"╠══════════════╬═════════════════════", &
      "╬══════════════╝")
 4167 FORMAT(4x,"║     TOTAL    ║",i11,10x,"║")
 4168 FORMAT(4x,"╚══════════════╩═════════════════════", &
      "╝")
      END SUBROUTINE lecturaECP


      SUBROUTINE dataECPelement(Z,elemento) !lee los valores del pseudopotencial para el "elemento" pedido
        USE ECP_mod, ONLY : Zcore, Lmax,expnumbersECP , nECP, bECP, aECP,tipeECP
! Zcore(Z) carga del core para el ECP elegido del atomo con carga nuclear Z
! Lmax(Z) L maximo del ECP elegido para el atomo con carga nuclear Z
! expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y momento angular l del ECP

!nECP, bECP y aECP contiene los coeficientes n,b y A respectivamente de cada contraccion de los pseudopotenciales
!                                                                                                                                  
!  Vl = Σ ai_ECP * r^ni_ECP * exp(-bi_ECP*r^2)                                                                                
!       i 

!tipeECP contiene el tipo de pseudopotencial

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Z !nuclear charge
        CHARACTER,INTENT(IN) :: elemento*3 !elemento es aquel que buscara en el archivo de pseudopotenciales
        CHARACTER :: simbolo*3!simbolo sirve para buscar el atomo en el archivo ECP

!auxiliar variables
        CHARACTER :: boba*6 
        INTEGER :: u,l
        CHARACTER(LEN=255) :: lio
        LOGICAL :: existeECP,found

        IF (tipeECP.EQ."USSERDEF") THEN
          existeECP=.true.
          INQUIRE(FILE="ECPPARAM", EXIST=existeECP)
        ELSE
          CALL get_environment_variable("LIOHOME", lio)!
          existeECP=.true.
          INQUIRE(FILE=TRIM(lio)//"/dat/ECP/"//tipeECP, EXIST=existeECP)
        ENDIF

        IF ( .NOT. existeECP) THEN !verifica la existencia del archivo de parametros
         WRITE(*,*) "Effective Core Potential parameters file not found"

         IF (tipeECP.EQ."USSERDEF") THEN
           WRITE(*,*) "Effective Core potential parameter file have", &
           " to be name ECPPARAM"
         ELSE       
           WRITE(*,*) "check ",TRIM(lio),"/dat/ECP/"
         ENDIF

         STOP

        ELSE
         IF (tipeECP.EQ."USSERDEF") THEN
           OPEN(UNIT=7,FILE="ECPPARAM")
         ELSE
           OPEN(UNIT=7,FILE=TRIM(lio)//"/dat/ECP/"//tipeECP)
         END IF

         found=.true.
         DO WHILE(found)
          READ(7,*) simbolo
          IF (simbolo.EQ.elemento) THEN !lee linea por linea hasta encontrar el atomo que busca 
           found=.false.       
           READ(7,*) boba, Lmax(Z), Zcore(Z) 
           READ(7,*)
           READ(7,*) expnumbersECP(Z,Lmax(Z))

           DO u=1, expnumbersECP(Z,Lmax(Z))
            READ(7,*) nECP(Z,Lmax(Z),u), bECP(Z,Lmax(Z),u), aECP(Z,Lmax(Z),u)
           END DO

           DO l=0, Lmax(Z)-1 !repite para l=0 hasta Lmax-1
            READ(7,*)
            READ(7,*) expnumbersECP(Z,l)
            DO u=1, expnumbersECP(Z,l)
             READ(7,*) nECP(Z,l,u), bECP(Z,l,u), aECP(Z,l,u)
            END DO
           END DO
          ELSEIF (simbolo .EQ. "END") THEN !corta el calculo si no encuentra el pseudopotencial
           WRITE (*,*) "element ",elemento ," not found in ECP file"
           IF (tipeECP.EQ."USSERDEF") THEN
             WRITE(*,*) "check your ECP param file"
           ELSE
             WRITE(*,*) "check ",TRIM(lio),"/dat/ECP/",tipeECP
           END IF
           STOP
          ENDIF
         ENDDO
         CLOSE(7)
       ENDIF
      END SUBROUTINE dataECPelement


