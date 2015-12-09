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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ECP name FORMAT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                                                                                  !
!Atom_symbol(in capitals)                                                                                                          !
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
!END                                                                                                                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! where V = ai_ECP * r^ni_ECP * exp(-bi_ECP r^2)                                                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Example:                                                                                                                          !
!                                                                                                                                  !
!V     0                                                                                                                           !
!V-ECP     2     10                                                                                                                !
!d-ul potential                                                                                                                    !
!  1                                                                                                                               !
!1     12.7965900             -3.7914600                                                                                           !
!s-ul potential                                                                                                                    !
!  3                                                                                                                               !
!0      1.5043900              3.4325800                                                                                           !
!2      4.2732100            178.0395700                                                                                           !
!2      3.8412100           -156.3060500                                                                                           !
!p-ul potential                                                                                                                    !
!  2                                                                                                                               !
!0     30.8372000              3.7648100                                                                                           !
!2      8.3445400             58.9650300                                                                                           !
!CR     0                                                                                                                          !
!CR-ECP     2     10                                                                                                               !
!d-ul potential                                                                                                                    !
!  1                                                                                                                               !
!1     13.9330700             -3.6436200                                                                                           !
!s-ul potential                                                                                                                    !
!  3                                                                                                                               !
!0      1.7213600              3.5967900                                                                                           !
!2      4.7568100            184.0766300                                                                                           !
!2      4.2534800           -161.4199600                                                                                           !
!p-ul potential                                                                                                                    !
!  2                                                                                                                               !
!0     32.2822300              3.6886200                                                                                           !
!2      9.2813700             65.4804600                                                                                           !
!END                                                                                                                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


      SUBROUTINE LecturaECP
       USE ECP_mod, ONLY : ecptypes, tipeECP, ZlistECP,asignacion,Zcore, Lmax, expnumbersECP, nECP, bECP, aECP,verbose_ECP
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
           ppotati=0
           CALL asignacion(ZlistECP(i),simb)
           CALL dataECPelement(ZlistECP(i),simb)
           DO jj=1,natom
              IF (Iz(jj).EQ.ZlistECP(i)) THEN
                 ppotati=ppotati+1
              END if
           END DO
           ppotat=ppotat+ppotati
        
           IF (verbose_ECP .GT. 0) THEN
              WRITE(*,4165) simb, ppotati,tipeECP
           ELSE
              IF (ppotati .GT. 0) WRITE(*,4165) simb, ppotati,tipeECP
           END if
           ppotat_acum=ppotat_acum+ppotati
           ppotati=0
        END DO

        WRITE(*,4166)
        WRITE(*,4167),ppotat_acum
        WRITE(*,4168)
        WRITE(*,*)

 4160 FORMAT(4x,"╔═════════════════════════════════════&
      ══════════════╗")
 4161 FORMAT(4x,"║    READING EFFECTIVE CORE POTENTIAL PARAMETERS    ║")
 4162 FORMAT(4x,"╠══════════════╦═════════════════════╦&
      ══════════════╣")
 4163 FORMAT(4x,"║   ATOM TYPE  ║   NUMBER IN SYSTEM  ║   ECP TYPE   ║")
 4164 FORMAT(4x,"╠══════════════╬═════════════════════╬&
      ══════════════╣")
 4165 FORMAT(4x,"║",a9,5x,"║",i11,10x,"║",5x,a9 ,"║")
 4166 FORMAT(4x,"╠══════════════╬═════════════════════╬&
      ══════════════╝")
 4167 FORMAT(4x,"║     TOTAL    ║",i11,10x,"║")
 4168 FORMAT(4x,"╚══════════════╩═════════════════════╝&
      ")
      END SUBROUTINE lecturaECP


      SUBROUTINE dataECPelement(Z,elemento)
!lee los valores del pseudopotencial para el "elemento" pedido
        USE ECP_mod, ONLY : Zcore, Lmax,expnumbersECP , nECP, bECP, aECP,tipeECP
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Z
!Z = nuclear charge
!tipeECP contiene el tipo de pseudopotencial
        CHARACTER :: boba*6
!boba para eliminar partes de la lectura del archivo que no sirven
        CHARACTER,INTENT(IN) :: elemento*3!
!elemento es aquel que buscara en el archivo de pseudopotenciales
        CHARACTER :: simbolo*3!
!simbolo sirve para buscar el atomo en el archivo ECP

!variables q contienen datos del pseudopotencial
!el pseudo potencial viene escrito como: Vl= A * r^n * e^(-B*r^2)
!nECP, bECP y aECP contiene los coeficientes n,b y A respectivamente de cada contraccion de los pseudopotenciales

!auxiliar variables
        INTEGER :: u,l
        CHARACTER(LEN=255) :: lio
        LOGICAL :: existeECP,found

        CALL get_environment_variable("LIOHOME", lio)!
        existeECP=.true.
        INQUIRE(FILE=TRIM(lio)//"/dat/ECP/"//tipeECP, EXIST=existeECP)


        IF ( .NOT. existeECP) THEN
!c.EQ.e que el archivo ECP este
         WRITE(*,*) "Effective Core potential parameters file.NOT.found"
         WRITE(*,*) "check ",TRIM(lio),"/dat/ECP/"
         STOP
        ELSE
         OPEN(UNIT=7,FILE=TRIM(lio)//"/dat/ECP/"//tipeECP)
         found=.true.
         DO WHILE(found)
          READ(7,*) simbolo
          IF (simbolo.EQ.elemento) THEN
!lee linea por linea hasta encontrar el atomo que busca 
           found=.false.       
           READ(7,*) boba, Lmax(Z), Zcore(Z)
! asigna Lmax y Ncore
!asigna los valores de nECP, bECP y aECP al leerlos desde el pseudopotencial
           READ(7,*)
           READ(7,*) expnumbersECP(Z,Lmax(Z))

           DO u=1, expnumbersECP(Z,Lmax(Z))
            READ(7,*) nECP(Z,Lmax(Z),u), bECP(Z,Lmax(Z),u), aECP(Z,Lmax(Z),u)
           END DO
!repite para l=0 hasta Lmax-1
           DO l=0, Lmax(Z)-1
            READ(7,*)
            READ(7,*) expnumbersECP(Z,l)
            DO u=1, expnumbersECP(Z,l)
             READ(7,*) nECP(Z,l,u), bECP(Z,l,u), aECP(Z,l,u)
            END DO
           END DO
          ELSEIF (simbolo .EQ. "END") THEN
!corta el calculo si no encuentra el pseudopotencial
           WRITE (*,*) "element ",elemento ,".NOT.found in ECP file"
           WRITE(*,*) "check ",TRIM(lio),"/dat/ECP/",tipeECP
           STOP
          ENDIF
         ENDDO
         CLOSE(7)
       ENDIF
      END SUBROUTINE dataECPelement


