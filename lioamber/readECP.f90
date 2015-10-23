!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Effective Core Potential Read Subroutines                          !
!                                                                    !
! This routines read ECP parameter located in                        !
! $LIOHOME/dat/ECP/"tipeECP"                                         !
!                                                                    ! 
! V 1.0 september 2015                                               !
!                                                                    !
! Nicolas Foglia                                                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%   ECP name format   %%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                    !
!Atom_symbol(in capitals)                                            !
!Junk    LMax_of_ECP   Z_CORE                                        !
!Junk                                                                !
!Number_of_functions_for_LMax                                        !
!n1_ECP b1_ECP a1_ECP                                                !
!n2_ECP b2_ECP a2_ECP                                                !
!...                                                                 !
!...                                                                 ! 
!...                                                                 !
!Number_of_functions_for_L=0                                         !
!n1_ECP b1_ECP a1_ECP                                                !
!n2_ECP b2_ECP a2_ECP                                                !
!...                                                                 !
!...                                                                 !
!...                                                                 !
!Number_of_functions_for_L=1                                         !
!...                                                                 !
!...                                                                 !
!...                                                                 !
!END                                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! where V = ai_ECP * r^ni_ECP * exp(-bi_ECP r^2)                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Example:                                                            !
!                                                                    !
!V     0                                                             !
!V-ECP     2     10                                                  !
!d-ul potential                                                      !
!  1                                                                 !
!1     12.7965900             -3.7914600                             !
!s-ul potential                                                      !
!  3                                                                 !
!0      1.5043900              3.4325800                             !
!2      4.2732100            178.0395700                             !
!2      3.8412100           -156.3060500                             !
!p-ul potential                                                      !
!  2                                                                 !
!0     30.8372000              3.7648100                             !
!2      8.3445400             58.9650300                             !
!CR     0                                                            !
!CR-ECP     2     10                                                 !
!d-ul potential                                                      !
!  1                                                                 !
!1     13.9330700             -3.6436200                             !
!s-ul potential                                                      !
!  3                                                                 !
!0      1.7213600              3.5967900                             !
!2      4.7568100            184.0766300                             !
!2      4.2534800           -161.4199600                             !
!p-ul potential                                                      !
!  2                                                                 !
!0     32.2822300              3.6886200                             !
!2      9.2813700             65.4804600                             !
!END                                                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



	subroutine LecturaECP
	use ECP_mod, only : ecptypes, tipeECP, ZlistECP,asignacion,Zcore, Lmax, expnumbersECP, nECP, bECP, aECP,verbose_ECP
	use garcha_mod, only : Iz, natom
	implicit none
	integer :: i,jj,ppotat,ppotati,ppotat_acum
	character :: simb*3

	write(*,*) 
	write(*,4160)
	write(*,4161)
	write(*,4162)
	write(*,4163)
	write(*,4164)

	ppotat=0
	ppotat_acum=0
	Lmax=0
	Zcore=0
	nECP=0
	bECP=0.d0
	aECP=0.d0

	do i=1,ecptypes
	   ppotati=0
	   call asignacion(ZlistECP(i),simb)
	   call dataECPelement(ZlistECP(i),simb)
	   do jj=1,natom
	      if (Iz(jj).eq.ZlistECP(i)) then
	         ppotati=ppotati+1
	      end if
	   end do
	   ppotat=ppotat+ppotati
	
	   if (verbose_ECP .gt. 0) then
	      write(*,4165) simb, ppotati,tipeECP
	   else
	      if (ppotati .gt. 0) write(*,4165) simb, ppotati,tipeECP
	   end if
	   ppotat_acum=ppotat_acum+ppotati
	   ppotati=0
	end do

	write(*,4166)
	write(*,4167),ppotat_acum
	write(*,4168)
	write(*,*)

	4160 format(4x,"╔═══════════════════════════════════════════════════╗")
	4161 format(4x,"║    READING EFFECTIVE CORE POTENTIAL PARAMETERS    ║")
	4162 format(4x,"╠══════════════╦═════════════════════╦══════════════╣")
	4163 format(4x,"║   ATOM TYPE  ║   NUMBER IN SYSTEM  ║   ECP TYPE   ║")
	4164 format(4x,"╠══════════════╬═════════════════════╬══════════════╣")
	4165 format(4x,"║",a9,5x,"║",i11,10x,"║",5x,a9 ,"║")
	4166 format(4x,"╠══════════════╬═════════════════════╬══════════════╝")
	4167 format(4x,"║     TOTAL    ║",i11,10x,"║")
	4168 format(4x,"╚══════════════╩═════════════════════╝")
	end subroutine lecturaECP


	Subroutine dataECPelement(Z,elemento)
!lee los valores del pseudopotencial para el "elemento" pedido
	use ECP_mod, only : Zcore, Lmax,expnumbersECP , nECP, bECP, aECP,tipeECP
	implicit none
	integer, intent(in) :: Z
!Z = nuclear charge
!tipeECP contiene el tipo de pseudopotencial
        character :: boba*6
!boba para eliminar partes de la lectura del archivo que no sirven
        character,intent(in) :: elemento*3!
!elemento es aquel que buscara en el archivo de pseudopotenciales
        character :: simbolo*3!
!simbolo sirve para buscar el atomo en el archivo ECP

!variables q contienen datos del pseudopotencial
!el pseudo potencial viene escrito como: Vl= A * r^n * e^(-B*r^2)
!nECP, bECP y aECP contiene los coeficientes n,b y A respectivamente de cada contraccion de los pseudopotenciales

!auxiliar variables
	integer :: u,l
	CHARACTER(len=255) :: lio
	logical :: existeECP,found

	CALL get_environment_variable("LIOHOME", lio)!
	existeECP=0
        inquire(file=TRIM(lio)//"/dat/ECP/"//tipeECP, exist=existeECP)

	if ( .not. existeECP) then
!cheque que el archivo ECP este
	   write(*,*) "Effective Core potential parameters file not found"
	   write(*,*) "check ",TRIM(lio),"/dat/ECP/"
	   stop
	else
	   open(unit=7,file=TRIM(lio)//"/dat/ECP/"//tipeECP)
	   found=.true.
	   do while(found)
	      read(7,*) simbolo
	      if (simbolo.eq.elemento) then
!lee linea por linea hasta encontrar el atomo que busca 
	         found=.false.       
	         read(7,*) boba, Lmax(Z), Zcore(Z)
! asigna Lmax y Ncore
!asigna los valores de nECP, bECP y aECP al leerlos desde el pseudopotencial
	         read(7,*)
	         read(7,*) expnumbersECP(Z,Lmax(Z))

	         do u=1, expnumbersECP(Z,Lmax(Z))
	            read(7,*) nECP(Z,Lmax(Z),u), bECP(Z,Lmax(Z),u), aECP(Z,Lmax(Z),u)
	         end do
!repite para l=0 hasta Lmax-1
	         do l=0, Lmax(Z)-1
	            read(7,*)
	            read(7,*) expnumbersECP(Z,l)
	            do u=1, expnumbersECP(Z,l)
	               read(7,*) nECP(Z,l,u), bECP(Z,l,u), aECP(Z,l,u)
	            end do
	         end do
	      elseif (simbolo .eq. "END") then
!corta el calculo si no encuentra el pseudopotencial
	         write (*,*) "element not found in ECP file"
	         write(*,*) "check ",TRIM(lio),"/dat/ECP/",tipeECP
	         stop
	      endif
	   enddo
	   close(7)
	endif
	end subroutine dataECPelement


	subroutine asignacion(Z,simb)
!Toma el numero atomico (Z)  y devuelve el simbolo del elemento en mayusculas (simb)
        implicit none
	integer,intent(in) :: Z
	character, intent(out) :: simb*3
	character (len=3), dimension(118) :: vec
	vec=(/'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'UUU', 'UUB', 'UUT', 'UUQ', 'UUP', 'UUH', 'UUS', 'UUO'/)
	simb=vec(Z)
	return
	end subroutine asignacion



