!Rutina para la lectura de coeficientes del ECP
	subroutine lecturaECP
	use ECP_mod, only : ecptypes, tipeECP, ZlistECP,asignacion,Zcore, Lmax, expnumbersECP, nECP, bECP, aECP
	use garcha_mod, only : Iz, natom
	implicit none
	integer :: i,jj,ppotat,ppotati
	character :: simb*3
	write(*,*) "Leyendo parametros de pseudopotencial"
	write(*,*) "tipos de atomos con ECP ", ecptypes
	write(*,*) "pseudopotencial elegido ", tipeECP
	ppotat=0
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
	write(*,*) ppotati, "atomo(s) de ", simb, "con ECP", tipeECP
                ppotati=0
        end do

       write(*,*) "atomo(s) con ecp", ppotat, tipeECP

	end subroutine lecturaECP


	Subroutine dataECPelement(Z,elemento)
!lee los valores de aECP, bECP, nECP, para el elemento pedido
	use ECP_mod, only : Zcore, Lmax,expnumbersECP , nECP, bECP, aECP,tipeECP
!el array data tiene la forma data(Z,i)
!i=-2 Zcore, i=-1 Lmaximo del ECP, i=0 terminos s del ECP, i=1 terminos p ....   i=5 terminos h
	implicit none
	integer, intent(in) :: Z
	integer :: u,l

!tipeECP contiene el tipo de pseudopotencial
!elemento es aquel que buscara en el archivo de pseudopotenciales, simbolo sirve para buscar el atomo el el archivo ECP
! y boba para eliminar partes de la lectura del archivo que no sirven
        character :: boba*6
        character,intent(in) :: elemento*3!
        character :: simbolo*3!

!variables q contienen datos del pseudopotencial
!el pseudo potencial viene escrito como: Vl,AREP= A * r^n * e^(-b*r^2)
!nECP, bECP y aECP contiene los coeficientes n,b y A respectivamente de cada contraccion de los pseudopotenciales

!declara variables extra
        logical :: existeECP,found!existe la uso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!! variables para lectura en carpeta de lio
        CHARACTER(len=255) :: lio!
            CALL get_environment_variable("LIOHOME", lio)!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	existeECP=0

        inquire(file=TRIM(lio)//"/libraries/ECP/"//tipeECP, exist=existeECP)

!cheque q el archivo ECP este, si no esta corta el programa y avisa q no esta el pseudopotencial
	if (.not.existeECP) then
		write(*,*) "falta el archivo de ECP o no tiene el nombre correcto"
		write(*,*) "revisar $LIOHOME/libraries/ECP/"
		stop

	else
!		open(unit=7,file=tECP)
		open(unit=7,file=TRIM(lio)//"/libraries/ECP/"//tipeECP)
!open abre el archivo en la unidad 7
		found=.true.
		do while(found)
			read(7,*) simbolo
!			write(*,*) simbolo,elemento
			if    (simbolo.eq.elemento) then
				write (*,*) simbolo
!lee linea por linea hasta encontrar el atomo que busca 
				found=.false.       
				read(7,*) boba, Lmax(Z), Zcore(Z)
				write(*,*) "carga",z,"lmax",Lmax(Z), "zcore", Zcore(Z)
! asigna Lmax y Ncore
!asigna los valores de nECP, bECP y aECP al leerlos desde el pseudopotencial

!LeeLmax
				read(7,*)
				read(7,*) expnumbersECP(Z,Lmax(Z))
				write(*,*) "terminos lmax",expnumbersECP(Z,Lmax(Z))
				do u=1, expnumbersECP(Z,Lmax(Z))
					read(7,*) nECP(Z,Lmax(Z),u), bECP(Z,Lmax(Z),u), aECP(Z,Lmax(Z),u)
!a					write(*,*) Z, dataECP(Z,dataECP(Z,-1)) , u
!					write(*,*) nECP(Z,dataECP(Z,dataECP(Z,-1)),u), bECP(Z,dataECP(Z,dataECP(Z,-1)),u), aECP(Z,dataECP(Z,dataECP(Z,-1)),u)
				end do
!repite para l=0 hasta Lmax-1
				do l=0, Lmax(Z)-1
					read(7,*)
					read(7,*) expnumbersECP(Z,l)
					do u=1, expnumbersECP(Z,l)
						read(7,*) nECP(Z,l,u), bECP(Z,l,u), aECP(Z,l,u)
!						write(*,*) Z, w, u
!						write(*,*) nECP(Z,w,u), bECP(Z,w,u), aECP(Z,w,u)
					end do
				end do
		        endif
!si no encuentra el peudopotencial explota, hay q agregar un try-except
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



