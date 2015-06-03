	subroutine intECP(tipodecalculo)
	use garcha_mod, only :nshell,nuc
	use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP,nECP,bECP, aECP,Zcore, Lmax, expnumbersECP,VAAAcuadrada,lxyz, VAAA, VAAAcuadrada, VAAB, VAABcuadrada,pi
	implicit none

	integer, intent(in) :: tipodecalculo
!tipodecalculo=1 alocatea variables y alcula terminos AAA
!tipodecalculo=2 calcula terminos ABB y ABC
	integer z,l,t
	integer m1,lx,ly,lz
        integer :: ns,np,nd,M,i,j,e,pos
        ns=nshell(0)
        np=nshell(1)
        nd=nshell(2)
        M=ns+np+nd

	if (tipodecalculo .eq. 1) then
!prepara variables y calcula los terminos AAA	
        call ReasignZ()
!reasigna las cargas
        call obtainls()
!obtiene una matriz con lx,ly y lz
	call allocateV()
!allocatea la matriz de fock de pseudopotenciales, la matrix cuadrara para pruebas y el verctor con los terminos de 1 electron sin corregir
	call intECPAAA()
!calcula terminos AAA
	elseif (tipodecalculo .eq. 2) then
	call obtaindistance()
!ontiene arrays con la diatncia en x, y y z entre cada par de atomos i y j
	call intECPAAB()
!	write(90,*) VAAB
!calculo VAAB
	write(*,*) "en proceso rutinas VAAB"
	elseif (tipodecalculo .eq. 3) then
!calculo VABC
	write(*,*) "en proceso rutinas VABC"
	elseif (tipodecalculo .eq. 4) then
		call deallocateV()
        else
                Write(*,*) "ERROR in tipe of ECP calculation"
	endif



	if (tipodecalculo .ne. 4) then
!bloque de test de ECP
!Escribe coeficientes como testeo
	if ( .false. ) then
        do z=1,118
                do l=0, Lmax(Z)
                        do t=1, expnumbersECP(Z,l)
                                write(*,9018) Z,l,t,nECP(Z,l,t), bECP(Z,l,t), aECP(Z,l,t)
                        end do
                end do
        end do
	end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if ( .false. ) then
!escribe los coeficientes de expansir Slm en x^lx * y^ly * z^lz
		write (32,*) " l m lx ly lz U"
		do l=0,4
			do m1=-l,l
				do lx=0,l
				do ly=0,l-lx
				lz=l-lx-ly
					write (32,9012) l,m1,lx,ly,lz,Ucoef(l,m1,lx,ly,lz)
				end do
				end do
			end do
		end do
	end if



	if ( .false. ) then
!escribe coeficientes de fock del ECP AAA y testea que la matriz sea simetrica
	write(*,*) "test Vij AAA"
		do i=1,M
			do j=i,M
				if (nuc(i) .eq. nuc(j)) then
!					write(*,9015) i,j,VAAA(i,j), VAAA(j,i), VAAA(i,j)-VAAA(j,i)
					if (abs(VAAAcuadrada(i,j)-VAAAcuadrada(j,i)) .gt. 1E-10) then
						do e=1,20
						write(*,*) "*********  error matriz VAAA no simetrica  **********"
						write(*,*) VAAAcuadrada(i,j),VAAAcuadrada(j,i),VAAAcuadrada(i,j)-VAAAcuadrada(j,i)
						end do
					end if

				end if
			end do
		end do
!		write(*,*) VAAA
	end if

	if (tipodecalculo .eq. 2) then
	if ( .true. ) then
!escribe coeficientes de fock del ECP AAB y AAA
	write(*,*) "VAAB                     VAAA"
	do i=1,M
                do j=i,M
			write(*,*) i,j,VAABcuadrada(i,j),VAAAcuadrada(i,j)
		end do
	end do
	end if

        if ( .true. ) then
!testea que la matriz VAAB sea simetrica
        write(*,*) "test Vij AAB"
                do i=1,M
                        do j=i,M
                                if (nuc(i) .ne. nuc(j)) then
!                                       write(*,9015) i,j,VAAB(i,j), VAAB(j,i), VAAB(i,j)-VAAB(j,i)
                                        if (abs(VAABcuadrada(i,j)-VAABcuadrada(j,i)) .gt. 1E-10) then
                                                do e=1,20
						write(*,*) "*********  error matriz VAAB no simetrica  **********"
                                                write(*,*) VAABcuadrada(i,j),VAABcuadrada(j,i),VAABcuadrada(i,j)-VAABcuadrada(j,i)
                                                end do
                                        end if

                                end if
                        end do
                end do
        end if
	end if


	if ( .true. ) then
!compara la matriz cuadrada con el vector
	write(*,*) "testeando array de VAAA"
	        do i=1,M
                        do j=1,i
                                if (nuc(i) .eq. nuc(j)) then
					pos=i+(1-j)*(j-2*M)/2
!					write(*,9013) VAAA(pos), VAAAcuadrada(i,j), VAAA(pos)-VAAAcuadrada(i,j)
					if ( abs(VAAA(pos)-VAAAcuadrada(i,j)) .gt. 0.0000000000000001 ) then
						do e=1,20
						write(*,*) "no coinciden la matriz cuadrada con el vector VAAA",i,j,pos
						end do
					end if
				end if
			end do
		end do
	end if


	if (.false.) then
!escribe VAAA
	write(73,*) "testeando array de VAAA"
		do pos=1,M*(M+1)/2
			write(73,*) pos, VAAA(pos)
		end do
	end if



	if ( .true. ) then
!escribe matriz de esponentes de la parte angular
		write(*,*) "escribe lx,ly,lz"
		do i=1, M
			write(*,9014) i,lxyz(i,1),lxyz(i,2),lxyz(i,3)
		end do
	end if



	end if


!		l,m1,lx,ly,lz,Ucoef(l,m1,lx,ly,lz)
	9012 format(/1x,i2,2x,i2,2x,i2,2x,i2,2x,i2,2x,f18.10)
	9013 format(/1x,"vector",f18.10,2x,"matriz",f18.10,2x,"diff",f18.10)
	9014 format(/1x,"i",i3,2x,"lx",i2,2x,"ly",i2,2x,"lz",i2)
	9015 format(/1x,"i",i3,2x,"j",i3,2x,"Vij",f10.5,2x,"Vji",f10.5,2x,     "diff",f18.15)
        9018 format(/1x,'Z =',i4,2x, 'L =',i4,2x,'coefnumber =',i3,2x,  'n =',i2,2x, 'b =', f15.5,2x,'c =',f15.5)



	contains

	subroutine allocateV
!Allocatea VAAA, que contendra los terminos del ECP <A|A|A>
	use garcha_mod, only :nshell, natom
	use ECP_mod, only :VAAAcuadrada,VAABcuadrada, VBACcuadrada,VAAA,VAAB,VBAC,term1e,distx, disty, distz
!VAAA terminos de Fock de pseudopotencial para ECP y bases en el mismo atomo
!VAAA cuadrada, idem
!term1e terminos de fock de 1 electron sin agregarles VAAA
	implicit none
	integer :: ns,np,nd,M,Mcuad
	ns=nshell(0)
        np=nshell(1)
        nd=nshell(2)
        M=ns+np+nd
	allocate (VAAAcuadrada(M,M),VAABcuadrada(M,M), VBACcuadrada(M,M))
	VAAAcuadrada=0.d0
	VAABcuadrada=0.d0
	VBACcuadrada=0.d0
!VAAAcuadrada terminos de Fock de pseudopotencial para ECP y bases en el mismo atomo
        Mcuad=M*(M+1)/2
        allocate (VAAA(Mcuad),VAAB(Mcuad),VBAC(Mcuad),term1e(Mcuad))
	VAAA=0.d0
	VAAB=0.d0
	VBAC=0.d0
	term1e=0.d0
	allocate (distx(natom,natom), disty(natom,natom), distz(natom,natom))
	distx=0.d0
	disty=0.d0
	distz=0.d0
	end subroutine allocateV

	subroutine deallocateV
!desalocatea variables
	use ECP_mod, only :VAAAcuadrada,VAABcuadrada, VBACcuadrada,VAAA,VAAB,VBAC,term1e,distx, disty, distz
	implicit none
	deallocate (VAAAcuadrada,VAABcuadrada, VBACcuadrada, VAAA,VAAB,VBAC, term1e,distx, disty, distz)
	end subroutine deallocateV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!rutinas para VAAB

	subroutine intECPAAB()
!calcula los terminos de Fock para <xA|VA|xB> 1 base en el mismo atomo que el ecp

        use garcha_mod, only : nshell,nuc,a,c,ncont
!a(i,ni) exponente de la funcion de base i, contrccion ni
!c(i,ni) coeficiente de la funcion de base i, contrccion ni
!ncont(i) cantidad de contracciones de la funcion de base i
!nshell(i) cantidad de funciones i=1 s, i=2, p, i=3, d
        use ECP_mod, only : ecptypes,cutECP,IzECP,cutecp2,distx, disty, distz,Lxyz, VAAB, VAABcuadrada
!nECP, bECP, aECP valores del pseudo potencial
! aECP*r^b * exp(-bECP r^2)
!estan escritos como: xECP(Z,l,i) Z carga del nucleo, l del ecp, i numero de funcion del ecp con Z,l
!coeficientes(aECP) y exponentes(bECP) del pseudopotencial
!ecptypes cantidad de atomos con ECP
!IzECP cargas nucleares sin corregir por el Zcore
! Lmax(Z) L maximo del ECP
! Lxyz(i,j) contiene los exponentes de la parte angular de la funcion de base i
!|x> = A x^lx y^ly z^lz *e^-ar^2, j=1 lx, j=2, ly, j=3 lz para la funcion i de la base
!VAAA contiene los terminos <A|A|A> del pseudo potencial
!VAAAcuadrada es solo para testeo de simetria
        implicit none
        integer :: i,j,k,M,ii,ji,lxi,lxj,lyi,lyj,lzi,lzj
!lmaxQ, l   !orrar esta linea
!M cantidad total de funciones de base
        double precision :: Distcoef, AAB, acum, dx,dy,dz
	M=nshell(0)+nshell(1)+nshell(2)
	Distcoef=0.d0
	AAB=0.d0
	acum=0.d0
!	write(*,*) "******entre AAB"
        do i=1, M
!	write(*,*) i
!barre funciones de la base 
                do j=1, M
!	write(*,*) i,j,nuc(i),nuc(j)
!cambiar por do j=i,M para barrer solo la mitad de la matriz
!barre el otro coef de la base j>=i ya que la matriz tiene q ser simetrica
                        if (nuc(i) .ne. nuc(j)) then
!	write(*,*) nuc(i),nuc(j)
!agarra bases de atomos distintos
!				write(32,*) nuc(i),nuc(j),nuc(i)-nuc(j)
                                do k=1, ecptypes
!barre atomos con ecp
!	write(*,*) "*******antes del if"
					if (IzECP(nuc(i)) .eq. ZlistECP(k) .or. IzECP(nuc(j)) .eq. ZlistECP(k)) then
!solo calcula si el atomo tiene ECP
!	write(*,*) "************************************************pase un if"
						dx=distx(nuc(i),nuc(j))
						dy=disty(nuc(i),nuc(j))
						dz=distz(nuc(i),nuc(j))
!testedo con distancia a mano, nick
!						dx=-1.d0
!						dy=0.d0
!						dz=0.d0

						Distcoef=(dx**2 + dy**2 + dz**2)
!distancia entre atomos
!						write(*,*) "distancia", Distcoef
                                                lxi=Lxyz(i,1)
                                                lxj=Lxyz(j,1)
                                                lyi=Lxyz(i,2)
                                                lyj=Lxyz(j,2)
                                                lzi=Lxyz(i,3)
                                                lzj=Lxyz(j,3)
!exponentes de la parte angular
						if (IzECP(nuc(i)) .eq. ZlistECP(k)) then
!calculo para ECP en i, hay q repetir para j
                                                do ji=1, ncont(j)
	                                                if ( .not. cutECP .or. (Distcoef*a(j,ji) .lt. cutecp2)) then
!solo calcula los terminos que luego se multipliquen por un factor q no sea demasiado pequeño,
!el cut va a ser obligatorio, ya que si los atomos quedan muy separados el termino VAAB se obtiene como el producto de un numero MUY grande por uno que es casi 0
!y se obtienen NAN para el numero muy grande
!a distancias grades gana el termino que es casi 0

								acum=0.d0
		                                                do ii=1, ncont(i)
!ji y ii barren contracciones de las funcion de base
!me conviene que barra primero los terminos de la base de i, ya que el factor ezponencial que multiplica solo depende de j

						AAB=AABlocal(i,j,k,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz) +AABNonLocal(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)
                                                acum=acum+AAB*c(i,ii)
!multiplica por el coeficiente de la base
                		                                end do
                                                VAABcuadrada(i,j) = VAABcuadrada(i,j) + acum*c(j,ji)*4*pi*exp(-Distcoef*a(j,ji))
!multiplica por el otro coef de la base

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                if (i .ge. j) then
!este if hay q sacarlo al fina, cuando cambie el rango en el q barre j , en vez de comenzar en 1 comience en i
                                                        pos=i+(1-j)*(j-2*M)/2   !chekeada
                                                        VAAB(pos) = VAAB(pos) + + acum*c(j,ji)*4*pi*exp(-Distcoef*a(j,ji))!esta linea es lo unico que quedaria!!!!!!!!!!!!!
                                                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							end if
                                                end do
						end if

!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%arreglando

                                                if (IzECP(nuc(j)) .eq. ZlistECP(k)) then
!calculo para ECP en j
                                                do ii=1, ncont(i)
                                                        if ( .not. cutECP .or. (Distcoef*a(i,ii) .lt. cutecp2)) then
!solo calcula los terminos que luego se multipliquen por un factor q no sea demasiado pequeño,
                                                                acum=0.d0
                                                                do ji=1, ncont(j)
!ji y ii barren contracciones de las funcion de base
!me conviene que barra primero los terminos de la base de j, ya que el factor ezponencial que multiplica solo depende de j

                                                AAB=AABlocal(j,i,k,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz) +AABNonLocal(j,i,ji,ii,k,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
                                                acum=acum+AAB*c(j,ji)
!multiplica por el coeficiente de la base
                                                                end do
                                                VAABcuadrada(i,j) = VAABcuadrada(i,j) + acum*c(i,ii)*4*pi*exp(-Distcoef*a(i,ii))
!multiplica por el otro coef de la base

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                if (i .ge. j) then
!este if hay q sacarlo al fina, cuando cambie el rango en el q barre j , en vez de comenzar en 1 comience en i
                                                        pos=i+(1-j)*(j-2*M)/2   !chekeada
                                                        VAAB(pos) = VAAB(pos) + acum*c(i,ii)*4*pi*exp(-Distcoef*a(i,ii)) !esta linea es lo unico que quedaria!!!!!!!!!!!!!
                                                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                        end if
                                                end do
                                                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                        end if
                                end do

                        end if
                end do
        end do




!	write(91,*) VAAB




	end subroutine intECPAAB




        DOUBLE PRECISION function AABNonLocal(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
!calcula el termino no local del pseudopotencial
!pseudopotencial centrado en i
!los coef de la base se multiplican en la rutina que llama a esta

!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
!lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
        use garcha_mod, only : a,c
        use ECP_mod, only :ZlistECP,Lmax,aECP,nECP,bECP, expnumbersECP,Qnl
        implicit none
        integer, intent(in) :: i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj
	double precision, intent(in) :: dx,dy,dz
	double precision, dimension(3) :: Kvector
        integer :: l,m, term, lx,ly,lz, lambda,lmaxbase
!auxiliades para ciclos
        integer :: Z,n
        double precision :: A2, Acoef, acumang, acumint, AABx, AABy, AABz, Kmod,Ccoef
!Z= carga del nucleo
        AABNonLocal=0.d0
!        A2=0.d0


        Z=ZlistECP(k)
        n=lxi+lxj+lyi+lyj+lzi+lzj
	Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
        Kmod=2.d0 * sqrt(dx**2 + dy**2 + dz**2) *a(j,ji)
        lmaxbase=lxj+lyj+lzj
        do l = 0 , Lmax(z)-1
!barre todos los l de la parte no local
                do term=1, expnumbersECP(z,l)
!barre contracciones del ECP para el atomo con carga z y l del ecp
			Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
			Qnl=0.d0
			call Qtype1(Kmod,Ccoef,lmaxbase,necp(Z,l,term))
!			write(*,*) "nolocal",Kmod,Ccoef,lmaxbase,necp(Z,l,term)
!			write(*,*) "nolocal", Qnl(0,0),Qnl(1,0),Qnl(2,0),Qnl(3,0)
!			write(*,*) "variables Q",Kmod,Ccoef,lmaxbase,necp(Z,l,term)
!			write(*,*) "i,j,ai,aj,ci,cj",i,j,a(i,1),a(j,1),c(i,1), c(j,1)
			do lx=0, lxj
			do ly=0,lyj
			do lz=0,lzj
				acumint=0.d0
	                        do lambda=lxj+lyj+lzj,0,-2
					acumang=0.d0
	        	                do m=-l,l
						acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
!						write(*,*) Aintegral(l,m,lxi,lyi,lzi),OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
	        	                end do
					acumint=acumint+acumang*Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda) 
!					write(*,*) "Q",Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda),necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda
!					write(*,*) "acums",acumint,acumang
				end do
			AABz=AABz+acumint * dz**(lzj-lz) * comb(lzj,lz)
!			write(*,*) "AABz",comb(lzj,lz), dz**(lzj-lz), dz,lz,lzj
			acumint=0.d0
			end do
			AABy=AABy+AABz*dy**(lyj-ly)*comb(lyj,ly)
!			write(*,*) AABy,"AABy",comb(ly,lyj)
			AABz=0.d0
			end do
			AABx=AABx+AABy*dx**(lxj-lx)*comb(lxj,lx)
			AABy=0.d0
!			write(*,*) AABx,"AABx",comb(lx,lxj)
			end do
			AABNonLocal=AABNonLocal+AABx*aECP(z,L,term)
		end do
!                        Acoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
!Acoef es el exponente de la integral radial
!                        AABNonLocal=AABNonLocal+A2*aECP(z,L,term)*Q0(n+nECP(z,l,term),Acoef)
!Q0 integral radial   Q0(n,alpha)= int r^n exp(-alpha * r^2) dr fron 0 to inf
!aECP coeficiente del ECP
!                        A2=0.d0
        end do
	return
	end function AABNonLocal







        DOUBLE PRECISION function AABlocal(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz)

!Calcula el termino local del pseudopotencial
!   [<xi|V(LM)|xj>] =  para el pseudopotencial en i
!los coef N de la base se multiplican en la rutina que llama a esta

!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
        use garcha_mod, only : a,c
        use ECP_mod, only :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, Qnl, angularint
!expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y l del ECP
        implicit none
        integer :: w,n,z,l,lxi, lyi, lzi, lambda, lmaxbase
        integer, intent(in) :: i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi
        double precision :: Ccoef, acum,dx,dy,dz,integral, Kmod
	double precision, dimension (3) :: Kvector
        z=ZlistECP(k)
        L=Lmax(ZlistECP(k))
        lmaxbase=lx+ly+lz+kxi+kyi+kzi
        AABlocal=0.d0
	acum=0.d0
	integral=0.d0
	Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
	Kmod= 2.d0 * sqrt(dx**2 + dy**2 + dz**2) *a(j,ji)


!	write(*,*) "values", a(1,1), c(1,1)
!	write(*,*) "K", Kmod
        do w =1, expnumbersECP(z,l)
!barre todos los terminos del Lmaximo
	        Qnl=0.d0
		Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
!		WRITE(*,*) "cCOEF",bECP(z,L,w),a(i,ii),a(j,ji),i,ii,j,ji
!parece q esta mal K 
!en C el problema parece ser bECP
!		WRITE(*,*) "bECP",z,L,w
		call Qtype1(Kmod,Ccoef,lmaxbase,necp(Z,l,w))
!                write(*,*) "local",Kmod,Ccoef,lmaxbase,necp(Z,l,w)
!                write(*,*) "local", Qnl(0,0),Qnl(1,0),Qnl(2,0),Qnl(3,0)



		do lxi=0,lx
		do lyi=0,ly
		do lzi=0,lz
!a,b,c barren desde 0 hasta lx,ly,lz
			do lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
!hay q chekear este rando de lambda
				integral=integral + OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) * Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
!				write(*,*) "integral", integral, OMEGA1(Kvector,lambda,lxi,lyi,lzi), Qnl(lxi+lyi+lzi+ nECP(Z,l,w),lambda)
!				write(*,*) "Q", Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda),lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda
			end do
			acum= acum + integral*dx**(lx-lxi) * dy**(ly-lyi) * dz**(lz-lzi) *comb(lx,lxi) *comb(ly,lyi) * comb(lz,lzi)
			integral=0.d0
		end do
		end do
		end do
		AABlocal=AABlocal+aECP(z,L,w)*acum
		acum=0.d0
!aca multiplica x el coef
	end do
        return
        end function AABlocal


        double precision function comb(a,b)
!devuelve el combinatorio a,b; a>=b
!en el futuro puedo pasarlo a un array
	use ECP_mod, only :fac
        integer, intent(in) :: a,b
        comb=0.d0
        comb=fac(a)/((fac(b)*fac(a-b)))
        return
        end function comb





	subroutine obtaindistance()
!obtiene matrices de distancias
	use garcha_mod, only :r, natom
	use ECP_mod, only :distx, disty, distz
!distx(i,j) = xi-xj
	implicit none
	integer :: i,j
	do i=1,natom
	do j=1,natom
	distx(i,j)=r(i,1)-r(j,1)
	disty(i,j)=r(i,2)-r(j,2)
	distz(i,j)=r(i,3)-r(j,3)
	end do
	end do
	end subroutine obtaindistance


	Subroutine intECPAAA
!calcula los terminos de Fock para bases y pseudopotenciales en el mismo atomo
	use garcha_mod, only : a,c,ncont, nshell, nuc, ncont
!a(i,ni) exponente de la funcion de base i, contrccion ni
!c(i,ni) coeficiente de la funcion de base i, contrccion ni
!ncont(i) cantidad de contracciones de la funcion de base i
!nshell(i) cantidad de funciones i=1 s, i=2, p, i=3, d
	use ECP_mod, only :nECP,bECP, aECP, ecptypes, IzECP, Lmax, Lxyz, VAAAcuadrada,VAAA
!nECP, bECP, aECP valores del pseudo potencial
! aECP*r^b * exp(-bECP r^2)
!estan escritos como: xECP(Z,l,i) Z carga del nucleo, l del ecp, i numero de funcion del ecp con Z,l
!coeficientes(aECP) y exponentes(bECP) del pseudopotencial
!ecptypes cantidad de atomos con ECP
!IzECP cargas nucleares sin corregir por el Zcore
! Lmax(Z) L maximo del ECP
! Lxyz(i,j) contiene los exponentes de la parte angular de la funcion de base i
!|x> = A x^lx y^ly z^lz *e^-ar^2, j=1 lx, j=2, ly, j=3 lz para la funcion i de la base
!VAAA contiene los terminos <A|A|A> del pseudo potencial
!VAAAcuadrada es solo para testeo de simetria
	implicit none
	integer :: i,j,k, ii,ji, lxi, lxj,lyi,lyj,lzi, lzj,pos,M
!lmaxQ, l   !orrar esta linea
!M cantidad total de funciones de base
	double precision :: local, nonlocal, Kbessel, exponente, AAA, acum
	local=0.d0
	nonlocal=0.d0
	Kbessel=0.d0
	exponente=0.d0
!antes habia puesto esta variable, la saco pero dejo el comentario por si me olvide q se usaba en algun lado. sacarla al final!!!!!
!	lmaxQ=0
	AAA=0.d0
	M=nshell(0)+nshell(1)+nshell(2)


	do i=1, M
!barre funciones de la base 
		do j=1, M
!cambiar por do j=i,M para barrer solo la mitad de la matriz
!barre el otro coef de la base j>=i ya que la matriz tiene q ser simetrica
			if (nuc(i) .eq. nuc(j)) then
!solo calcula si los terminos corresponden al mismo atomo
				do k=1, ecptypes
!barre atomos con ecp
					if (IzECP(nuc(i)) .eq. ZlistECP(k))then
!solo calcula si el atomo tiene ECP
!						write(*,*) "1 coincidencia", ZlistECP(k)
						do ii=1, ncont(i)
						do ji=1, ncont(j)
!ji y ii barren contracciones de las funcion de base
						lxi=Lxyz(i,1)
						lxj=Lxyz(j,1)
						lyi=Lxyz(i,2)
						lyj=Lxyz(j,2)
						lzi=Lxyz(i,3)
						lzj=Lxyz(j,3)
!exponentes de la parte angular
						AAA=AAAlocal(i,j,k,ii,ji,lxi+lxj,lyi+lyj,lzi+lzj) + AAANonLocal(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj)
!suma los terminos locales y no locales del pseudopotencial
						acum=acum+AAA*c(j,ji)
!multiplica por el coeficiente de la base
						end do
						VAAAcuadrada(i,j) = VAAAcuadrada(i,j) + acum*c(i,ii)
!multiplica por el otro coef de la base

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						if (i .ge. j) then
!este if hay q sacarlo al fina, cuando cambie el rango en el q barre j , en vez de comenzar en 1 comience en i
							pos=i+(1-j)*(j-2*M)/2   !chekeada
!							write(31,*) pos
							VAAA(pos) = VAAA(pos) + acum*c(i,ii) !esta linea es lo unico que quedaria!!!!!!!!!!!!!
!							write(*,*) VAAAcuadrada(i,j),VAAA(pos),VAAA(pos)-VAAAcuadrada(i,j)
							if (abs(VAAA(pos)-VAAAcuadrada(i,j)) .gt. 0.000000000001) then
								do e=1,20
									write(*,*) i,j,pos,"arma mal el vector!!!!"
								end do
							end if
!							write(*,*) VAAA(pos),VAAAcuadrada(i,j)
						end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

						acum=0.d0
						end do
					end if
				end do

			end if
		end do
	end do

        end subroutine intECPAAA

	DOUBLE PRECISION function AAANonLocal(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj)
!calcula el termino no local del pseudopotencial

!suma m=-l hasta l  [<xi|lm> V(l-LM) <lm|xj>] = 
!=Ni Nj C(l-LM) suma m=-l, l int(r^n exp(-alpha *r^2) dr, 0, inf) * int ((x/r)^ni (y/r)^li (z/r)^mi Ylm d(angulo solido) *
! * int ((x/r)^nj (y/r)^lj (z/r)^mj Ylm d(angulo solido)
!los coef de la base se multiplican en la rutina que llama a esta

!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
!lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
	use garcha_mod, only : a,c
	use ECP_mod, only :ZlistECP,Lmax,aECP,nECP,bECP, expnumbersECP
	implicit none
	integer, intent(in) :: i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj
	integer :: l,m, term
!auxiliades para ciclos
	integer :: Z,n
	double precision :: A2, Acoef
!Z= carga del nucleo
	AAANonLocal=0.d0
	A2=0.d0
	Z=ZlistECP(k)
	n=lxi+lxj+lyi+lyj+lzi+lzj

	do l = 0 , Lmax(z)-1
!barre todos los l de la parte no local
		do term=1, expnumbersECP(z,l)
!barre contracciones del ECP para el atomo con carga z y l del ecp
			do m=-l,l
!barre m
				A2=A2+Aintegral(l,m,lxi,lyi,lzi)*Aintegral(l,m,lxj,lyj,lzj)
!A2 contiene la parte angular de la integral
			end do
			Acoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
!Acoef es el exponente de la integral radial
			AAANonLocal=AAANonLocal+A2*aECP(z,L,term)*Q0(n+nECP(z,l,term),Acoef)
!Q0 integral radial   Q0(n,alpha)= int r^n exp(-alpha * r^2) dr fron 0 to inf
!aECP coeficiente del ECP
			A2=0.d0
		end do
	end do
	return
	end function AAANonLocal

	DOUBLE PRECISION function Aintegral(l,m,lx,ly,lz)
!calcula int ((x/r)^lx (y/r)^ly (z/r)^lz * Ylm d(angulo solido)
!expande el armonico esferico en cartesianas e integra
!Ucoef es el coeficiente de la expansion
	use ECP_mod, only :angularint
	implicit none
	integer, intent(in) :: l,m,lx,ly,lz	
	integer :: ix,iy,iz
	Aintegral=0.d0
		do ix=0,l
			do iy=0,l-ix
				iz=l-ix-iy
				Aintegral=Aintegral+Ucoef(l,m,ix,iy,iz)*angularint(lx+ix,ly+iy,lz+iz)
!				write(*,*) "integral"
!angularint(nx,ny,nz) calcula int ((x/r)^nx (y/r)^ny (z/r)^nz d(angulo solido)
			end do
		end do
	return
	end function Aintegral



	DOUBLE PRECISION function AAAlocal(i,j,k,ii,ji,lx,ly,lz)
!Calcula el termino local del pseudopotencial
!   [<xi|V(LM)|xj>] = 
!=Ni Nj C(LM) int(r^n exp(-alpha *r^2) dr, 0, inf) * int ((x/r)^ni (y/r)^li (z/r)^mi d(angulo solido)
!los coef de la base se multiplican en la rutina que llama a esta

!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
        use garcha_mod, only : a,c
        use ECP_mod, only :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, Qnl, angularint
!expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y l del ECP
	implicit none
	integer :: w,n,z,l
	integer, intent(in) :: i,j,k,ii,ji,lx,ly,lz
	double precision :: Ccoef
	z=ZlistECP(k)
	L=Lmax(ZlistECP(k))
	n=lx+ly+lz
	AAAlocal=0.d0
	do w =1, expnumbersECP(z,l)
!barre todos los terminos del Lmaximo
		Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
		AAAlocal=AAAlocal+aECP(z,L,w)*angularint(lx,ly,lz)*Q0(n+nECP(z,l,w),Ccoef)
	end do
	return
	end function AAAlocal

	DOUBLE PRECISION function Q0(n,alpha)
!calcula  int r^n * exp(-alpha * r^2) dr entre 0 e infinito
	use ECP_mod, only :doublefac, pi12, fac
	implicit none
	integer, intent(in) :: n
	double precision, intent(in) :: alpha
	if ( .not. mod(n,2)) then
                Q0=0.5d0*pi12/sqrt(alpha) * doublefac(n-1)/((2*alpha)**(n/2))
		return
	else
		Q0=fac((n-1)/2)/(2*alpha**((n+1)/2))
		return
	end if
	end function Q0



!subrutinas angulares

	DOUBLE PRECISION FUNCTION OMEGA2(K,lambda,l,m,a,b,c)
!Esta funcion devuelve el valor de omega2 evaluado en el vector K
!OMEGA2(K,lambda,l,m,a,b,c) = Suma(o=-lambda a o=lambda) ( Y lambda o(K)*  int (x/r)^a * (y/r)^b * (z/r)^c Ylambda o() Ylm()  d angulo solido )

!no esta testeada a full ya q es muy engorroso
	use ECP_mod, only : intangular,angularint
        implicit none
        DOUBLE PRECISION, intent(in), dimension(3) :: K
        integer, intent(in) :: lambda,l,m,a,b,c
        DOUBLE PRECISION,dimension(3) :: Kun
        integer :: o,r,s,t,u,v,w
        Double precision :: SUM1, SUM2
	Kun=K/sqrt(K(1)**2 + K(2)**2 + K(3)**2)
	SUM1=0.d0
        SUM2=0.d0
	OMEGA2=0.d0
!	write(15,*) "lkkk"
!	write(15,*) "lambda,l,m,a,b,c",lambda,l,m,a,b,c
	do o=-lambda,lambda
		do r=0,lambda
			do s=0,lambda-r
				t=lambda-r-s
                                SUM1=SUM1+Ucoef(lambda,o,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
!				write(*,*) "OMEGA2"
!				write(15,*) "sum1",sum1
				do u=0,l
					do v = 0,l-u
						w=l-u-v
			                        SUM2=SUM2+Ucoef(lambda,o,r,s,t)*Ucoef(l,m,u,v,w)*angularint(a+r+u,b+s+v,c+t+w)
!						write(*,*) "OMEGA2.2"
!						write(15,*) "a+r+u,b+s+v,c+t+w",a+r+u,b+s+v,c+t+w
!						write(15,*) "Ucoef(lambda,o,r,s,t)*Ucoef(l,m,u,v,w)*angularint(a+r+u,b+s+v,c+t+w)",Ucoef(lambda,o,r,s,t),Ucoef(l,m,u,v,w),angularint(a+r+u,b+s+v,c+t+w)
!						write(15,*) "sum2",sum2,"u,v,w",u,v,w
					end do
				end do
			end do
		end do
!		write(15,*) "sum1",sum1,"sum2",sum2
		OMEGA2=OMEGA2+SUM1*SUM2
!		write(15,*) o,omega2
                SUM1=0.d0
                SUM2=0.d0
	end do
	end function OMEGA2

	DOUBLE PRECISION FUNCTION OMEGA1(K,l,a,b,c)
!Esta funcion devuelve el valor de omega evaluado en el vector K
!OMEGA1(K,l,a,b,c) = Suma(u=-l a u=l) ( Ylu(K)* int(x/r)^a * (y/r)^b * (z/r)^c Ylu() d angulo solido )
	use ECP_mod, only : intangular,angularint
	implicit none
	DOUBLE PRECISION, intent(in), dimension(3) :: K
	integer, intent(in) :: l,a,b,c 
	DOUBLE PRECISION,dimension(3) :: Kun
	integer :: r,s,t,u
	Double precision :: SUM1, SUM2
	SUM1=0.d0
	SUM2=0.d0
	OMEGA1=0.d0
	if ( all(K .eq. (/0.d0,0.d0,0.d0/))) then
!caso especial para los terminos <A|A|A>
	return
	end if


	Kun=K/sqrt(K(1)**2 + K(2)**2 + K(3)**2)

!	write(*,*) "l",l
	do u=-l,l
!		write(*,*) "u",u
		do r=0,l
!			write(*,*) "r",r
			do s=0,l-r
!				write(*,*) "s",s
				t=l-r-s
!				write(*,*) "t",t
				SUM1=SUM1+Ucoef(l,u,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
!				write(*,*) "OMEGA1"
				SUM2=SUM2+Ucoef(l,u,r,s,t)*angularint(a+r,b+s,c+t)
!Usa la funcion angular int en lugar del array int angular. Hay que agrandar el array
!				write(*,*) SUM1, SUM2, angularint(a+r,b+s,c+t),a+r,b+s,c+t
			end do
		end do
!		write(*,*) SUM1, SUM2, intangular(a+r,b+s,c+t),a+r,b+s,c+t
		OMEGA1=OMEGA1+SUM1*SUM2
		SUM1=0.d0
		SUM2=0.d0
	end do
	return
	end function OMEGA1


	DOUBLE PRECISION FUNCTION Ucoef(l,m,lx,ly,lz)
!lz esta al pedo, puedo usarlo para chekear que lx+ly+lz=l
	use ECP_mod, only :  l0,l1,l2,l3,l4
	implicit none
	integer, intent(in) :: l,m,lx,ly,lz

	if (lx+ly+lz .ne. l) then
		write(*,*) "problem with Ucoef: lx+ly+lz not equal l"
	end if

	if (l .eq. 0) then
		Ucoef=l0(1)
	elseif (l .eq. 1) then
		Ucoef=l1(2*lx+ly+1,m)
	elseif (l .eq. 2) then
		Ucoef=l2(0.5*lx*(7-lx)+ly+1,m)
        elseif (l .eq. 3) then
		Ucoef=l3(0.5*lx*(9-lx)+ly+1,m)
	elseif (l .eq. 4) then
		Ucoef=l4(0.5*lx*(11-lx)+ly+1,m)
	else
	write(*,*) "ECP error l is greater than 4",l
	end if
	return
	end function Ucoef



        end subroutine intECP


!subrutinas Radiales

        subroutine Qtype1(K,Ccoef,lmax,necp)
!agrega a la matriz Qnl los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl ya que hay q llamar a esta rutina por cada termino del pseudopotencial
!luego hay que optimizar la rutina para que calcule Qnl a partir de otros Qnl
	use ECP_mod, only :  alpha, betha, Bn, Cn, Qnl
	implicit none
	double precision, intent(in) :: K,Ccoef
!lmax  = 0 para <s||s>, 1 para <s||p>, 2 para <s||d>, ... , 6 para <d||d>	
!necp corresponde al exponente del paseudopotencial en r^n
	integer, intent(in) :: necp,lmax
!nmin y nmax dan el rango en que tiene q calcular Q
	integer :: nmin,nmax
!variables auxiliares
	integer :: n,l,i
        double precision :: acoef, gam, acum
!	write(*,*) "Q recibe:",  K,Ccoef,lmax,necp
!Caso 1-calcula todas las integrales
	nmin=necp
	nmax=necp+lmax
!	write(*,*) "Q, k", K
	acoef=K/(2.d0*Ccoef)
	gam=0.5d0*exp(K**2/(4*Ccoef))
	Bn=0.d0
	Cn=0.d0
	call ByC(acoef,ccoef,nmin,nmax,Bn,Cn)


	do n=nmin,nmax
	do l=0,lmax

	acum=0.d0
	do i=l,1,-2
		acum=acum+alpha(l,i)*Bn(n-i)/k**i
!		write(*,*) "acum",acum
	end do
	do i=l+1,1,-2
		acum=acum+betha(l+1,i)*Cn(n-i)/k**i
!		write(*,*) "acum",acum
	end do
		Qnl(n,l)=Qnl(n,l)+acum*gam
!		write(*,*) "Q",Qnl(n,l),acum,gam
		acum=0.d0
	end do
	end do
!	write(*,*) "Q00",Qnl(0,0),"Q10", Qnl(1,0),"Q20", Qnl(2,0)
!escribo todo como test
!	write(*,*) "entro en QNL"
!	do n=nmin,nmax
!        do l=0,lmax
!	write(*,*) n,l,Qnl(n,l)
!	end do
!	end do
        End subroutine Qtype1


        subroutine ByC(acoef,ccoef,nmin,nmax,Barray,Carray)
!calcula los coeficientes B y C 
	use ECP_mod, only : DAW,DAWERF,NEXTCOEF,pi,pi12
	IMPLICIT NONE
!acoef,ccoef son los coeficientes para el calculo de B y C
! A(ccoef,acoef)= int exp(-ccoef*x^2)(x+acoef)^n dx from -acoef to inf
!nmin,nmax delimitan el rango de n
!Bn=An(ccoef,acoef) + An(ccoef,-acoef)
!Cn=An(ccoef,acoef) - An(ccoef,-acoef)
	DOUBLE PRECISION, intent(in) :: acoef,ccoef
	integer, intent(in) :: nmin,nmax
	DOUBLE PRECISION, intent(inout), dimension (-12:14) :: Barray,Carray
!variables auxiliares
	DOUBLE PRECISION :: c0sq,ncos,ca
	integer :: i
!	write(*,*) "BYC", acoef,ccoef,nmin,nmax!,Barray,Carray
	c0sq=sqrt(ccoef)
	Barray(0)=pi12/c0sq
	Carray(0)=Barray(0)*erf(acoef*c0sq)

	if (nmax>0) then
	        Barray(1)= exp(-ccoef*acoef**2)/ccoef + acoef*Carray(0)
	        Carray(1)= Barray(0)*acoef
		do i=2,nmax
			ncos=(i-1)/(2*ccoef)
			Barray(i)=ncos*Barray(i-2)+acoef*Carray(i-1)
			Carray(i)=ncos*Carray(i-2)+acoef*Barray(i-1)
		end do
	end if

        if (nmin<0) then
		Barray(-1)=2*pi12*DAWERF(acoef*c0sq)
		Carray(-1)=2*pi12*DAW(acoef*c0sq)
		if (nmin<-1) then
			ca=2*ccoef*acoef
			Barray(-2)=ca*Carray(-1)-2*ccoef*Barray(0)
			Carray(-2)=2*ca*exp(-ccoef*acoef**2)+ca*Barray(-1)-2*ccoef*Carray(0)
			do i=-3,nmin,-1
				Barray(i)=NEXTCOEF(1,i,ca,exp(-ccoef*acoef**2),ccoef,Carray(i+1), Barray(i+2))
				Carray(i)=NEXTCOEF(-1,i,ca,exp(-ccoef*acoef**2),ccoef,Barray(i+1), Carray(i+2))
			end do
		end if
	end if
        end subroutine ByC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!this routine obtain Q(l1,l2,n,k1,k2,a) where
!Q(l1,l2,n,k1,k2,a)= int r^n exp(-ar^2) Ml1(k1r) Ml2(k2r) dr form r=0 to inf
!where Ml are the modified spherical Bessel function of the first kind.

        subroutine Qtype2(Ka,Kb,Ccoef,l1max,l2max,necp)
!agrega a la matriz Qnl1l2 los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl1l2 ya que hay q llamar a esta rutina por cada termino del pseudopotencial
!luego hay que optimizar la rutina para que calcule Qnl1l2 a partir de otros Qnl1l2
	use ECP_mod, only :  alpha, betha, rho, tau, sigma, sigmaR, Qnl1l2
	implicit none
	double precision, intent(in) :: Ka,Kb,Ccoef
!l1max y l2max = 0 para s, 1 para p, 2 para d, etc	
!n corresponde al exponente del paseudopotencial en r^n
	integer, intent(in) :: necp,l1max, l2max
	integer :: nmin,nmax
!variables auxiliares
	integer :: i,j,n,l1,l2
	double precision :: alfok, betok, acum1
	
!Caso 1-calcula todas las integrales
	nmin=necp
	nmax=necp+l1max+l2max
!	write(*,*) nmin,nmax
	call integrals(Ka,Kb,Ccoef,necp-2-l1max-l2max,necp+l1max+l2max-2)
	acum1=0.d0

	do n=nmin,nmax
	do l1=0,l1max
	do l2=0,l2max
	do i=l1,1,-2
		alfok=alpha(l1,i)/Ka**i
		do j=l2,1,-2
			if (tau(n-i-j) == 0.d0) write(*,*) "Error, no calculo tau",n-i-j
			acum1=acum1+alpha(l2,j)*tau(n-i-j)/Kb**j
		end do
		do j=l2+1,1,-2
			acum1=acum1+betha(l2+1,j)*sigmaR(n-i-j)/Kb**j
			if (sigmaR(n-i-j) == 0.d0) write(*,*) "Error, no calculo sigmaR",n-i-j
		end do
		Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*alfok
		acum1=0.d0
	end do
	do i=l1+1,1,-2
		betok=betha(l1+1,i)/Ka**i
		do j=l2,1,-2
			acum1=acum1+alpha(l2,j)*sigma(n-i-j)/Kb**j
			if (sigma(n-i-j) == 0.d0) write(*,*) "Error, no calculo sigma",n-i-j
		end do
		do j=l2+1,1,-2
			acum1=acum1+betha(l2+1,j)*rho(n-i-j)/Kb**j
			if (rho(n-i-j) == 0.d0) write(*,*) "Error, no calculo rho",n-i-j
		end do
		Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*betok
                acum1=0.d0
	end do

	end do
	end do
	end do

!escribo todo como chekeo
	if (1) then
	write(*,*) "   n    l1    l2    Qnl1l2"
        do n=nmin,nmax
        do l1=0,l1max
        do l2=0,l2max
	write(*,*) n,l1,l2,Qnl1l2(n,l1,l2)
	end do
        end do
        end do
	end if

        End subroutine Qtype2

	subroutine integrals(Ka,Kb,Ccoef,nmin,nmax)
!obtain integrals rho, tau, sigma and sigmaR from n between nmin and nmax
! rho(n) = int exp(-cr^2) * sinh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! sigma(n) = int exp(-cr^2) * sinh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf
! sigmaR(n) = int exp(-cr^2) * cosh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! tau(n) = int exp(-cr^2) * cosh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf

	use ECP_mod, only : Bn1,Bn2,Cn1,Cn2,rho, tau, sigma, sigmaR
	implicit none
	integer, intent(in) :: nmin,nmax
	double precision, intent(in) :: Ka,Kb,Ccoef
	double precision, dimension(2) :: acoef,gammacoef
	double precision :: signo
	acoef(1)=0.5d0*(Ka+Kb)/Ccoef
	acoef(2)=0.5d0*abs(Ka-Kb)/Ccoef
	gammacoef(1)=0.25d0*exp(Ccoef*acoef(1)**2)
	gammacoef(2)=0.25d0*exp(Ccoef*acoef(2)**2)
!	call ByCdoble(acoef(1),acoef(2),Ccoef,nmin,nmax)
	Bn1=0.d0
	Bn2=0.d0
	Cn1=0.d0
	Cn2=0.d0
	rho=0.d0
	tau=0.d0
	sigma=0.d0
	sigmaR=0.d0
	call ByC(acoef(1),Ccoef,nmin,nmax,Bn1,Cn1)
	call ByC(acoef(2),Ccoef,nmin,nmax,Bn2,Cn2)
	rho=gammacoef(1)*Bn1-gammacoef(2)*Bn2
	tau=gammacoef(1)*Bn1+gammacoef(2)*Bn2
	sigma=gammacoef(1)*Cn1+sign(1d0,Ka-Kb)*gammacoef(2)*Cn2
	sigmaR=gammacoef(1)*Cn1-sign(1d0,Ka-Kb)*gammacoef(2)*Cn2

	end subroutine integrals

        subroutine ByCdoble(a1coef,a2coef,c0coef,nmin,nmax)
!esta subrutina ya no la usa!!!!!
!calcula los coeficientes B y C para las integrales rho tau sigma y sigma raya
	use ECP_mod, only : Bn1,Bn2,Cn1,Cn2,DAW,DAWERF,NEXTCOEF,pi,pi12
	IMPLICIT NONE
!a1coef,a2coef,c0coef son los coeficientes para el calculo de B y C
! Ai(c0coef,aicoef)= int exp(-c0coef*x^2)(x+aicoef)^n dx from -aicoef to inf
!nmin,nmax delimitan el rango de n
!Bni=An(c0coef,aicoef) + An(c0coef,-aicoef)
!Cni=An(c0coef,aicoef) - An(c0coef,-aicoef)
	DOUBLE PRECISION, intent(in) :: a1coef,a2coef,c0coef
	integer, intent(in) :: nmin,nmax
!variables auxiliares
	DOUBLE PRECISION :: c0sq,ncos,ca1,ca2
	integer :: i
!	write(*,*) "byc",nmin,nmax
	c0sq=sqrt(c0coef)
	Bn1(0)=pi12/c0sq
	Bn2(0)=Bn1(0)
	Cn1(0)=Bn1(0)*erf(a1coef*c0sq)
	Cn2(0)=Bn1(0)*erf(a2coef*c0sq)

	if (nmax>0) then
	        Bn1(1)= exp(-c0coef*a1coef**2)/c0coef + a1coef*Cn1(0)
	        Bn2(1)= exp(-c0coef*a2coef**2)/c0coef + a2coef*Cn2(0)
	        Cn1(1)= Bn1(0)*a1coef
	        Cn2(1)= Bn1(0)*a2coef
		do i=2,nmax
			ncos=(i-1)/(2*c0coef)
			Bn1(i)=ncos*Bn1(i-2)+a1coef*Cn1(i-1)
			Bn2(i)=ncos*Bn2(i-2)+a2coef*Cn2(i-1)
			Cn1(i)=ncos*Cn1(i-2)+a1coef*Bn1(i-1)
			Cn2(i)=ncos*Cn2(i-2)+a2coef*Bn2(i-1)
		end do
	end if

        if (nmin<0) then
		Bn1(-1)=2*pi12*DAWERF(a1coef*c0sq)
		Bn2(-1)=2*pi12*DAWERF(a2coef*c0sq)
!		write(*,*) "DAWERF(a2coef*c0sq)",DAWERF(a2coef*c0sq)
		Cn1(-1)=2*pi12*DAW(a1coef*c0sq)
		Cn2(-1)=2*pi12*DAW(a2coef*c0sq)
		if (nmin<-1) then
			ca1=2*c0coef*a1coef
			ca2=2*c0coef*a2coef
			Bn1(-2)=ca1*Cn1(-1)-2*c0coef*Bn1(0)
			Bn2(-2)=ca2*Cn2(-1)-2*c0coef*Bn2(0)
			Cn1(-2)=2*ca1*exp(-c0coef*a1coef**2)+ca1*Bn1(-1)-2*c0coef*Cn1(0)
			Cn2(-2)=2*ca2*exp(-c0coef*a2coef**2)+ca2*Bn2(-1)-2*c0coef*Cn2(0)
			do i=-3,nmin,-1
				Bn1(i)=NEXTCOEF(1,i,ca1,exp(-c0coef*a1coef**2),c0coef,Cn1(i+1), Bn1(i+2))
				Bn2(i)=NEXTCOEF(1,i,ca2,exp(-c0coef*a2coef**2),c0coef,Cn2(i+1), Bn2(i+2))
				Cn1(i)=NEXTCOEF(-1,i,ca1,exp(-c0coef*a1coef**2),c0coef,Bn1(i+1), Cn1(i+2))
				Cn2(i)=NEXTCOEF(-1,i,ca2,exp(-c0coef*a2coef**2),c0coef,Bn2(i+1), Cn2(i+2))
			end do
		end if
	end if

        end subroutine ByCdoble

	subroutine ReasignZ
!cambia la carga del nucleo sacandole la carga del core y guarda las cargas originales en IzECP
!tambien corrige la cantidad de electrones sacando los que van al core
	use garcha_mod, only :Iz,nuc,nshell, natom, NCO
	use ECP_mod, only : ZlistECP,IzECP,Zcore,ecptypes
	implicit none
	integer :: i,j
	allocate (IzECP(natom))
	do i=1, natom
!barre atomos
		IzECP(i)=Iz(i)
		do j=1, ecptypes
!barre atomos con ecp
			if (IzECP(i) .eq. ZlistECP(j)) then
				write(*,*) "editando atomo", i
				Iz(i)=Iz(i)-Zcore(ZlistECP(j))
				write(*,*) "carga nuclear efectiva = ", Iz(i)
!cambia la carga del nucleo
				NCO=NCO-Zcore(ZlistECP(j))/2
!saca e-
				write(*,*) "electrones del sistema removidos ", Zcore(ZlistECP(j))
!				write(*,*) "carga nueva" , Iz(i)
!				write(*,*) "caraga guardada", IzECP(i)
			end if
		end do
	end do

	end subroutine ReasignZ

!	9018 format(/1x,'Z =',i4,2x, 'L =',i4,2x,'coefnumber =',i3,2x, 'n =',i2,2x, 'b =', f15.5,2x,'c =',f15.5)
!        end subroutine intECP


        subroutine obtainls
!arma una matriz que contenga los exponentes de la parte ngular de la base
! |x> = A x^lx y^ly z^lz *e^-ar^2
! devuelve angularL(i,j) j=1 lx, j=2, ly, j=3 lz para la funcion i de la base
        use garcha_mod, only : nshell
        use ECP_mod, only : Lxyz
        implicit none
        integer :: i
        integer :: D, lx,ly,lz
        integer :: resto
        D = nshell(0)+nshell(1)+nshell(2)
        allocate (Lxyz(D, 3))
        Lxyz=0
        do i=nshell(0)+1,D
                if (i .le. nshell(0)+nshell(1)) then
!funciones p
                        resto=modulo(i-nshell(0),3)
                        if (resto .eq. 1) Lxyz(i,1)=1
                        if (resto .eq. 2) Lxyz(i,2)=1
                        if (resto .eq. 0) Lxyz(i,3)=1
                else if (i .le. D) then
                        resto=modulo(i-nshell(0)+nshell(1),6)
                        if (resto .eq. 1) Lxyz(i,1)=2
                        if (resto .eq. 2) then
                                Lxyz(i,1)=1
                                Lxyz(i,2)=1
                        end if
                        if (resto .eq. 3) Lxyz(i,2)=2
                        if (resto .eq. 4) then
                                Lxyz(i,1)=1
                                Lxyz(i,3)=1
                        end if
                        if (resto .eq. 5) then
                                Lxyz(i,2)=1
                                Lxyz(i,3)=1
                        end if
                        if (resto .eq. 0) Lxyz(i,3)=2
!		else
!			write(*,*) "sobra este i",i
                end if
                
        end do

!testeo escribiendo
!        do i=1,D
!                write(*,*) i,Lxyz(i,1),Lxyz(i,2),Lxyz(i,3)
!        end do

!        write(*,*) nshell(0),nshell(1),nshell(2)


        end subroutine obtainls

