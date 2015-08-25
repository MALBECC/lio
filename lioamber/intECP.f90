
	subroutine intECP(tipodecalculo)
	use garcha_mod, only :nshell,nuc,a,c, ncont, natom
	use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP,nECP,bECP, aECP,Zcore, Lmax, expnumbersECP,VAAAcuadrada,lxyz, VAAA, VAAAcuadrada, VAAB, VAABcuadrada,pi,doublefac,VXXXgamess,distx,disty,distz,VBAC,VBACcuadrada,verbose_ECP,ecp_debug
	implicit none

	integer, intent(in) :: tipodecalculo
!tipodecalculo=1 alocatea variables y alcula terminos AAA
!tipodecalculo=2 calcula terminos ABB y ABC
	integer z,l,t,ji,ii,jk,ik,jii,iii
	integer m1,lx,ly,lz
        integer :: ns,np,nd,M,i,j,e,pos,n
	double precision :: alpha
        ns=nshell(0)
        np=nshell(1)
        nd=nshell(2)
        M=ns+np+nd
	if (tipodecalculo .eq. 1) then
        call allocate_ECP()
!allocatea la matriz de fock de pseudopotenciales, la matrix cuadrara para pruebas y el verctor con los terminos de 1 electron sin corregir
!edita normalizacion base d, luego lo cambiare aguadabdo los coeficientes en un array
	call correctbasis(1)
!prepara variables y calcula los terminos AAA	
        call ReasignZ()
!reasigna las cargas
        call obtainls()
!obtiene una matriz con lx,ly y lz
	call intECPAAA()
	call correctbasis(-1)
	write(*,*) "termino ECP tipo 1"
	go to 325
!calcula terminos AAA
	elseif (tipodecalculo .eq. 2) then
        call correctbasis(1)
	call obtaindistance()
!ontiene arrays con la diatncia en x, y y z entre cada par de atomos i y j
	call intECPAAB()
!	write(90,*) VAAB
!calculo VAAB
        call correctbasis(-1)
	write(*,*) "termino ECP tipo 2"
	go to 325
	elseif (tipodecalculo .eq. 3) then
!calculo VABC
        call correctbasis(1)
	call intECPABC()
        call correctbasis(-1)
	write(*,*) "termino ECP tipo 3"
		go to 324
	elseif (tipodecalculo .eq. 4) then
		call deallocateV()
		go to 325
        else
                Write(*,*) "ERROR in tipe of ECP calculation"
	endif




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




        if (tipodecalculo .eq. 2) then
        if ( .false. ) then
!escribe coeficientes de fock del ECP AAB y AAA
        write(*,*) "VAAA+VAAB"
!        do i=1,M
!                do j=i,M
!barrido de gamess
	write(*,*) nshell(0),nshell(1),nshell(2),M, "chequeardimm, Nick"
	VXXXgamess=0.d0
!copia VAAAcuadrada a VXXXgamess cambiando el orden para que sea el mismo que en gamess. solo para testeto
	do j=1,M
		do i=1,j
			if (j .le. nshell(0)+nshell(1) .and. i .le. nshell(0)+nshell(1)) then
				VXXXgamess(i,j)=VAAAcuadrada(i,j)+VAABcuadrada(i,j)
			else
				ji=j
				ii=i
				jk=0
				ik=0
				if (j .gt. nshell(0)+nshell(1)) then
					
					ji=j-nshell(0)-nshell(1)
					jk=ji/6
!					write(*,*) "*********************",j,ji,jk
					ji=modulo(ji,6)
!					write(*,*) "****************", ji
					if (ji .eq. 0) jii=-3
					if (ji .eq. 1) jii=1
                                        if (ji .eq. 2) jii=4
                                        if (ji .eq. 3) jii=2
                                        if (ji .eq. 4) jii=5
                                        if (ji .eq. 5) jii=6
					if (ji .gt. 5) stop "error ji gt 5"
					ji=jii+jk*6+nshell(0)+nshell(1)
!					write(*,*) "**** nuevo ji", ji
!					write(*,*) "***",j,ji,jii,jk
				end if
				if (i .gt. nshell(0)+nshell(1)) then
                                        ii=i-nshell(0)-nshell(1)
					ik=ii/6
                                        ii=modulo(ii,6)
                                        if (ii .eq. 0) iii=-3
                                        if (ii .eq. 1) iii=1
                                        if (ii .eq. 2) iii=4
                                        if (ii .eq. 3) iii=2
                                        if (ii .eq. 4) iii=5
                                        if (ii .eq. 5) iii=6
					ii=iii+ik*6+nshell(0)+nshell(1)
				end if
				VXXXgamess(ii,ji)=VAAAcuadrada(i,j)+VAABcuadrada(i,j)
!				write(*,*) "viejo",i,j,"nuevo",ii,ji
			end if
		end do
	end do
	

	n=1
	do j=1,M
		do i=1,j
                        write(*,*) VXXXgamess(i,j),",",i,",",j,",",n
			n=n+1
                end do
        end do
        end if
	end if






 324    write(*,*) "ahora rutinas de debug y verbose"

	if ( verbose_ECP .gt. 0) then
		call WRITE_BASIS()
		call WRITE_ECP_PARAMETERS()
		call WRITE_FOCK_ECP_TERMS()
		call WRITE_FOCK_ECP()
		call WRITE_ANG_EXP()
		call WRITE_DISTANCE()
	end if

	if ( ecp_debug ) then
		call SEARCH_NAN()
		call TEST_SYMMETRY()
		call TEST_COPY_FOCK()
	end if

 325    write(*,*) "terminaron calculos de ECP version betha"

	9012 format(/1x,i2,2x,i2,2x,i2,2x,i2,2x,i2,2x,f18.10)
	9013 format(/1x,"vector",f18.10,2x,"matriz",f18.10,2x,"diff",f18.10)
	9014 format(/1x,"i",i3,2x,"lx",i2,2x,"ly",i2,2x,"lz",i2)
	9015 format(/1x,"i",i3,2x,"j",i3,2x,"Vij",f10.5,2x,"Vji",f10.5,2x,     "diff",f18.15)
        9018 format(/1x,'Z =',i4,2x, 'L =',i4,2x,'coefnumber =',i3,2x,  'n =',i2,2x, 'b =', f15.5,2x,'c =',f15.5)
	9019 format(/1x,i3,2x,i3,2x,f18.10)
	1111 format(/1x,"i",i4,1x,"j",i4,1x,"l",i2,1x,"m",i2,1x,"a2",f18.10,1x,"q",f18.10)

	contains



!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!terminos de un centro <Xa|Va|Xa>


        Subroutine intECPAAA
!calcula los terminos de Fock para bases y pseudopotenciales en el mismo atomo
        use garcha_mod, only : a,c,ncont, nshell, nuc, ncont
!a(i,ni) exponente de la funcion de base i, contrccion ni
!c(i,ni) coeficiente de la funcion de base i, contrccion ni
!ncont(i) cantidad de contracciones de la funcion de base i
!nshell(i) cantidad de funciones i=1 s, i=2, p, i=3, d
        use ECP_mod, only :nECP,bECP, aECP, ecptypes, IzECP, Lmax, Lxyz, VAAAcuadrada,VAAA,local_nonlocal, ecp_debug
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
!       lmaxQ=0
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
                          AAA=AAA_LOCAL(i,j,k,ii,ji,lxi+lxj,lyi+lyj,lzi+lzj) + AAA_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj)

                          if (local_nonlocal .eq. 1 .and. ecp_debug) then
! local_nonlocal 1 y 2 son solo para debugueo, este pedazo hay q sacarlo al final
                             AAA=AAA_LOCAL(i,j,k,ii,ji,lxi+lxj,lyi+lyj,lzi+lzj)
                          else if (local_nonlocal .eq. 2 .and. ecp_debug) then
                             AAA= AAA_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj)
                          end if

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
                          VAAA(pos) = VAAA(pos) + acum*c(i,ii) !esta linea es lo unico que quedaria!!!!!!!!!!!!!
                             if (abs(VAAA(pos)-VAAAcuadrada(i,j)) .gt. 0.000000000001) then
                                do e=1,20
                                write(*,*) i,j,pos,"arma mal el vector!"
                                end do
                             end if
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





        DOUBLE PRECISION function AAA_LOCAL(i,j,k,ii,ji,lx,ly,lz)
!Calcula el termino local del pseudopotencial
!   [<xi|V(LM)|xj>] = 
!=Ni Nj C(LM) int(r^n exp(-alpha *r^2) dr, 0, inf) * int ((x/r)^ni (y/r)^li (z/r)^mi d(angulo solido)
!los coef de la base se multiplican en la rutina que llama a esta

        use garcha_mod, only : a,c,nshell
        use ECP_mod, only :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, Qnl, angularint
!expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y l del ECP
        implicit none
        integer :: w,n,Z,l
        integer, intent(in) :: i,j,k,ii,ji,lx,ly,lz
!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
!Z carga del nucleo
        double precision :: Ccoef
        Z=ZlistECP(k)
        L=Lmax(Z)
        n=lx+ly+lz
        AAA_LOCAL=0.d0
        do w =1, expnumbersECP(z,l)
!barre todos los terminos del Lmaximo
           Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
           AAA_LOCAL=AAA_LOCAL+aECP(z,L,w)*angularint(lx,ly,lz)*Q0(n+nECP(z,l,w),Ccoef)
        end do
        return
        end function AAA_LOCAL

        DOUBLE PRECISION function AAA_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj)
!calcula el termino semilocal del pseudopotencial
!suma m=-l hasta l  [<xi|lm> V(l-LM) <lm|xj>] = 
!=Ni Nj C(l-LM) suma m=-l, l int(r^n exp(-alpha *r^2) dr, 0, inf) * int ((x/r)^ni (y/r)^li (z/r)^mi Ylm d(angulo solido) *
! * int ((x/r)^nj (y/r)^lj (z/r)^mj Ylm d(angulo solido)
!los coef de la base se multiplican en la rutina que llama a esta
        use garcha_mod, only : a,c
        use ECP_mod, only :ZlistECP,Lmax,aECP,nECP,bECP, expnumbersECP
        implicit none
        integer, intent(in) :: i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
!lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
        integer :: l,m,term
!auxiliades para ciclos
        integer :: Z,n
!Z= carga del nucleo
        double precision :: A2, Acoef
        integer :: suma1, suma2
        AAA_SEMILOCAL=0.d0
        A2=0.d0
        Z=ZlistECP(k)
        n=lxi+lxj+lyi+lyj+lzi+lzj

        do l = 0 , Lmax(z)-1
!barre todos los l de la parte no local
!           do m=-l,l     
!barre m
!              A2=A2+Aintegral(l,m,lxi,lyi,lzi)*Aintegral(l,m,lxj,lyj,lzj)
!A2 contiene la parte angular de la integral
!           end do

           do term=1, expnumbersECP(z,l)
!barre contracciones del ECP para el atomo con carga z y l del ecp

!%%%%%%%%%%%%%%%%%%%%asi estaba antes!!!!!!!!!!!!!
              do m=-l,l     !este do se tiene que poder poner antes del do de term=1, expnumbersECP(z,l)
!barre m
                 A2=A2+Aintegral(l,m,lxi,lyi,lzi)*Aintegral(l,m,lxj,lyj,lzj)
!A2 contiene la parte angular de la integral
              end do
!%%%%%%%%%%%%%%%%%%%%asi estaba antes!!!!!!!!!!!!!

              Acoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
!Acoef es el exponente de la integral radial
              AAA_SEMILOCAL=AAA_SEMILOCAL+A2*aECP(z,L,term)*Q0(n+nECP(z,l,term),Acoef)
!Q0 integral radial   Q0(n,alpha)= int r^n exp(-alpha * r^2) dr fron 0 to inf
!aECP coeficiente del ECP

!%%%%%%%%%%%%%%%%%%%%asi estaba antes!!!!!!!!!!!!!
              A2=0.d0
!%%%%%%%%%%%%%%%%%%%%asi estaba antes!!!!!!!!!!!!!

           end do
!	   A2=0.d0
        end do
        return
        end function AAA_SEMILOCAL



!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!terminos de 2centros <Xa|Va|Xb>,<Xb|Va|Xa> NO incluye el caso <Xb|Va|Xb>


	subroutine intECPAAB()
        use garcha_mod, only : nshell,nuc,a,c,ncont
!a(i,ni) exponente de la funcion de base i, contrccion ni
!c(i,ni) coeficiente de la funcion de base i, contrccion ni
!ncont(i) cantidad de contracciones de la funcion de base i
!nshell(i) cantidad de funciones i=1 s, i=2, p, i=3, d
        use ECP_mod, only : ecptypes,cutECP,IzECP,cutecp2,distx, disty, distz,Lxyz, VAAB, VAABcuadrada,local_nonlocal, ecp_debug
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
        do i=1, M
!barre funciones de la base 
        do j=1, M
!cambiar por do j=i,M para barrer solo la mitad de la matriz
!barre el otro coef de la base j>=i ya que la matriz tiene q ser simetrica
           if (nuc(i) .ne. nuc(j)) then
!agarra bases de atomos distintos
              do k=1, ecptypes
!barre atomos con ecp
	         if (IzECP(nuc(i)) .eq. ZlistECP(k) .or. IzECP(nuc(j)) .eq. ZlistECP(k)) then
!solo calcula si el atomo tiene ECP
	            dx=distx(nuc(i),nuc(j))
	            dy=disty(nuc(i),nuc(j))
	            dz=distz(nuc(i),nuc(j))
	            Distcoef=(dx**2.d0 + dy**2.d0 + dz**2.d0)
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
			if (Distcoef*a(j,ji) .ge. cutecp2) write(*,*) "cut 2 sobro omit int"
	                  if ( .not. cutECP .or. (Distcoef*a(j,ji) .lt. cutecp2)) then
!solo calcula los terminos que luego se multipliquen por un factor q no sea demasiado pequeño,
!el cut va a ser obligatorio, ya que si los atomos quedan muy separados el termino VAAB se obtiene como el producto de un numero MUY grande por uno que es casi 0
!y se obtienen NAN para el numero muy grande
!a distancias grades gana el termino que es casi 0
	                     acum=0.d0
	                     do ii=1, ncont(i)
!								acum=0.d0
!ji y ii barren contracciones de las funcion de base
!me conviene que barra primero los terminos de la base de i, ya que el factor ezponencial que multiplica solo depende de j

	                        AAB=AAB_LOCAL(i,j,k,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz) +AAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)

                                if (local_nonlocal .eq. 1 .and. ecp_debug) then
! local_nonlocal 1 y 2 son solo para debugueo
	                           AAB=AAB_LOCAL(i,j,k,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
                                else if (local_nonlocal .eq. 2 .and. ecp_debug) then
	                           AAB=AAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
                                end if
                             acum=acum+AAB*c(i,ii)
	                     AAB=0.d0
!multiplica por el coeficiente de la base
	                     end do
                             VAABcuadrada(i,j) = VAABcuadrada(i,j) + acum*c(j,ji)*4*pi*exp(-Distcoef*a(j,ji))
!multiplica por el otro coef de la base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             if (i .ge. j) then
!este if hay q sacarlo al fina, cuando cambie el rango en el q barre j , en vez de comenzar en 1 comience en i
                                pos=i+(1-j)*(j-2*M)/2   
                                VAAB(pos) = VAAB(pos) + + acum*c(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))!esta linea es lo unico que quedaria!!!
                             end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	                  end if
                       end do
	            end if
                    if (IzECP(nuc(j)) .eq. ZlistECP(k)) then
!calculo para ECP en j
	               do ii=1, ncont(i)
			 if (Distcoef*a(i,ii) .ge. cutecp2) write(*,*) "cut 2 sobro omit int"
                          if ( .not. cutECP .or. (Distcoef*a(i,ii) .lt. cutecp2)) then
!solo calcula los terminos que luego se multipliquen por un factor q no sea demasiado pequeño,
                             acum=0.d0
                             do ji=1, ncont(j)
!ji y ii barren contracciones de las funcion de base
!me conviene que barra primero los terminos de la base de j, ya que el factor ezponencial que multiplica solo depende de j
                                AAB=AAB_LOCAL(j,i,k,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz) +AAB_SEMILOCAL(j,i,ji,ii,k,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz)

                                if (local_nonlocal .eq. 1 .and. ecp_debug) then
! local_nonlocal 1 y 2 son solo para debugueo
                                   AAB=AAB_LOCAL(j,i,k,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)
                                else if (local_nonlocal .eq. 2 .and. ecp_debug) then
                                   AAB=AAB_SEMILOCAL(j,i,ji,ii,k,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz)
                                end if

                                acum=acum+AAB*c(j,ji)
	                        AAB=0.d0
!multiplica por el coeficiente de la base
                             end do
                             VAABcuadrada(i,j) = VAABcuadrada(i,j) + acum*c(i,ii)*4*pi*exp(-Distcoef*a(i,ii))
!multiplica por el otro coef de la base

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             if (i .ge. j) then
!este if hay q sacarlo al fina, cuando cambie el rango en el q barre j , en vez de comenzar en 1 comience en i
                                pos=i+(1-j)*(j-2*M)/2   !chekeada
                                VAAB(pos) = VAAB(pos) + acum*c(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii)) !esta linea es lo unico que quedaria!!!!!!!!!!!!!
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
	end subroutine intECPAAB




        DOUBLE PRECISION function AAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
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
        do l = 0 , Lmax(z)-1
!barre todos los l de la parte no local
	   do term=1, expnumbersECP(z,l)
!barre contracciones del ECP para el atomo con carga z y l del ecp
	      Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
	      Qnl=0.d0
	      call Qtype1(Kmod,Ccoef,lmaxbase,necp(Z,l,term))
	      do lx=0,lxj
	      do ly=0,lyj
	      do lz=0,lzj
	         acumint=0.d0
	         do lambda=lxj+lyj+lzj,0,-1
	            acumang=0.d0
	            do m=-l,l
	               acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
	            end do
	            acumint=acumint+acumang*Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda)*aECP(z,L,term) 
	            if (Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda) .eq. 0.d0)  stop " q = 0 "
	         end do

	         AABz=AABz+acumint * dz**(lzj-lz) * comb(lzj,lz)
	         acumint=0.d0
	      end do
	      AABy=AABy+AABz * dy**(lyj-ly) * comb(lyj,ly)
	      AABz=0.d0
	      end do
	      AABx=AABx+AABy * dx**(lxj-lx) * comb(lxj,lx)
	      AABy=0.d0
	      end do
	      AAB_SEMILOCAL=AAB_SEMILOCAL+AABx
	      AABx=0.d0

	   end do
        end do
	return
	end function AAB_SEMILOCAL

        DOUBLE PRECISION function AAB_LOCAL(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz)
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
        integer :: w,n,z,l,lxi, lyi, lzi, lambda, Lmaxbase
        integer, intent(in) :: i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi
        double precision :: Ccoef, acum,dx,dy,dz,integral, Kmod
	double precision, dimension (3) :: Kvector

!variables temporales
	double precision :: ome, qrad


        Z=ZlistECP(k)
        L=Lmax(Z)
        Lmaxbase=lx+ly+lz+kxi+kyi+kzi
        AAB_LOCAL=0.d0
	acum=0.d0
	integral=0.d0
	Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
	Kmod= 2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)

        do w =1, expnumbersECP(z,l)
!barre todos los terminos del Lmaximo
	   Qnl=0.d0
	   Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
	   call Qtype1(Kmod,Ccoef,lmaxbase,necp(Z,l,w))

	   do lxi=0,lx
	   do lyi=0,ly
	   do lzi=0,lz
	      do lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
	         integral=integral + OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) * Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
!para debugueo
		ome=OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi)
		qrad=Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
!		write(*,*) i,j,ome,qrad,lambda,lxi+kxi,lyi+kyi,lzi+kzi
!		write(*,*) "KVECTOR", Kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      end do
	      acum= acum + integral*dx**(lx-lxi) * dy**(ly-lyi) * dz**(lz-lzi) *comb(lx,lxi) *comb(ly,lyi) * comb(lz,lzi)
	      integral=0.d0
	   end do
	   end do
	   end do
	   AAB_LOCAL=AAB_LOCAL+aECP(z,L,w)*acum
	   acum=0.d0
!aca multiplica x el coef
	end do
        return
        end function AAB_LOCAL




!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!terminos de 3centros <Xb|Va|Xc>, incluye el caso de 2 centros <Xb|Va|Xb>



	subroutine intECPABC()
	use garcha_mod, only : nshell,nuc,ncont,natom,a,c
	use ECP_mod, only : pi,ecptypes,cutECP,cutecp3,Lxyz,IzECP,VBACcuadrada,VBAC,local_nonlocal,ecp_debug,distx,disty,distz
	implicit none
	integer :: i,j,ii,ji,M,k,ki
	Double Precision :: Distcoef,dxi,dxj,dyi,dyj,dzi,dzj,ABC,acum,pos
	integer :: lxi,lxj,lyi,lyj,lzi,lzj
	M=nshell(0)+nshell(1)+nshell(2)
	acum=0.d0
	ABC=0.d0
	do i=1,M
	do j=1,M
!barre la base
	   do k=1, natom 


!		write(*,*) "calculando abc para", k

!necesito que este do barra por todos los nucleos del sistema
	      if (nuc(i) .ne. k .and. nuc(j) .ne.k) then
!solo calcula los terminos en q las 2 funciones de base NO corres ponden al atomo con el ECP
	         do ki=1, ecptypes
!barre atomos con ecp
	            if ( IzECP(k) .eq. ZlistECP(ki)) then
!solo calcula si el nucleo tiene ecp


!			write(*,*) "k,ki,izecp,zlist",k,ki,IzECP(k),ZlistECP(ki)

                       dxi=-distx(nuc(i),k)
                       dyi=-disty(nuc(i),k)
                       dzi=-distz(nuc(i),k)
                       dxj=-distx(nuc(j),k)
                       dyj=-disty(nuc(j),k)
                       dzj=-distz(nuc(j),k)
!distancias
                       lxi=Lxyz(i,1)
                       lxj=Lxyz(j,1)
                       lyi=Lxyz(i,2)
                       lyj=Lxyz(j,2)
                       lzi=Lxyz(i,3)
                       lzj=Lxyz(j,3)

	               do ii=1, ncont(i)
	               do ji=1, ncont(j)
!barre contracciones de la base
	                  Distcoef=a(i,ii)*(dxi**2.d0 + dyi**2.d0 + dzi**2.d0) + a(j,ji)*(dxj**2.d0 + dyj**2.d0 + dzj**2.d0)

			 if (Distcoef .ge. cutecp3) write(*,*) "cut 3 sobrepasado omit int"
                          if ( .not. cutECP .or. (Distcoef .lt. cutecp3)) then
!solo calcula los terminos que luego se multipliquen por un factor q no sea demasiado pequeño,
!el cut va a ser obligatorio, ya que si los atomos quedan muy separados el termino de Fock se obtiene como el producto de un numero MUY grande por uno que es casi 0
!y se obtienen NAN para el numero muy grande
!a distancias grades gana el termino que es casi 0

	                     ABC=ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)+4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)


!				if ( i .eq. j) write(*,*) "dist coef",Distcoef
!				if ( i .eq. j) write(*,*) "i,j,abc",i,j,abc

                                if (local_nonlocal .eq. 1 .and. ecp_debug) then
! local_nonlocal 1 y 2 son solo para debugueo
                                   ABC=ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
                                else if (local_nonlocal .eq. 2 .and. ecp_debug) then
                                   ABC=4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
                                end if
		                        acum=acum + ABC*c(j,ji)*exp(-Distcoef)
!			if ( i .eq. j) write(*,*) "i,j,abc",i,j,abc
!			if ( i .eq. j) write(*,*) "ii,ji,k", ii,ji,k
!			if ( i .eq. j) write(*,*) "lxi,lyi,lzi,lxj,lyj,lzj",lxi,lyi,lzi,lxj,lyj,lzj
!			if ( i .eq. j) write(*,*) "-dxi,-dyi,-dzi,-dxj,-dyj,-dz",-dxi,-dyi,-dzi,-dxj,-dyj,-dzj

	                     ABC=0.d0
	                  end if
	               end do

	               VBACcuadrada(i,j)= VBACcuadrada(i,j) + acum*c(i,ii)*4.d0*pi

	               if (i .ge. j) then
!este if hay q sacarlo al fina, cuando cambie el rango en el q barre j , en vez de comenzar en 1 comience en i
	                  pos=i+(1-j)*(j-2*M)/2   !chekeada
	                  VBAC(pos) = VBAC(pos) + acum*c(i,ii)*4.d0*pi !esta linea es lo unico que quedaria!!!!!!!!!!!!!
!			  write(*,*) i,j,pos,VBAC(pos)
	               end if
                       acum=0.d0
	               end do
	            end if 
	         end do
	      end if
	   end do
	end do
	end do
	end subroutine intECPABC


	DOUBLE PRECISION function ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2)
	use garcha_mod, only : a,c
	use ECP_mod, only :expnumbersECP, Qnl,bECP,IzECP
	implicit none
        integer, intent(in) :: i,j,ii,ji
!terminos de la base
	integer, intent(in) :: lxi,lyi,lzi,lxj,lyj,lzj
!potencias de la parte angular
	double precision, intent(in) :: dx1,dy1,dz1,dx2,dy2,dz2
!distancias del centro con ecp a cada nucleo
        integer, intent(in) :: k
!numero de atomo,
        integer :: L, Z
!L maximo del ecp, carga sin modificar
        double precision, dimension (3) :: Kvector
	double precision :: Kmod,Ccoef, integral,auxcomb,auxdist,acum
	integer :: lmaxbase

	integer :: ac,bc,cc,dc,ec,fc,w,lambda
!variables auxiliares


!temporales para debug
	double precision :: radial1,omegal


	ABC_LOCAL=0.d0
	z=IzECP(k)
	L=Lmax(Z)
	lmaxbase=lxi+lyi+lzi+lxj+lyj+lzj
	Kvector=(/a(i,ii)*dx1+a(j,ji)*dx2,a(i,ii)*dy1+a(j,ji)*dy2,a(i,ii)*dz1+a(j,ji)*dz2/)
	Kvector=-2.d0*Kvector
	Kmod=sqrt(Kvector(1)**2.d0+Kvector(2)**2.d0+Kvector(3)**2.d0)
	integral=0.d0
	acum=0.d0

!	write(*,*) "testeando local"
!	write(*,*) "kvector",Kvector
!	if (Kvector(1) .ge. 0.d0) Kvector(1)=-Kvector(1)

	do w =1, expnumbersECP(z,l)
!barre terminos del ECP para el atomo con carga nuclear Z y l del ECP
	   Qnl=0.d0
           Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
		if (i .eq. j) then
!			write(*,*) "llamando a Q", i, j
!			write(*,*) "b,a1,a2",bECP(z,L,w),a(i,ii),a(j,ji)
			write(*,*) "becp,z,l,w",bECP(z,L,w),z,L,w
!			write(*,*) Kmod,Ccoef,lmaxbase,necp(Z,l,w)
		end if
	   call Qtype1(Kmod,Ccoef,lmaxbase,necp(Z,l,w))
!calcula integral radial
	   do ac=0,lxi
	   do bc=0,lyi
	   do cc=0,lzi
	   do dc=0,lxj
	   do ec=0,lyj
	   do fc=0,lzj
!barre los coef de la expansion en el binomio de Newton
              do lambda=ac+bc+cc+dc+ec+fc,0,-2
                 integral=integral + OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc) * Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)

		if (i .eq. j) then
!			write(*,*) "*************************************"
!			write(*,*) "omega, int",OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc),Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)
!	                write(*,*) "*************************************"
		endif

!		if ( i .eq. j) write(*,*) "omega,Q",OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc),Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)
		if (Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda) .ne. Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)) stop "Q-1314=0"
              end do
	      auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
	      auxdist=dx1**(lxi-ac)*dy1**(lyi-bc)*dz1**(lzi-cc)*dx2**(lxj-dc)*dy2**(lyj-ec)*dz2**(lzj-fc)
	      acum=acum + auxcomb*auxdist*integral
		if (i .eq. j) then
!		write(*,*) "auxcomb*auxdist*integral",auxcomb,auxdist,integral
		end if
	      integral=0.d0
	   end do
           end do
           end do
           end do
           end do
           end do
           ABC_LOCAL=ABC_LOCAL+aECP(z,L,w)*acum
           acum=0.d0
	end do
	return
	end function ABC_LOCAL

	DOUBLE PRECISION function ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
        use garcha_mod, only : a,c
        use ECP_mod, only : Qnl1l2,necp,bECP,IzECP
	implicit none
	integer, intent(in) :: i,j,ii,ji,k
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
	integer, intent(in) :: lxi,lyi,lzi,lxj,lyj,lzj
!potencias de la parte angular
	double precision, intent(in) :: dxi,dyi,dzi,dxj,dyj,dzj
!distancias de las nucleos con bases al nucleo con ecp
        integer :: Z,l1max,l2max
!Z= carga del nucleo, y momento angular de la base i y j
        double precision, dimension (3) :: Kivector,Kjvector
        double precision :: Kimod, Kjmod,Ccoef
	integer :: l,term,ac,bc,cc,dc,ec,fc,lambdai, lambdaj,m
!auxiliares ciclos
	double precision :: acumang,integral,auxcomb,auxdist,acum
!auxiliaresa


 

	ABC_SEMILOCAL=0.d0
	integral=0.d0
	acum=0.d0
	Qnl1l2=0.d0
	acumang=0.d0
	auxcomb=0.d0
	auxdist=0.d0

        Z=IzECP(k)
	l1max=lxi+lyi+lzi
	l2max=lxj+lyj+lzj

        Kivector=-2.d0*a(i,ii)*(/dxi,dyi,dzi/)
        Kjvector=-2.d0*a(j,ji)*(/dxj,dyj,dzj/)

        Kimod=sqrt(Kivector(1)**2.d0+Kivector(2)**2.d0+Kivector(3)**2.d0)
	Kjmod=sqrt(Kjvector(1)**2.d0+Kjvector(2)**2.d0+Kjvector(3)**2.d0)


        do l = 0 , Lmax(z)-1
!barre todos los l de la parte no local
           do term=1, expnumbersECP(z,l)
!barre contracciones del ECP para el atomo con carga z y l del ecp
	      Qnl1l2=0.d0
	      Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
	      call Qtype2(Kimod,Kjmod,Ccoef,l1max,l2max,necp(Z,l,term))


!agrega a la matriz Qnl1l2 los terminos correspondientes a un termino radiales.
	      do ac=0,lxi
	      do bc=0,lyi
	      do cc=0,lzi
	      do dc=0,lxj
	      do ec=0,lyj
	      do fc=0,lzj

	         do lambdai=ac+bc+cc,0,-2
	         do lambdaj=dc+ec+fc,0,-2
                    acumang=0.d0
                    do m=-l,l
                       acumang=acumang+OMEGA2(Kivector,lambdai,l,m,ac,bc,cc)*OMEGA2(Kjvector,lambdaj,l,m,dc,ec,fc)
                    end do
	            integral=integral+acumang*Qnl1l2(ac+bc+cc+dc+ec+fc+necp(Z,l,term),lambdai,lambdaj)
                    if (Qnl1l2(ac+bc+cc+dc+ec+fc+necp(Z,l,term),lambdai,lambdaj) .eq. 0.d0) stop " qnl1l2 = 0 in ABC_SEMILOCAL"
	         end do
	         end do
                 auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
                 auxdist=dxi**(lxi-ac)*dyi**(lyi-bc)*dzi**(lzi-cc)*dxj**(lxj-dc)*dyj**(lyj-ec)*dzj**(lzj-fc)
                 acum=acum + auxcomb*auxdist*integral

                 integral=0.d0

	      end do
	      end do
	      end do
	      end do
	      end do
	      end do
              ABC_SEMILOCAL=ABC_SEMILOCAL+aECP(z,L,term)*acum
              acum=0.d0
	   end do
	end do
	return
	end function ABC_SEMILOCAL


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!funciones para integrales radiales


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
           if ( (n-1)/2 .lt. 0) stop "Er factorial de un negativo en Q0"
           Q0=fac((n-1)/2)/(2*alpha**((n+1)/2))
           return
        end if
        end function Q0



!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! funciones para integrales angulares

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
!angularint(nx,ny,nz) calcula int ((x/r)^nx (y/r)^ny (z/r)^nz d(angulo solido)
           end do
        end do
        return
        end function Aintegral


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
        if ( all(K .eq. (/0.d0,0.d0,0.d0/))) return
!caso especial para los terminos <A|A|A>

        Kun=K/sqrt(K(1)**2.d0 + K(2)**2.d0 + K(3)**2.d0)
        do u=-l,l
           do r=0,l
              do s=0,l-r
                 t=l-r-s
                 SUM1=SUM1+Ucoef(l,u,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
                 SUM2=SUM2+Ucoef(l,u,r,s,t)*angularint(a+r,b+s,c+t)
!Usa la funcion angular int en lugar del array int angular. Hay que agrandar el array
              end do
           end do
           OMEGA1=OMEGA1+SUM1*SUM2
           SUM1=0.d0
           SUM2=0.d0
        end do
        return
        end function OMEGA1


        DOUBLE PRECISION FUNCTION OMEGA2(K,lambda,l,m,a,b,c)
!Esta funcion devuelve el valor de omega2 evaluado en el vector K
!OMEGA2(K,lambda,l,m,a,b,c) = Suma(o=-lambda a o=lambda) ( Y lambda o(K)*  int (x/r)^a * (y/r)^b * (z/r)^c Ylambda o() Ylm()  d angulo solido )
        use ECP_mod, only : intangular,angularint
        implicit none
        DOUBLE PRECISION, intent(in), dimension(3) :: K
        integer, intent(in) :: lambda,l,m,a,b,c
        DOUBLE PRECISION,dimension(3) :: Kun
        integer :: o,r,s,t,u,v,w
        Double precision :: SUM1, SUM2
        Kun=K/sqrt(K(1)**2.d0 + K(2)**2.d0 + K(3)**2.d0)
        SUM1=0.d0
        SUM2=0.d0
        OMEGA2=0.d0
        do o=-lambda,lambda
           do r=0,lambda
              do s=0,lambda-r
                 t=lambda-r-s
                 SUM1=SUM1+Ucoef(lambda,o,r,s,t)*Kun(1)**r * Kun(2)**s * Kun(3)**t
                 do u=0,l
                    do v = 0,l-u
                       w=l-u-v
                       SUM2=SUM2+Ucoef(lambda,o,r,s,t)*Ucoef(l,m,u,v,w)*angularint(a+r+u,b+s+v,c+t+w)
                    end do
                 end do
              end do
           end do
           OMEGA2=OMEGA2+SUM1*SUM2
           SUM1=0.d0
           SUM2=0.d0
        end do
        end function OMEGA2



        DOUBLE PRECISION FUNCTION Ucoef(l,m,lx,ly,lz)
!lz esta al pedo, puedo usarlo para chekear que lx+ly+lz=l
        use ECP_mod, only :  l0,l1,l2,l3,l4
        implicit none
        integer, intent(in) :: l,m,lx,ly,lz

        if (lx+ly+lz .ne. l) then
           stop "problem with Ucoef: lx+ly+lz not equal l"
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
           stop "ECP error l is greater than 4"
        end if
        return
        end function Ucoef



!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! funciones para calculos auxiliares

        double precision function comb(a,b)
!devuelve el combinatorio a,b; a>=b
!en el futuro puedo pasarlo a un array
        use ECP_mod, only :fac
        integer, intent(in) :: a,b
        if (a .lt. 0 .or. b .lt. 0) stop "comb con numeros negativos"
        if (a .lt. b) stop "b mayor que a en comb"
        comb=0.d0
        comb=fac(a)/((fac(b)*fac(a-b)))
        return
        end function comb


        end subroutine intECP





!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subrutinas para integrales Radiales

        subroutine Qtype1(K,Ccoef,lmax,necp)
!this routine obtain Q(l,n,k,a) where
!Q(l,n,k,a)= int r^n exp(-ar^2) Ml(kr) dr form r=0 to inf
!where Ml are the modified spherical Bessel function of the first kind.

!cuidado esta mal el rango en que se calcula Q, hay q chequear
!agrega a la matriz Qnl los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl ya que hay q llamar a esta rutina por cada termino del pseudopotencial
!luego hay que optimizar la rutina para que calcule Qnl a partir de otros Qnl
	use ECP_mod, only :  alpha, betha, Bn, Cn, Qnl, ecp_full_range_int
	implicit none
	double precision, intent(in) :: K,Ccoef
!lmax  = 0 para <s||s>, 1 para <s||p>, 2 para <s||d>, ... , 6 para <d||d>	
!necp corresponde al exponente del paseudopotencial en r^n
!	integer, intent(in) :: necp,lmax
	integer :: necp,lmax
!nmin y nmax dan el rango en que tiene q calcular Q
	integer :: nmin,nmax
!variables auxiliares
	integer :: n,l,i
        double precision :: acoef, gam, acum
	nmin=necp
	nmax=necp+lmax
	
	if (ecp_full_range_int) then
!solucion temporal hasta ver el rango correcto de calculo
	   nmin=0
	   nmax=10
	   lmax=4
	end if

	acoef=K/(2.d0*Ccoef)
	gam=0.5d0*exp(K**2.d0/(4.d0*Ccoef))
	Bn=0.d0
	Cn=0.d0
	call ByC(Acoef,Ccoef,nmin-1,nmax,Bn,Cn)

	do n=nmin,nmax
	do l=0,lmax

	acum=0.d0
	do i=l,1,-2
	   acum=acum+alpha(l,i)*Bn(n-i)/k**i
	end do
	do i=l+1,1,-2
	   acum=acum+betha(l+1,i)*Cn(n-i)/k**i
	end do
	if (acum .eq. 0.d0) then
	   write(*,*) "error en Qtype 1, integral radial 0. n,l",n,l
	   stop
	end if 
		Qnl(n,l)=Qnl(n,l)+acum*gam
		acum=0.d0
	end do
	end do
        End subroutine Qtype1


        subroutine ByC(Acoef,Ccoef,nmin,nmax,Barray,Carray)
!calcula los coeficientes B y C 
	use ECP_mod, only : DAW,DAWERF,NEXTCOEF,pi,pi12
	IMPLICIT NONE
!acoef,ccoef son los coeficientes para el calculo de B y C
! A(ccoef,acoef)= int exp(-ccoef*x^2)(x+acoef)^n dx from -acoef to inf
!nmin,nmax delimitan el rango de n
!Bn=An(ccoef,acoef) + An(ccoef,-acoef)
!Cn=An(ccoef,acoef) - An(ccoef,-acoef)
	DOUBLE PRECISION, intent(in) :: Acoef,Ccoef
	integer, intent(in) :: nmin,nmax
	DOUBLE PRECISION, intent(inout), dimension (-12:14) :: Barray,Carray
!variables auxiliares
	DOUBLE PRECISION :: C0sq,ncos,ca
	integer :: i
!	write(*,*) "BYC", acoef,ccoef,nmin,nmax!,Barray,Carray

	C0sq=sqrt(Ccoef)
	Barray(0)=pi12/c0sq
	Carray(0)=Barray(0)*erf(Acoef*C0sq)

	if (nmax>0) then
	   Barray(1)= exp(-Ccoef*Acoef**2)/Ccoef + Acoef*Carray(0)
	   Carray(1)= Barray(0)*Acoef
	   do i=2,nmax
	      ncos=(i-1)/(2*Ccoef)
	      Barray(i)=ncos*Barray(i-2)+Acoef*Carray(i-1)
	      Carray(i)=ncos*Carray(i-2)+Acoef*Barray(i-1)
	   end do
	end if

        if (nmin<0) then
	   Barray(-1)=2*pi12*DAWERF(Acoef*C0sq)
	   Carray(-1)=2*pi12*DAW(Acoef*C0sq)
	   if (nmin<-1) then
	      ca=2*Ccoef*Acoef
	      Barray(-2)=ca*Carray(-1)-2*Ccoef*Barray(0)
	      Carray(-2)=2*ca*exp(-Ccoef*Acoef**2)+ca*Barray(-1)-2*Ccoef*Carray(0)
	      do i=-3,nmin,-1
	         Barray(i)=NEXTCOEF(1,i,ca,exp(-Ccoef*Acoef**2),Ccoef,Carray(i+1), Barray(i+2))
	         Carray(i)=NEXTCOEF(-1,i,ca,exp(-Ccoef*Acoef**2),Ccoef,Barray(i+1), Carray(i+2))
	      end do
	   end if
	end if
        end subroutine ByC



        subroutine Qtype2(Ka,Kb,Ccoef,l1max,l2max,necp)
!this routine obtain Q(l1,l2,n,k1,k2,a) where
!Q(l1,l2,n,k1,k2,a)= int r^n exp(-ar^2) Ml1(k1r) Ml2(k2r) dr form r=0 to inf
!where Ml are the modified spherical Bessel function of the first kind.

!agrega a la matriz Qnl1l2 los terminos correspondientes a un termino del pseudopotencial.
!CUIDADO no borra Qnl1l2 ya que hay q llamar a esta rutina por cada termino del pseudopotencial
!luego hay que optimizar la rutina para que calcule Qnl1l2 a partir de otros Qnl1l2
	use ECP_mod, only :  alpha, betha, rho, tau, sigma, sigmaR, Qnl1l2,ecp_full_range_int
	implicit none
	double precision, intent(in) :: Ka,Kb,Ccoef
!l1max y l2max = 0 para s, 1 para p, 2 para d, etc	
!n corresponde al exponente del paseudopotencial en r^n
!	integer, intent(in) :: necp,l1max, l2max
	integer :: necp,l1max, l2max
	integer :: nmin,nmax
!variables auxiliares
	integer :: i,j,n,l1,l2
	double precision :: alfok, betok, acum1
	
!Caso 1-calcula todas las integrales
	nmin=necp
	nmax=necp+l1max+l2max

	if (ecp_full_range_int) then
!solucion temporal hasta ver el rango correcto de calculo
	   nmin=0
	   nmax=10
	   l1max=4
	   l2max=4
	end if

	call integrals(Ka,Kb,Ccoef,necp-2-l1max-l2max,necp+l1max+l2max-2)
	acum1=0.d0
	do n=nmin,nmax
	do l1=0,l1max
	do l2=0,l2max
	   do i=l1,1,-2
	      alfok=alpha(l1,i)/Ka**i
	      do j=l2,1,-2
	         if (tau(n-i-j) == 0.d0) stop "Error, no calculo tau"
	         acum1=acum1+alpha(l2,j)*tau(n-i-j)/Kb**j
	      end do
	      do j=l2+1,1,-2
	         acum1=acum1+betha(l2+1,j)*sigmaR(n-i-j)/Kb**j
		 if (sigmaR(n-i-j) == 0.d0) stop "Error, no calculo sigmaR"
	      end do
	   Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*alfok
	   acum1=0.d0
	   end do
	   do i=l1+1,1,-2
	      betok=betha(l1+1,i)/Ka**i
	      do j=l2,1,-2
	         acum1=acum1+alpha(l2,j)*sigma(n-i-j)/Kb**j
	         if (sigma(n-i-j) == 0.d0) stop "Error, no calculo sigma"
	      end do
	      do j=l2+1,1,-2
	         acum1=acum1+betha(l2+1,j)*rho(n-i-j)/Kb**j
	         if (rho(n-i-j) == 0.d0) stop "Error, no calculo rho"
	      end do
	      Qnl1l2(n,l1,l2)=Qnl1l2(n,l1,l2)+acum1*betok
	      acum1=0.d0
	   end do
	end do
	end do
	end do
        End subroutine Qtype2

	subroutine integrals(Ka,Kb,Ccoef,nmin,nmax)
!obtain integrals rho, tau, sigma and sigmaR from n between nmin and nmax
! rho(n) = int exp(-cr^2) * sinh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! sigma(n) = int exp(-cr^2) * sinh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf
! sigmaR(n) = int exp(-cr^2) * cosh(Ka*r)* sinh(Kb*r) r^n  dr  from 0 to inf
! tau(n) = int exp(-cr^2) * cosh(Ka*r)* cosh(Kb*r) r^n  dr  from 0 to inf

	use ECP_mod, only : Bn1,Bn2,Cn1,Cn2,rho, tau, sigma, sigmaR,ecp_full_range_int
	implicit none
!	integer, intent(in) :: nmin,nmax
	integer :: nmin,nmax
	double precision, intent(in) :: Ka,Kb,Ccoef
	double precision, dimension(2) :: acoef,gammacoef
	double precision :: signo
	acoef(1)=0.5d0*(Ka+Kb)/Ccoef
	acoef(2)=0.5d0*abs(Ka-Kb)/Ccoef
	gammacoef(1)=0.25d0*exp(Ccoef*acoef(1)**2.d0)
	gammacoef(2)=0.25d0*exp(Ccoef*acoef(2)**2.d0)
!	call ByCdoble(acoef(1),acoef(2),Ccoef,nmin,nmax)
	Bn1=0.d0
	Bn2=0.d0
	Cn1=0.d0
	Cn2=0.d0
	rho=0.d0
	tau=0.d0
	sigma=0.d0
	sigmaR=0.d0

        if (ecp_full_range_int) then
!solucion temporal hasta ver el rango correcto de calculo
           nmin=-12
           nmax=14
        end if
	call ByC(acoef(1),Ccoef,nmin,nmax,Bn1,Cn1)
	call ByC(acoef(2),Ccoef,nmin,nmax,Bn2,Cn2)
	rho=gammacoef(1)*Bn1-gammacoef(2)*Bn2
	tau=gammacoef(1)*Bn1+gammacoef(2)*Bn2
	sigma=gammacoef(1)*Cn1+sign(1.d0,Ka-Kb)*gammacoef(2)*Cn2
	sigmaR=gammacoef(1)*Cn1-sign(1.d0,Ka-Kb)*gammacoef(2)*Cn2
	end subroutine integrals


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!rutinas de preparacion y calculo de variables previas

        subroutine allocate_ECP
!alocatea todas las variables que van a necesitar los pseudopotenciales
        use garcha_mod, only :nshell, natom
        use ECP_mod, only :VAAAcuadrada,VAABcuadrada, VBACcuadrada,VAAA,VAAB,VBAC,term1e,distx, disty, distz, IzECP,Lxyz,Cnorm,VXXXgamess
!	use param, only :nl
!term1e terminos de fock de 1 electron sin agregarles VAAA
        implicit none
        integer :: ns,np,nd,M,Mcuad
        ns=nshell(0)
        np=nshell(1)
        nd=nshell(2)
        M=ns+np+nd
        allocate (VAAAcuadrada(M,M),VAABcuadrada(M,M), VBACcuadrada(M,M),VXXXgamess(M,M))
        VAAAcuadrada=0.d0
        VAABcuadrada=0.d0
        VBACcuadrada=0.d0
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
	allocate (IzECP(natom))
	allocate (Lxyz(M, 3))
!	allocate (Cnorm(ngnu,nl))
        end subroutine allocate_ECP

        subroutine deallocateV
!desalocatea variables de ECP
        use ECP_mod, only :VAAAcuadrada,VAABcuadrada, VBACcuadrada,VAAA,VAAB,VBAC,term1e,distx, disty, distz,IzECP,Lxyz
        implicit none
        deallocate (VAAAcuadrada,VAABcuadrada, VBACcuadrada)
        deallocate (VAAA,VAAB,VBAC)
        deallocate (term1e)
        deallocate (distx, disty, distz)
	deallocate (IzECP,Lxyz)
        end subroutine deallocateV

	subroutine norm_C()
	use ECP_mod, only :Cnorm
	implicit none
	end subroutine norm_C

        subroutine correctbasis(r)
!cambia el coeficiente de las funciones d x^2, y^2 y z^2
        use garcha_mod, only :nshell,c,ncont
        implicit none
        integer :: i,j,ns,np,nd
        integer, intent(in) :: r
        double precision :: factor
        ns=nshell(0)
        np=nshell(1)
        nd=nshell(2)

        if (r.eq.1) then
                factor=sqrt(1.d0/3.d0)
        else if (r.eq.-1)then
                factor=sqrt(3.d0)
        end if

	write(*,*) " EDITANDO BASE D *****************************************"
        do i=ns+np+1,ns+np+nd,6
        write(*,*) "edita i=",i
                do j=1,ncont(i)
                        c(i,j)=c(i,j)*factor
                end do
                do j=1,ncont(i+2)
                        c(i+2,j)=c(i+2,j)*factor
                end do
                do j=1,ncont(i+5)
                        c(i+5,j)=c(i+5,j)*factor
                end do
        end do
        end subroutine correctbasis


        subroutine ReasignZ
!cambia la carga de los nucleos con pseudopotenciales  sacandole la carga del core y guarda las cargas originales en IzECP
!tambien corrige la cantidad de electrones sacando los que van al core
        use garcha_mod, only :Iz,nuc,nshell, natom, NCO
        use ECP_mod, only : ZlistECP,IzECP,Zcore,ecptypes
        implicit none
        integer :: i,j
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
          end if
         end do
        end do
        end subroutine ReasignZ


        subroutine obtainls
!arma una matriz que contenga los exponentes de la parte angular de la base
! |x> = A x^lx y^ly z^lz *e^-ar^2
! angularL(i,j) j=1 => lx, j=2 => ly, j=3 => lz para la funcion i de la base
        use garcha_mod, only : nshell
        use ECP_mod, only : Lxyz
        implicit none
        integer :: i,resto
        integer :: M,lx,ly,lz
        M = nshell(0)+nshell(1)+nshell(2)
        Lxyz=0
        do i=nshell(0)+1,M
           if (i .le. nshell(0)+nshell(1)) then
!funciones p
              resto=modulo(i-nshell(0),3)
              if (resto .eq. 1) Lxyz(i,1)=1
              if (resto .eq. 2) Lxyz(i,2)=1
              if (resto .eq. 0) Lxyz(i,3)=1
           else if (i .le. M) then
!funciones d
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
           end if
        end do
        end subroutine obtainls

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




!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! rutinas para verbose_ECP

	subroutine WRITE_BASIS()
!luego agregar para q escriba el atomo, y el momento angular, y que use los coeficiente d bien normalizados
	use garcha_mod, only : nshell,a,c,ncont
	implicit none
	integer :: i,j,M
	M=nshell(0)+nshell(1)+nshell(2)
!escribe a y c (solo el 1er valor)
        write(*,4010)
	write(*,4011)
	write(*,4012)
	write(*,4013)
	write(*,4014)
        do i=1,M
                do j=1,ncont(i)
                        write(*,4015) i,j,a(i,j),c(i,j)
                end do
        end do
	write(*,4016)

	4010 format("╔═══════════════════════════════════════════╗")
	4011 format("║            NORMALIZED BASIS SET           ║")
	4012 format("╠═══════╦═══════╦════════════╦══════════════╣")
	4013 format("║ Basis ║ Cont. ║  Exponent  ║ Coefficient  ║")
	4014 format("╠═══════╬═══════╬════════════╬══════════════╣")
	4015 format("║",1x,i3,3x,"║",1x,i3,3x,"║",1x,f10.7,1x,"║",1x,f12.9,1x,"║")
	4016 format("╚═══════╩═══════╩════════════╩══════════════╝ ")
	end subroutine WRITE_BASIS


	subroutine WRITE_FOCK_ECP_TERMS()
!escribe los terminos de fock del pseudopotencial
	use ECP_mod, only : VAAA, VAAB, VBAC
	use garcha_mod, only : nshell
	implicit none
	integer :: i,j,k,M
	logical :: no_zero
	no_zero = .true.
	M=nshell(0)+nshell(1)+nshell(2)
	write(*,4020)
        write(*,4021)
        write(*,4022)
        write(*,4023)
        write(*,4024)
	do i=1,M
	do j=1,M
	   if (i .ge. j) then
	      k=i+(1-j)*(j-2*M)/2
	         if (no_zero) then
	            if ( abs(VAAA(k)+VAAB(k)+VBAC(k)) .gt. 1E-30) then
	               write(*,4025) i,j,VAAA(k),VAAB(k),VBAC(k)
	            end if
	         else
	            write(*,4025) i,j,VAAA(k),VAAB(k),VBAC(k)
	         end if
            end if
	end do
	end do
	write(*,4026)

	4020 format("╔═══════════════════════════════════════════════════════"&
		,"════════════╗") 
	4021 format("║                       FOCK Pseudopotencials           "&
		,"            ║")
	4022 format("╠═══╦═══╦═══════════════════╦═══════════════════╦═══════"&
		,"════════════╣")
	4023 format("║ i ║ j ║      <A|A|A>      ║      <A|A|B>      ║      <"&
		,"B|A|C>      ║")
	4024 format("╠═══╬═══╬═══════════════════╬═══════════════════╬═══════"&
		,"════════════╣")
	4025 format("║",i2,1x,"║",i2,1x,"║",f18.15,1x,"║",f18.15,1x,"║",f18.15,1x,"║",f18.15,1x,"║")
	4026 format("╚═══╩═══╩═══════════════════╩═══════════════════╩═══════"&
		,"════════════╝ ")
	end subroutine WRITE_FOCK_ECP_TERMS


	subroutine WRITE_FOCK_ECP()
	        use ECP_mod, only : VAAA, VAAB, VBAC
        use garcha_mod, only : nshell
        implicit none
        integer :: i,j,k,M
        logical :: no_zero
        no_zero = .true.
        M=nshell(0)+nshell(1)+nshell(2)
        write(*,4032)
        write(*,4033)
        write(*,4034)
        do i=1,M
        do j=1,M
           if (i .ge. j) then
              k=i+(1-j)*(j-2*M)/2
                 if (no_zero) then
                    if ( abs(VAAA(k)+VAAB(k)+VBAC(k)) .gt. 1E-30) then
                       write(*,4035) i,j,VAAA(k)+VAAB(k)+VBAC(k)
                    end if
                 else
                    write(*,4035) i,j,VAAA(k)+VAAB(k)+VBAC(k)
                 end if
            end if
        end do
        end do
        write(*,4036)

	4032 format(5x,"╔═══╦═══╦═══════════════════════╗")
	4033 format(5x,"║ i ║ j ║ FOCK Pseudopotencials ║")
	4034 format(5x,"╠═══╬═══╬═══════════════════════╣")
	4035 format(5x,"║",i2,1x,"║",i2,1x,"║",2x,f18.15,3x,"║")
	4036 format(5x,"╚═══╩═══╩═══════════════════════╝ ")

	end subroutine WRITE_FOCK_ECP


	subroutine WRITE_ANG_EXP()
!escribe el exponente de la parte angular de la funcion de base:
!xi(x,y,z)= x^nx * y^ny * z^nz *f(r)
	use ECP_mod, only :lxyz
	use garcha_mod, only : nshell
	implicit none
	integer :: i

	write(*,4040)
	write(*,4041)
	write(*,4042)
        write(*,4043)
        write(*,4044)

	do i=1,nshell(0)+nshell(1)+nshell(2)
	   write(*,4045) i,lxyz(i,1),lxyz(i,2),lxyz(i,3)
	end do

        write(*,4046)

	4040 format("╔═══════════════════╗")
	4041 format("║ Angular Exponents ║")
	4042 format("╠═══════╦═══╦═══╦═══╣")
	4043 format("║ Basis ║ x ║ y ║ z ║")
	4044 format("╠═══════╬═══╬═══╬═══╣")
	4045 format("║",2x,i2,3x,"║",i2,1x,"║",i2,1x,"║",i2,1x,"║")
	4046 format("╚═══════╩═══╩═══╩═══╝ ")
	end subroutine WRITE_ANG_EXP

	subroutine WRITE_DISTANCE
	use ECP_mod, only :distx,disty,distz
	use garcha_mod, only : natom
	implicit none
	integer :: i,j

	write(*,4050)
	write(*,4051)
	write(*,4052)
	write(*,4053)
	write(*,4054)

	do i=1,natom
	do j=1,natom
	   write(*,4055) i,j,distx(i,j),disty(i,j),distz(i,j)
	end do
	end do
	write(*,4056)

        4050 format("╔════════════════════════════════════════════════" &
		,"═════════════════════════════╗")
        4051 format("║                               Distances (Bohr) " &
		,"                             ║")
        4052 format("╠════════╦════════╦═══════════════════╦══════════" &
		,"═════════╦═══════════════════╣")
        4053 format("║ atom i ║ atom j ║     distance x    ║     dista" &
		,"nce y    ║     distance z    ║")
        4054 format("╠════════╬════════╬═══════════════════╬══════════" &
		,"═════════╬═══════════════════╣")
        4055 format("║",2x,i3,3x,"║",2x,i3,3x,"║",f18.15,1x,"║",f18.15,1x,"║",f18.15,1x,"║")
        4056 format("╚════════╩════════╩═══════════════════╩══════════" &
		,"═════════╩═══════════════════╝ ")

	end subroutine WRITE_DISTANCE

	subroutine WRITE_ECP_PARAMETERS()
	use ECP_mod, only : ecptypes,Lmax,ZlistECP,expnumbersECP,ZlistECP,nECP,bECP,aECP,asignacion,Zcore
	implicit none
	integer :: k,l,w,z
        character :: simb*3
	write(*,*)
	write(*,*)
	write(*,4060)
	write(*,4061)
	write(*,4062)

	do k=1, ecptypes
	z=ZlistECP(k)
	   call asignacion(Z,simb)
	   write(*,*)
	   write(*,4063)
	   write(*,4064) simb,Zcore(z),Lmax(z)
	   write(*,4065)
!barre todos los atomos con ECP
	   do l=0,Lmax(z)
	      write(*,4066) l
	      write(*,4067)
	      write(*,4068)
	      write(*,4069)
!barre momento angular del ecp
	         do w =1, expnumbersECP(z,l)
!barre funciones
	            write(*,4070) necp(Z,l,w),bECP(z,L,w),aECP(z,L,w)
	         end do
	      if ( l .lt. Lmax(z) ) write(*,4071)
	      if ( l .eq. Lmax(z) ) write(*,4072)
	   end do
	end do


        4060 format(4x,"╔═════════════════════════════════════════════╗")
        4061 format(4x,"║     EFFECTIVE CORE POTENTIAL PARAMETERS     ║")
	4062 format(4x,"╚═════════════════════════════════════════════╝")
	4063 format(8x,"╔════════╦════════════════╦═══════════╗")
        4064 format(8x,"║",2x,A3,3x,"║",2x,"ZCore=",i3,5x,"║",2x,"Lmax=",i1,3x,"║")
	4065 format(8x,"╠════════╩════════════════╩═══════════╣")
	4066 format(8x,"║   L=",i1,"                               ║") 
	4067 format(8x,"╠═══╦════════════════╦════════════════╣")
	4068 format(8x,"║ n ║    exponent    ║  coefficient   ║")
	4069 format(8x,"╠═══╬════════════════╬════════════════║")
	4070 format(8x,"║",1x,i1,1x,"║",f12.6,4x,"║",f12.6,4x,"║")
	4071 format(8x,"╠═══╩════════════════╩════════════════╣")
        4072 format(8x,"╚═══╩════════════════╩════════════════╝")

	end subroutine WRITE_ECP_PARAMETERS

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! rutinas para debug_ecp 

	subroutine SEARCH_NAN()
	use ECP_mod, only : VAAA,VAAB,VBAC
	use garcha_mod, only : nshell
	implicit none
	integer :: M,MM,i
	M=nshell(0)+nshell(1)+nshell(2)
	MM=M*(M+1)/2
	do i=1,MM
	   if (VAAA(i) .ne. VAAA(i)) stop "NAN in VAAA"
	   if (VAAB(i) .ne. VAAB(i)) stop "NAN in VAAB"
	   if (VBAC(i) .ne. VBAC(i)) stop "NAN in VBAC"
	end do
	end subroutine SEARCH_NAN


	subroutine TEST_SYMMETRY()
        use ECP_mod, only : VAAAcuadrada,VAABcuadrada,VBACcuadrada
        use garcha_mod, only : nshell,nuc
	implicit none
	integer :: i,j,M
	logical :: Er
	Er=.false.
	M=nshell(0)+nshell(1)+nshell(2)
	do i=1,M
	do j=i,M

	   if (abs(VAAAcuadrada(i,j) - VAAAcuadrada(j,i)) .gt. 1D-08) then
			write(*,*) "symmetry error in ecp-fock matrix for <xA|VA|xA> terms"
			write(*,4031) i,j,VAAAcuadrada(i,j), VAAAcuadrada(j,i)
			Er=.true.
	   end if

	   if (abs(VAABcuadrada(i,j) - VAABcuadrada(j,i)) .gt. 1D-08) then
		write(*,*) "symmetry error in ecp-fock matrix for <xA|VA|xB> terms",i,j
		write(*,4031) i,j,VAABcuadrada(i,j), VAABcuadrada(j,i)
		Er=.true.
	   end if

	   if (abs(VBACcuadrada(i,j) - VBACcuadrada(j,i)) .gt. 1D-08) then
		write(*,*) "symmetry error in ecp-fock matrix for <xB|VA|xC> terms",i,j
		write(*,4031) i,j,VBACcuadrada(i,j), VBACcuadrada(j,i)
		Er=.true.
	   end if

	end do
	end do

	if ( Er ) then
           open(UNIT=92,FILE="ECP-FOCK",STATUS="UNKNOWN",ACTION="WRITE")
	   write(92,4032)
	   write(92,*)
	   do i=1,M
           do j=1,M
	      write(92,4033) i,j,nuc(i), nuc(j),VAAAcuadrada(i,j)+VAABcuadrada(i,j)+VBACcuadrada(i,j)
	   end do
	   end do
	   close(92)
	   stop "check ECP routines, ECP FOCK MATRIX writed in ECP-FOCK"
	end if
	4031 format("i",i2,1x,"j",i2,1x,"V(i,j)",f18.15,1x,"V(j,i)",f18.15,1x)
	4032 format(7x,"i",5x,"j",3x,"nuci",2x,"nucj",7x,"ECP-FOCK")
!le agregue a 4033 para q escriba los nuc
	4033 format(5x,i3,3x,i3,3x,i3,3x,i3,4x,f18.15)
	end subroutine TEST_SYMMETRY

	subroutine TEST_COPY_FOCK
!testea que la matriz cuadrada y la matriz guardada como vector sean iguales
	use ECP_mod, only : VAAAcuadrada,VAABcuadrada,VBACcuadrada,VAAA,VAAB,VBAC
	use garcha_mod, only : nshell
	implicit none
	integer :: pos,i,j,M
	logical :: Er1,Er2,Er3
	Er1 = .false.
	Er2 = .false.
	Er3 = .false.
	M=nshell(0)+nshell(1)+nshell(2)
	do i=1,M
	do j=1,M
	   if (i .ge. j) then
	      pos=i+(1-j)*(j-2*M)/2 
	      if (VAAA(pos) .ne. VAAAcuadrada(i,j)) Er1 = .true.
	      if (VAAB(pos) .ne. VAABcuadrada(i,j)) Er2 = .true.
	      if (VBAC(pos) .ne. VBACcuadrada(i,j)) Er3 = .true.
	   end if
	end do
	end do

	if ( Er1 ) stop " Error in ECP-FOCK vector VAAA"
	if ( Er2 ) stop " Error in ECP-FOCK vector VAAB"
	if ( Er3 ) stop " Error in ECP-FOCK vector VBAC"
	end subroutine

