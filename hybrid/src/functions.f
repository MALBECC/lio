c funciones q utiliza el programa

        function dist(x1,y1,z1,x2,y2,z2)
        implicit none
        double precision dist,x1,y1,z1,x2,y2,z2
        dist = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        dist = sqrt(dist)
        return
        end function dist

        function dist2(x1,y1,z1,x2,y2,z2)
        implicit none
        double precision dist2,x1,y1,z1,x2,y2,z2
        dist2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        return
        end function dist2

        function angle(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit none
        double precision scalar,dist,angle,pi
        double precision x1,y1,z1,x2,y2,z2,x3,y3,z3
	pi=DACOS(-1.d0)
c       calculo del producto escalar
        scalar = (x1-x2)*(x3-x2)+(y1-y2)*(y3-y2)+(z1-z2)*(z3-z2)
        angle = dist(x1,y1,z1,x2,y2,z2)*dist(x3,y3,z3,x2,y2,z2)
        angle = scalar/angle
	if (angle.ge.1.d0) then
	angle = 1.d0
	elseif(angle.le.-1.d0) then
        angle=-1.d0
        endif
        angle = dACOS(angle)*180d0/pi
	return
        end function angle

        function scalar(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit none
        double precision scalar,pi
        double precision x1,y1,z1,x2,y2,z2,x3,y3,z3
        pi=DACOS(-1.d0)
c       calculo del producto escalar
        scalar = (x1-x2)*(x3-x2)+(y1-y2)*(y3-y2)+(z1-z2)*(z3-z2)
        end function scalar

        function dihedro(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
        implicit none
        double precision dihedro,dist,angle,l1,l2,pi,arg
        double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
        double precision l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
        pi=DACOS(-1.d0)
	mx = (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2)
        my = -((x1-x2)*(z3-z2) -(z1-z2)*(x3-x2))
        mz = (x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)
        nx = (y2-y3)*(z4-z3) - (z2-z3)*(y4-y3)
        ny = -((x2-x3)*(z4-z3) -(z2-z3)*(x4-x3))
        nz = (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)
c       calculo del prod scalar n*m
        scalar = mx*nx + my*ny + mz*nz
        m = mx**2d0 + my**2d0 + mz**2d0
        m = m**(0.5d0)
        n = nx**2d0 + ny**2d0 + nz**2d0
        n = n**(0.5d0)
c si el argumento da mal (m*n)=0.0 avisa y sale el dihe vale 500.0
	if(m*n.eq.0) then
	dihedro=500.d0
	STOP "problemas en el calculo de un diedro, Nick"
	go to 10
	endif
        arg=scalar/(m*n)
        if(arg.ge.1.0d0) then
	dihedro= 0.d0
	elseif (arg.le.-1.d0) then
	dihedro=180.d0
	else
        dihedro = dACOS(arg)
	endif
        dihedro = dihedro*180d0/pi
 10	continue
        end function dihedro

        function dihedro2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
        implicit none
        double precision dihedro2,dist,angle,l1,l2,pi
        double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
        double precision l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
 
c       vectors n and m
c       m = vectorial product 21 and 23
c       n = vectorial product 23 and 34
 
        pi=DACOS(-1.d0)
	mx = (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2)
        my = -((x1-x2)*(z3-z2) -(z1-z2)*(x3-x2))
        mz = (x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)
 
        nx = (y2-y3)*(z4-z3) - (z2-z3)*(y4-y3)
        ny = -((x2-x3)*(z4-z3) -(z2-z3)*(x4-x3))
        nz = (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)
 
c       scalar product n*m
        scalar = mx*nx + my*ny + mz*nz
        m = mx**2d0 + my**2d0 + mz**2d0
        m = m**(0.5d0)
        n = nx**2d0 + ny**2d0 + nz**2d0
        n = n**(0.5d0)
        dihedro2 = ACOS(scalar/(m*n))
        dihedro2 = dihedro2*180d0/pi
 
c       which cuadrant? 
c       plane generation with points 1 2 3
        A = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
        B = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
        C = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
        D = -x1*A-y1*B-z1*C
 
c	distance(l) at4 to ABCD=0 plane
        l1 = A*x4+B*y4+C*z4+D
        l2 = (A**2d0+B**2d0+C**2d0)
        l2 = l2**0.5d0
        l= l1/l2
 
c	if l>0 -> dihe<0 , else >0
        if(l.lt.0) then
                dihedro2 = 360d0-dihedro2
        endif
 
        end function dihedro2

c************************************************************************************
c subroutine that calculates the force of a dihedral angle in amber

	subroutine diheforce(nac,ramber,i1,i2,i3,i4,
     .                       atom,kd,eq,per,mult,fce)
c variables asoc a funcion como subrutina
        integer i1,i2,i3,i4
        double precision ramber(3,nac),fce(12),dtot
c i1-i4 numero de atomos 1 a 4  atom es que derivada tiene que calcular 
c atom=1 indica el primer atomo del dihe atomo=2 el segundo etc...
c fce(12) es la fuerza (derivada) para at1x,at1y,at1z,at2x....at4z
        integer  mult
        double precision  kd,eq,per
c variables asoc al calculo de la derivada
        double precision dih,dihedro,dist,rm,rn,mx,my,mz,nx,ny,nz,scal
        double precision dscalar,dm,dn,dmx,dmy,dmz,dnx,dny,dnz,dmn,
     .                  fx,fy,fz,E(3),fpx,fpxmasd,fpxmenosd,fr,step,prue
c variables generales
        integer i,j,k,l,m,n,nac,vez,coord,atom
        character exp
        double precision pi
	pi=DACOS(-1.d0)

c	write(*,*) "i1 i2 i3 i4", i1, i2, i3, i4

       call dihevars(   ramber(1,i1),ramber(2,i1),ramber(3,i1),
     .                  ramber(1,i2),ramber(2,i2),ramber(3,i2),
     .                  ramber(1,i3),ramber(2,i3),ramber(3,i3),
     .                  ramber(1,i4),ramber(2,i4),ramber(3,i4),
     .                  mx,my,mz,rm,nx,ny,nz,rn,dih)
        fce=0.0       
	scal=mx*nx + my*ny + mz*nz    
        dtot = (kd/mult)*(-dSIN((pi/180d0)*
     .  (per*dih-eq)))*(per)
        prue=scal/(rn*rm)
        prue=(1.d0-(prue)**2d0) 
        if (prue.lt.1.0d-15.and.prue.gt.-1.0d-15) go to 10  
        prue=dsqrt(prue) 
	dtot = -dtot/prue
	do j=1,3
	i=(atom-1)*3+j
	dmx=0.d0
        dmy=0.d0
        dmz=0.d0
        dnx=0.d0
        dny=0.d0
        dnz=0.d0
	if(i.eq.1) then
        dmy=ramber(3,i2)-ramber(3,i3)
        dmz=ramber(2,i3)-ramber(2,i2)
        elseif(i.eq.4) then
        dmy=ramber(3,i3)-ramber(3,i1)
        dmz=ramber(2,i1)-ramber(2,i3)
        dny=ramber(3,i3)-ramber(3,i4)
        dnz=ramber(2,i4)-ramber(2,i3)
        elseif(i.eq.7) then
        dmy=ramber(3,i1)-ramber(3,i2)
        dmz=ramber(2,i2)-ramber(2,i1)
        dny=ramber(3,i4)-ramber(3,i2)
        dnz=ramber(2,i2)-ramber(2,i4)
        elseif(i.eq.10) then
        dny=ramber(3,i2)-ramber(3,i3)
        dnz=ramber(2,i3)-ramber(2,i2)
        elseif(i.eq.2) then
        dmx=ramber(3,i3)-ramber(3,i2)
        dmz=ramber(1,i2)-ramber(1,i3)
        elseif(i.eq.5) then
        dmx=ramber(3,i1)-ramber(3,i3)
        dmz=ramber(1,i3)-ramber(1,i1)
        dnx=ramber(3,i4)-ramber(3,i3)
        dnz=ramber(1,i3)-ramber(1,i4)
        elseif(i.eq.8) then
        dmx=ramber(3,i2)-ramber(3,i1)
        dmz=ramber(1,i1)-ramber(1,i2)
        dnx=ramber(3,i2)-ramber(3,i4)
        dnz=ramber(1,i4)-ramber(1,i2)
        elseif(i.eq.11) then
        dnx=ramber(3,i3)-ramber(3,i2)
        dnz=ramber(1,i2)-ramber(1,i3)
        elseif(i.eq.3) then
        dmx=ramber(2,i2)-ramber(2,i3)
        dmy=ramber(1,i3)-ramber(1,i2)
        elseif(i.eq.6) then
        dmx=ramber(2,i3)-ramber(2,i1)
        dmy=ramber(1,i1)-ramber(1,i3)
        dnx=ramber(2,i3)-ramber(2,i4)
        dny=ramber(1,i4)-ramber(1,i3)
        elseif(i.eq.9) then
        dmx=ramber(2,i1)-ramber(2,i2)
        dmy=ramber(1,i2)-ramber(1,i1)
        dnx=ramber(2,i4)-ramber(2,i2)
        dny=ramber(1,i2)-ramber(1,i4)
        elseif(i.eq.12) then
        dnx=ramber(2,i2)-ramber(2,i3)
        dny=ramber(1,i3)-ramber(1,i2)
        endif
	dm=(mx*dmx+my*dmy+mz*dmz)/rm
        dn=(nx*dnx+ny*dny+nz*dnz)/rn
        dmn=rm*dn+rn*dm
        dscalar=nx*dmx+mx*dnx+ny*dmy+my*dny+nz*dmz+mz*dnz
	fce(i)= dtot*(dscalar*rm*rn-dmn*scal)/(rn*rm)**2
        enddo
  10	end
c************************************************************************* 
c subroutine that calculates the force of a dihedral angle in general

        subroutine diheforce2(nac,ramber,i1,i2,i3,i4,
     .                       atom,kd,fce)
 
c variables asoc a funcion como subrutina
        integer i1,i2,i3,i4
        double precision ramber(3,nac),fce(12),dtot
c i1-i4 numero de atomos 1 a 4  atom es que derivada tiene que calcular
c atom=1 indica el primer atomo del dihe atomo=2 el segundo etc...
c fce(12) es la fuerza (derivada) para at1x,at1y,at1z,at2x....at4z
        integer  mult
        double precision  kd,eq,per
c variables asoc al calculo de la derivada
        double precision dih,dihedro,dist,rm,rn,mx,my,mz,nx,ny,nz,scal
        double precision dscalar,dm,dn,dmx,dmy,dmz,dnx,dny,dnz,dmn,
     .                  fx,fy,fz,E(3),fpx,fpxmasd,fpxmenosd,fr,step,prue
c variables generales
        integer i,j,k,l,m,n,nac,vez,coord,atom
        character exp
        double precision pi
        pi=DACOS(-1.d0)
       call dihevars(   ramber(1,i1),ramber(2,i1),ramber(3,i1),
     .                  ramber(1,i2),ramber(2,i2),ramber(3,i2),
     .                  ramber(1,i3),ramber(2,i3),ramber(3,i3),
     .                  ramber(1,i4),ramber(2,i4),ramber(3,i4),
     .                  mx,my,mz,rm,nx,ny,nz,rn,dih)
        fce=0.d0
        scal=mx*nx + my*ny + mz*nz
	dtot=kd
        prue=scal/(rn*rm)
        prue=(1.d0-(prue)**2d0)
        if (prue.lt.1.0d-15.and.prue.gt.-1.0d-15) go to 10
        prue=dsqrt(prue)
        dtot = -dtot/prue
        do j=1,3
        i=(atom-1)*3+j
        dmx=0.d0
        dmy=0.d0
        dmz=0.d0
        dnx=0.d0
        dny=0.d0
        dnz=0.d0
        if(i.eq.1) then
        dmy=ramber(3,i2)-ramber(3,i3)
        dmz=ramber(2,i3)-ramber(2,i2)
        elseif(i.eq.4) then
        dmy=ramber(3,i3)-ramber(3,i1)
        dmz=ramber(2,i1)-ramber(2,i3)
        dny=ramber(3,i3)-ramber(3,i4)
        dnz=ramber(2,i4)-ramber(2,i3)
        elseif(i.eq.7) then
        dmy=ramber(3,i1)-ramber(3,i2)
        dmz=ramber(2,i2)-ramber(2,i1)
        dny=ramber(3,i4)-ramber(3,i2)
        dnz=ramber(2,i2)-ramber(2,i4)
        elseif(i.eq.10) then
        dny=ramber(3,i2)-ramber(3,i3)
        dnz=ramber(2,i3)-ramber(2,i2)
        elseif(i.eq.2) then
        dmx=ramber(3,i3)-ramber(3,i2)
        dmz=ramber(1,i2)-ramber(1,i3)
        elseif(i.eq.5) then
        dmx=ramber(3,i1)-ramber(3,i3)
        dmz=ramber(1,i3)-ramber(1,i1)
        dnx=ramber(3,i4)-ramber(3,i3)
        dnz=ramber(1,i3)-ramber(1,i4)
        elseif(i.eq.8) then
        dmx=ramber(3,i2)-ramber(3,i1)
        dmz=ramber(1,i1)-ramber(1,i2)
        dnx=ramber(3,i2)-ramber(3,i4)
        dnz=ramber(1,i4)-ramber(1,i2)
        elseif(i.eq.11) then
        dnx=ramber(3,i3)-ramber(3,i2)
        dnz=ramber(1,i2)-ramber(1,i3)
        elseif(i.eq.3) then
        dmx=ramber(2,i2)-ramber(2,i3)
        dmy=ramber(1,i3)-ramber(1,i2)
        elseif(i.eq.6) then
        dmx=ramber(2,i3)-ramber(2,i1)
        dmy=ramber(1,i1)-ramber(1,i3)
        dnx=ramber(2,i3)-ramber(2,i4)
        dny=ramber(1,i4)-ramber(1,i3)
        elseif(i.eq.9) then
        dmx=ramber(2,i1)-ramber(2,i2)
        dmy=ramber(1,i2)-ramber(1,i1)
        dnx=ramber(2,i4)-ramber(2,i2)
        dny=ramber(1,i2)-ramber(1,i4)
        elseif(i.eq.12) then
        dnx=ramber(2,i2)-ramber(2,i3)
        dny=ramber(1,i3)-ramber(1,i2)
        endif
	dm=(mx*dmx+my*dmy+mz*dmz)/rm
        dn=(nx*dnx+ny*dny+nz*dnz)/rn
        dmn=rm*dn+rn*dm
        dscalar=nx*dmx+mx*dnx+ny*dmy+my*dny+nz*dmz+mz*dnz
        fce(i)= dtot*(dscalar*rm*rn-dmn*scal)/(rn*rm)**2d0
        enddo
  10    end
c*******************************************************************
c subroutine that calculates the variables of a dihedral angle

        subroutine dihevars(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     .                      mx,my,mz,m,nx,ny,nz,n,dihedro)
        implicit none
        double precision dihedro,dist,angle,l1,l2,pi,arg
        double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
        double precision l,A,B,C,D,m,n,nx,ny,nz,mx,my,mz,scalar
 	pi=DACOS(-1.d0)
        mx = (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2)
        my = -((x1-x2)*(z3-z2) -(z1-z2)*(x3-x2))
        mz = (x1-x2)*(y3-y2)-(y1-y2)*(x3-x2)
        nx = (y2-y3)*(z4-z3) - (z2-z3)*(y4-y3)
        ny = -((x2-x3)*(z4-z3) -(z2-z3)*(x4-x3))
        nz = (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)
c       calculo del prod scalar n*m
        scalar = mx*nx + my*ny + mz*nz
        m = mx**2d0 + my**2d0 + mz**2d0
        m = m**(0.5d0)
        n = nx**2d0 + ny**2d0 + nz**2d0
        n = n**(0.5d0)
        arg=scalar/(m*n)
        if(arg.ge.1.d0) then 
	dihedro=0.d0 
	elseif(arg.le.-1.d0) then
 	dihedro=180.d0
	else
        dihedro = dACOS(arg)
	endif
        dihedro = dihedro*180d0/pi 
 10	continue 
        end

c *************************************************************************
c subroutine that calculates position of 4 from 1-3 

        subroutine pos4(x1,y1,z1,x2,y2,z2,x3,y3,z3,r4,a4,d4,x4,y4,z4)
        implicit none
        double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
        double precision r4,a4,d4,xejx,yejx,zejx,xejz,yejz,zejz
        double precision pi,rejx,r23,rejz,x,y,z
        double precision l1,m1,n1,l2,m2,n2,l3,m3,n3
 
        pi=DACOS(-1.d0)
 
c       directrices cosines l1,m1,n1,etc
c       l1: x' axis respect to x 
c       m1 x'axis respect to y
c       n1 x' axis respect to z
c       2 and 3 respect to y and z axis
c       x' axis ejx
 
        xejx = (y3-y2)*(z1-z2) - (z3-z2)*(y1-y2)
        yejx = -((x3-x2)*(z1-z2) -(z3-z2)*(x1-x2))
        zejx = (x3-x2)*(y1-y2)-(y3-y2)*(x1-x2)
 
        rejx = dsqrt(xejx**2d0+yejx**2d0+zejx**2d0)
 
        l1 = xejx/rejx
 
        r23=dsqrt((x3-x2)**2d0+(y3-y2)**2d0+(z3-z2)**2d0)
        m1 = yejx/rejx
        n1= zejx/rejx
 
        l2 = (x3 - x2)/r23
        m2 = (y3 - y2)/r23
        n2 = (z3 - z2)/r23
 
 
      xejz = (yejx)*(z3-z2) - (zejx)*(y3-y2)
      yejz = -((xejx)*(z3-z2) -(zejx)*(x3-x2))
      zejz = (xejx)*(y3-y2)-(yejx)*(x3-x2)
 
        rejz = dsqrt(xejz**2d0+yejz**2d0+zejz**2d0)
 
        l3 = xejz/rejz
        m3 = yejz/rejz
        n3 = zejz/rejz

c	at4 coords in x' y' z' system
c	dihedral equivalent to angle in z' x' plane
c	180-angle equivalent to angle with y' axis
c	r4=distance at3 - at4 
        d4 = d4*pi/180.d0
        a4 =180.d0-a4
        a4 = a4*pi/180.d0
 
        z = r4*dSIN(a4)*dCOS(d4)
        x = r4*dSIN(a4)*dSIN(d4)
        y = r4*dCOS(a4)
 
        y = y + r23
 
c       translating 
 
        x4 = (l1*x + l2*y + l3*z + x2)
        y4 = (m1*x + m2*y + m3*z + y2)
        z4 = (n1*x + n2*y + n3*z + z2)
 
        return
        end subroutine pos4
 
c******************************************************************************
