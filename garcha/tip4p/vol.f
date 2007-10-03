c-------------------------------------------------------------------
c
c A cube of points containing the molecule using Van der Waals distances
c for defining the max size is generated.
c Points inside a fraction of the Van der Waals radius are
c automatically included in the volume, others need a density
c calculation.
c The thershold used is 10-3 a.u
c The volume is estimated by numerically integrating the space
c inside this isodensity surface. A 'crude' numerical integration
c is performed, since results are only an estimation of the cavity
c size.
c If a value is given in the input, this calculation is not
c performed.
c 
c After the volume evaluation, the value is scaled (1.33) and
c a radius is evaluated, 0.50 angstroms are finally added to the
c obtained value , in order to consider the first solvation layer.
c
c Dario Estrin
c Buenos Aires, ARGENTINA 18-01-94
c
c 
c loop over all basis functions
c-------------------------------------------------------------------
      subroutine vol(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,M18,NCOa,NCOb,OPEN,RMM,a0)
c
      implicit real*8 (a-h,o-z)
      INCLUDE 'param'
      logical NORM,OUT,OPEN
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(natom)
      dimension r(nt,3),nshell(0:3)
      dimension RMM(*),xi(3),aux(ng),ds(ntq),vdw(54)
c
      COMMON /TABLE/ STR(880,0:21)
c
      MM=M*(M+1)/2
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, F also uses the same position after S was used
      M5=M3+MM
c
c data in angstroms
      data vdw /1.50,1.50,1.50,2.00,2.00,2.00,2.00,2.00,2.00,2.00,
     > 2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,
     > 2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,
     > 2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,
     > 2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50/
c
c
c changes to a.u and also scales slightly 
      do i=1,54
       vdw(i)=vdw(i)/0.529177*1.20D0
      enddo
c
      step=0.25D0
c
c finds minimum and maximum for X,Y,Z in order to define the box
       Xmin=r(1,1)
       Xmax=r(1,1)
       Ymin=r(1,2)
       Ymax=r(1,2)
       Zmin=r(1,3)
       Zmax=r(1,3)

       do n=1,natom
        if (r(n,1).ge.Xmax) then
         Xmax=r(n,1)
         nxmax=n
        endif
        if (r(n,1).le.Xmin) then
         Xmin=r(n,1)
         nxmin=n
        endif
c
        if (r(n,2).ge.Ymax) then
         Ymax=r(n,2)
         nymax=n
        endif
        if (r(n,2).le.Ymin) then
         Ymin=r(n,2)
         nymin=n
        endif
c
        if (r(n,3).ge.Zmax) then
         Zmax=r(n,3)
         nzmax=n
        endif
        if (r(n,3).le.Zmin) then
         Zmin=r(n,3)
         nzmin=n
        endif
c 
        enddo

c
       Xmin=Xmin-vdw(Iz(nxmin))
       Xmax=Xmax+vdw(Iz(nxmax))
       Ymin=Ymin-vdw(Iz(nymin))
       Ymax=Ymax+vdw(Iz(nymax))
       Zmin=Zmin-vdw(Iz(nzmin))
       Zmax=Zmax+vdw(Iz(nzmax))
c
       tv=step**3
c GRID GENERATION ------------------------------------------------
c
       Npoint=0
       do 10  x=Xmin,Xmax,step
       do 10  y=Ymin,Ymax,step
       do 10  z=Zmin,Zmax,step
c
       xi(1)=x
       xi(2)=y
       xi(3)=z
c
         iout=0
         OUT=.false.
        do 21 k=1,natom
         ds(k)=(xi(1)-r(k,1))**2+(xi(2)-r(k,2))**2+(xi(3)-r(k,3))**2
         
         if (sqrt(ds(k)).gt.vdw(Iz(k))) then
          iout=iout+1
         endif
c
         if (sqrt(ds(k)).lt.vdw(Iz(k))*0.60D0) then
          OUT=.true.
         endif 
 21     continue

c
        dens=0.0D0
        if ((iout.ne.natom).and.(.not.OUT)) then
c
        if (OPEN) then
c------------------OPEN SHELL CASE
        if (Iexch.le.3) then
c local density functionals, no gradients needed
         call DNSOP(DENSA,DENSB,aux,Xi,ds,NORM,Nuc,ncont,nshell,
     >     a,c,r,M,M18,NCOa,NCOb,RMM)
         DENS=DENSA+DENSB
        else
c non local density functionals, gradients and 2nd derivatives needed
         call DNSGOP(DENSA,DENSB,aux,aDx,bDx,aDy,bDy,aDz,bDz,aDxx,bDxx,
     >        aDyy,bDyy,aDzz,bDzz,aDxy,bDxy,aDyz,bDyz,aDxz,bDxz,
     >        Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,M,M18,NCOa,NCOb,RMM)
         DENS=DENSA+DENSB
        endif
c -----------------CLOSED SHELL CASE
        else
        if (Iexch.le.3) then
c local density functionals, no gradients needed
         call DNS(DENS,aux,Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,
     >            M,M18,NCOa,RMM)
        else
c non local density functionals, gradients and 2nd derivatives needed
         call DNSG(DENS,aux,Dx,Dy,Dz,Dxx,Dyy,Dzz,Dxy,Dyz,Dxz,
     >           Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,M,M18,NCOa,RMM)
        endif
c
       Npoint=Npoint+1
        endif
c------------------------------------------------------------
        endif
c--------------------------------------------------------------
c
        if (dens.gt.1.0D-03) then
c        Vlm=Vlm + tv
         Nvol=Nvol+1
        endif
c
        if (OUT) then
c        Vlm=Vlm + tv
         Nin=Nin+1
        endif
 10    continue
c
c volume calculation, simply number of points time size of delta
       Vlm=(Nin+Nvol)*tv
c
c this is a0 in a.u
      a0=(0.75D0/pi*Vlm*1.33)**0.33333333333333D0+0.50D0/0.529177
c
c prints out a0 in angstroms
c     write(*,*) 'VOLUME, RADIUS',Vlm,a0*0.529177
c     write(*,*) Npoint,Nvol,Nin
c
c     do i=1,54
c      vdw(i)=vdw(i)*0.529177/1.20D0
c     enddo

      end
