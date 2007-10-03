c-----------------------------------------------------------------
c This subroutine deals with the exchange-correlation energy term
c Numerical integration , given by Becke's grid
c
c generates grid points
c calls to a function that at a certain space point gives density
c functional for the particular density functional chosen
c
c Output: Exchange-correlation energy
c 11-2-93
c-----------------------------------------------------------------
      subroutine exchnumop(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >                   M18,NCOa,NCOb,Exc,nopt)
      implicit real*8 (a-h,o-z)
      logical NORM,integ,dens1
      integer igrid,igrid2
      INCLUDE 'param'
      parameter (pi=3.14159265358979312D0,pi2=6.28318530717958623D0)
c
c input
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng),Iz(ntq)
      dimension r(nt,3),nshell(0:3)
      dimension xi(3),aux(ng)
      dimension RMM(*),ds(ntq)
c Rm: 1/2 Slater's radius, Nr : # shells
      dimension RR(ntq,ntq),P(ntq)
c
c
c
      common /Ll/ Ll(3)
      common /fit/ Nang,dens1,integ,Iexch,igrid,igrid2
      common /intg2/ e(116,3),wang(116),Nr(54),e3(194,3),wang3(194)
      common /radii/ Rm(54)
      common /propt/ idip,ipol,ispin
c
c-----------------
c
c
      Exc=0.0D0
      excha=0.0D0
      ecorr=0.0D0
      ss0=0.0D0
c
c ------------------------------------------------
c
      do 43 l=1,3
 43    Ll(l)=l*(l-1)/2
c
      do 44 i=1,natom
      do 44 j=1,natom
       t0=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+(r(i,3)-r(j,3))**2
       RR(i,j)=sqrt(t0)
 44    continue
c
      ns=nshell(0)
      npp=nshell(1)
      nd=nshell(2)
c
c
c-------------------------------------------------------------
c loop 12 , over all grid  -----------------------------
c loop 99, over all fitting functions
c everything goes here, the rest is the same
c xi vector of 3 dimensions
c yi density functional for xi ( all necessary data in common)
c
c # angular points in the grid
c
      if(igrid.eq.1) then
        npoint=116
      else
        npoint=194
      endif
c------------------------------
      DO 12 na=1,natom
c
c     write(*,*) na,Iz(na),Nr(Iz(na)),Rm(Iz(na))
       do 16 n=1,Nr(Iz(na))
c
       t0=pi/(Nr(Iz(na))+1)
       t1=t0*n
       x=cos(t1)
c w: weight
       w=t0*abs(sin(t1))
c multiply also by radial part
c
       r1=Rm(Iz(na))*(1.D0+x)/(1.D0-x)
c      write(*,*) n,x,r1
       wrad=w * r1**2
c
       wrad=wrad*Rm(Iz(na)) * 2.D0 /(1.D0-x)**2
c
c Angular part now
c Grid  given by V.I.Lebedev's paper
c  110 points or 194 points, according to igrid
c
      do 15 i=1,npoint
c
      if(igrid.eq.1) then
         xi(1)=r(na,1)+r1*e(i,1)
         xi(2)=r(na,2)+r1*e(i,2)
         xi(3)=r(na,3)+r1*e(i,3)
         tmp0=wrad*wang(i)
       else
         xi(1)=r(na,1)+r1*e3(i,1)
         xi(2)=r(na,2)+r1*e3(i,2)
         xi(3)=r(na,3)+r1*e3(i,3)
         tmp0=wrad*wang3(i)
       endif
c
c
c
c
        do 21 k=1,natom
         ds(k)=(xi(1)-r(k,1))**2+(xi(2)-r(k,2))**2+(xi(3)-r(k,3))**2
 21     continue
c
c
        if (Iexch.le.3) then
c local density functionals, no gradients needed
         call DNSOP(DENSA,DENSB,aux,Xi,ds,NORM,Nuc,ncont,nshell,
     >     a,c,r,M,M18,NCOa,NCOb,RMM)
         call potop(Iexch,DENSA,DENSB,yiex,yiec,y2a,y2b)
         DENS=DENSA+DENSB
        else
c non local density functionals, gradients and 2nd derivatives needed
         call DNSGOP(DENSA,DENSB,aux,aDx,bDx,aDy,bDy,aDz,bDz,aDxx,bDxx,
     >        aDyy,bDyy,aDzz,bDzz,aDxy,bDxy,aDyz,bDyz,aDxz,bDxz,
     >        Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,M,M18,NCOa,NCOb,RMM)
         DENS=DENSA+DENSB
         call potgop(Iexch,densa,densb,adx,bdx,ady,bdy,adz,bdz,adxx,
     >        bdxx,adyy,bdyy,adzz,bdzz,adxy,bdxy,adyz,bdyz,adxz,bdxz,
     >        yiex,yiec,y2a,y2b)
        endif
c
        yi = yiex + yiec
c
c NUMERICAL INTEGRATION PART  --------------------------------
c weight for the numerical integration is tmp0
c 
        PP=0.0D0
        do 119 nb=1,natom
         P(nb)=1.D0
c
         rnb=(xi(1)-r(nb,1))**2+ (xi(2)-r(nb,2))**2+
     >       (xi(3)-r(nb,3))**2
         rnb=sqrt(rnb)
c
         do 120 nc=1,natom
         if (nc.eq.nb) goto 121
c
         rnc=(xi(1)-r(nc,1))**2+ (xi(2)-r(nc,2))**2+
     >       (xi(3)-r(nc,3))**2
         rnc=sqrt(rnc)
c
         u=(rnb-rnc)/RR(nb,nc)
c
c heteronuclear correction
c
         x=Rm(Iz(nb))/Rm(Iz(nc))
         x1=(x-1.D0)/(x+1.D0)
         aij=x1/(x1**2-1.0D0)
         u=u+aij*(1.D0-u**2)
c
         p1=1.5D0*u-0.5D0*u**3
         p2=1.5D0*p1-0.5D0*p1**3
         p3=1.5D0*p2-0.5D0*p2**3
c        p4=1.5D0*p3-0.5D0*p3**3
c        p5=1.5D0*p4-0.5D0*p4**3
         s=0.5D0*(1.D0-p3)
         P(nb)=P(nb)*s
c
 121   continue
 120   continue
        PP=PP+P(nb)
 119    continue
c
c
        PF=P(na)/PP
c ss : integrated density, check
        t1=PF*dens*tmp0
        excha = excha + t1*yiex
        ecorr = ecorr + t1*yiec 
        ss0=ss0 + t1
c
c---------------------------------------------------------
c
 15    CONTINUE
c
 161  continue
 16   continue
 12   continue
c-------------------------------------------------------
*
      Exc=excha+ecorr
*
      if (nopt.eq.0) then
        write(*,610)
        write(*,620) excha,ecorr,ss0
      endif
*
c SPIN POLARIZATION CALCULATION
c
       if (ispin.eq.1) then
        write(*,*)
        write(*,800)
        write(*,810)
       do k=1,natom
c
        xi(1)=r(k,1)
        xi(2)=r(k,2)
        xi(3)=r(k,3)
c
        do j=1,natom
         ds(j)=(xi(1)-r(j,1))**2+(xi(2)-r(j,2))**2+(xi(3)-r(j,3))**2
        enddo
c
        call DNSOP(densa,densb,aux,xi,ds,NORM,Nuc,ncont,nshell,a,c,r,
     >             M,M18,NCOa,NCOb,RMM)
        write(*,820) k,Iz(k),densa+densb,densb-densa
c
        enddo
        endif
c
 800  format('  SPIN POLARIZATION ON NUCLEI')
 810  format('ATOM #',4x,'ATOM TYPE',4x,'TOT DENSITY',6x,'POLARIZATION')
 820  format (I3,9x,I3,6x,F13.4,5x,F10.4)
 610  format(2x,'EXCHANGE',13x,'CORRELATION',7x,'INTEGRATED DENSITY')
 620  format(F14.7,4x,F14.7,4x,F14.7)
      return
      END
c-------------------------------------------------------------
