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
      subroutine exchnum(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >                   M18,NCO,Exc,nopt)
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
      common /intg2/ e(116,3),wang(116),Nr(0:54),e3(194,3),wang3(194)
      common /radii/ Rm(0:54)
c
c
c
      excha=0.0D0
      ecorr=0.0D0
      Exc=0.0D0
      ss0=0.0D0
c
c -------------------------------------------------------------
c
      do 43 l=1,3
 43    Ll(l)=l*(l-1)/2
c
c     Loop entre todos los atomos: RR(i,j): dist entre atomo i y j
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
c       write(*,*) 'zz',na,n,Rm(Iz(na)),t0,t1,x
       wrad=w * r1**2
c       write(*,*) 'r1',r1
c
       wrad=wrad*Rm(Iz(na)) * 2.D0 /(1.D0-x)**2
c
c Angular part now
c Grid  given by V.I.Lebedev's paper
c  110 points or 194 points, according to igrid
c
      do 15 i=1,npoint
c
c      write(*,*) 'Grilla ',igrid
      if(igrid.eq.1) then
         xi(1)=r(na,1)+r1*e(i,1)
         xi(2)=r(na,2)+r1*e(i,2)
         xi(3)=r(na,3)+r1*e(i,3)
         tmp0=wrad*wang(i)
c         write(*,*) 'wang(',i-1,')=',wang(i),' e(i)=(',e(i,1),',',
c     >     e(i,2),',',e(i,3),')'
c         write(*,*) 'atomo: ',na,' punto: ',i,'Iz:',Iz(na),
c     > 'Nr:',Nr(Iz(na)),'tmp',tmp0,'pos:',xi(1),xi(2),xi(3)
       else
         xi(1)=r(na,1)+r1*e3(i,1)
         xi(2)=r(na,2)+r1*e3(i,2)
         xi(3)=r(na,3)+r1*e3(i,3)
         tmp0=wrad*wang3(i)
c         write(*,*) 'wang(',i-1,')=',wang3(i),' e(i)=(',e3(i,1),',',
c     >     e3(i,2),',',e3(i,3),')'         
c         write(*,*) 'atomo: ',na,' punto: ',i,'Iz:',Iz(na),
c     > 'Nr:',Nr(Iz(na)),'tmp',tmp0,'pos:',xi(1),xi(2),xi(3)
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
         call DNS(DENS,aux,Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,
     >            M,M18,NCO,RMM)
         dxi=DENS

c         do l=1,m
c           write(958,*) 'func',npoint*((na-1)*60+n-1)+i-1,aux(l)
c         enddo
         
         call pot(Iexch,dxi,yiex,yiec,y2i)
c         write(*,*) 'dens:',DENS,'ex',yiex,'cor',yiec
        else
c non local density functionals, gradients and 2nd derivatives needed
         call DNSG(DENS,aux,Dx,Dy,Dz,Dxx,Dyy,Dzz,Dxy,Dyz,Dxz,
     >           Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,M,M18,NCO,RMM)
         dxi=DENS
         call potg(Iexch,dxi,dx,dy,dz,dxx,dyy,dzz,dxy,dyz,dxz,yiex,yiec,
     >             y2i)
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
c        write(*,*) 't1',t1
c        if (igrid.eq.1) then
c        write(*,*) 'wang',wang(i),'wrad',wrad
c        else
c          write(*,*) 'wang3',wang3(i),'wrad',wrad
c        endif
        excha = excha + t1*yiex
        ecorr = ecorr + t1*yiec
        PP=t1*yiex + t1*yiec
c        write(*,*) 'zz',npoint*((na-1)*60+n-1)+i-1,dens,PP
        ss0=ss0 + t1
        Npt=Npt+1
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
      if (nopt.eq.0) then
        write(*,610)
        write(*,620) excha,ecorr,ss0
      endif
*     
      return
 610  format(2x,'EXCHANGE',13x,'CORRELATION',7x,'INTEGRATED DENSITY')
 620  format(F14.7,4x,F14.7,4x,F14.7)
 777  format (3(F14.6,2x))
      END
c-------------------------------------------------------------
