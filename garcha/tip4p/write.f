c-----------------------------------------------------------------
c Subroutine for storing basis functions and auxiliar basis
c set
c
c
c generates grid points
c evaluates both type of functions at every point and
c writes the results down on disk, so it can be used later
c In this case, the grid is generated randomly at this step
c and is not changed throughout the SCF procedure.
c 30-10-92, Dario Estrin
c-----------------------------------------------------------------
      subroutine write(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,
     >                NCOa,NCOb,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)
      implicit real*8 (a-h,o-z)
      logical NORM,dens1,SVD,integ,OPEN
      integer iconst,igrid,igrid2
      INCLUDE 'param'
      parameter (pi=3.14159265358979312D0,pi2=6.28318530717958623D0)
c input
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng),Iz(nt)
      dimension r(nt,3),nshell(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),nshelld(0:3)
      dimension ncontd(Md),xi(3),RR(ntq,ntq),P(ntq)
      dimension AFUNC(ngd),ds(ntq),RMM(*),F(ngd)
c
c
      common /Ll/ Ll(3)
      common /fit/ Nang,dens1,integ,Iexch,igrid,igrid2
      common /intg1/ e(50,3),wang(50),Nr(54)
      common /Sys/ SVD,iconst
      common /radii/ Rm(54)
      common /Nc/ Ndens
c
c open files ------------
c
      open(unit=8,file='scratch/fit1',form='unformatted')
      open(unit=9,file='scratch/fit2',form='unformatted')
c-------------------------
c
c
      Ndens=1
      NCO=NCOa
      M18=M17+Md*(Md+1)/2
c
      do 44 i=1,natom
      do 44 j=1,natom
       t0=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+(r(i,3)-r(j,3))**2
       RR(i,j)=sqrt(t0)
 44    continue
c
c
      do 43 l=1,3
 43    Ll(l)=l*(l-1)/2
c
      ns=nshell(0)
      npp=nshell(1)
      nd=nshell(2)
      M2=2*m
c
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
c
      sq3=1.D0
      if (NORM) then
      sq3=sqrt(3.D0)
      endif
c
         MMd=Md*(Md+1)/2
        do 2 k=1,MMd
 2       RMM(M17+k-1)=0.D0
c------------------------------------------------------------------------
c loop 12 , over all grid  -----------------------------
c loop 99, over all fitting functions
c everything goes here, the rest is the same
c xi vector of 3 dimensions
c yi density functional for xi ( all necessary data in common)
c
      DO 12 na=1,natom
c
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
c
c--------------------------------------------
       do 15  i=1,Nang
c
       xi(1)=r(na,1)+r1*e(i,1)
       xi(2)=r(na,2)+r1*e(i,2)
       xi(3)=r(na,3)+r1*e(i,3)
c
        do 21 k=1,natom
         ds(k)=(xi(1)-r(k,1))**2+(xi(2)-r(k,2))**2+(xi(3)-r(k,3))**2
 21     continue
c
         if (OPEN) then
         call DNSOP(DENSA,DENSB,F,Xi,ds,NORM,Nuc,ncont,nshell,
     >     a,c,r,M,M18,NCOa,NCOb,RMM)
         call potop(Iexch,DENSA,DENSB,yiex,yiec,y2a,y2b)
         DENS=DENSA+DENSB
         else
c
          call DNS(DENS,F,Xi,ds,NORM,Nuc,ncont,nshell,
     >     a,c,r,M,M18,NCO,RMM)
        call pot(Iexch,dens,yiex,yiec,y2i)
          endif
c
          yi = yiex + yiec
*
          call ABC(F,M,9)
        tmp0=wrad*wang(i)
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
        tmp=abs(dens/yi)*tmp0*PF
c
        if (dens.eq.0.0) then
         tmp=0.0D0
         yi=0.0D0
         y2i=0.D0
         y2a=0.0D0
         y2b=0.0D0
        endif
c-------------------------------
c calculation of all fitting functions at xi
c
        do 18 k=1,Md
 18      AFUNC(k)=0.D0
c s type
        do 19 k=1,nsd
c
         dd=ds(Nucd(k))
        do 19 nk=1,ncontd(k)
c
         t1=pi2/ad(k,nk)
         uu=ad(k,nk)*dd
         term=cd(k,nk)*t1*FUNCT(0,uu)
 19      AFUNC(k)=AFUNC(k)+term
c
c p type
        do 29 k=nsd+1,nsd+npd,3
c
         dd=ds(Nucd(k))
        do 29 nk=1,ncontd(k)
c       
         t1=pi2/ad(k,nk)
         uu=ad(k,nk)*dd
         t2=cd(k,nk)*t1*FUNCT(1,uu)
c
         do 29 l1=1,3
         term=(xi(l1)-r(Nucd(k),l1))*t2
         kk=k+l1-1
 29      AFUNC(kk)=AFUNC(kk)+term
c
c d type
        do 39 k=nsd+npd+1,Md,6
c
         dd=ds(Nucd(k))
        do 39 nk=1,ncontd(k)
c       
         t=pi2/ad(k,nk)*cd(k,nk)
         uu=ad(k,nk)*dd
         t0=FUNCT(0,uu)
         t1=FUNCT(1,uu)
         t2=FUNCT(2,uu)
         tt=(t0-t1)/(2.*ad(k,nk))
c
         do 39 l1=1,3
          t3=xi(l1)-r(Nucd(k),l1)
c
         do 39 l2=1,l1
          t4=xi(l2)-r(Nucd(k),l2)
          t5=t2*t3*t4
c
         fac=1.D0
          if(l1.eq.l2) then
           t5=t5+tt
           fac=sq3
          endif
c
          term=t5*t/fac
         l12=Ll(l1)+l2
         kk=k+l12-1
 39      AFUNC(kk)=AFUNC(kk)+term
c
c
c------------------------------------------------
c writes down the matrix for least squares procedure
c
#ifdef essl
c ESSL requires lower packed form
        kk=M17-1
        DO 11 j=1,Md
         tt=AFUNC(j)*tmp
c
         DO 11 i2=j,Md
         kk=kk+1
         RMM(kk)=RMM(kk)+AFUNC(i2)*tt
 11      CONTINUE
c
#endif
#ifdef pack
c LAPACK requires upper packed form
        kk=M17-1
        DO 111 j=1,Md
         tt=AFUNC(j)*tmp
c
         DO 111 i2=1,j
         kk=kk+1
         RMM(kk)=RMM(kk)+AFUNC(i2)*tt
 111     CONTINUE
#endif
c
         do 118 i2=1,Md
 118     AFUNC(i2)=AFUNC(i2)*tmp
c
c
c writes down results for auxiliar basis set
c
      call ABC(AFUNC,Md,8)
c
c-----------------------------------------------------
c
 151  continue
 15   continue
 16   continue
 12   continue
c
c ESSL OPTION
#ifdef essl
      CALL DPPF(RMM(M17),Md,1)
#endif
c
c LAPACK OPTION
#ifdef pack
      call dppco(RMM(M17),Md,rcond,aux,info)
#endif
c
c
      close(8)
      close(9)
      return
      END
c-------------------------------------------------------------
c subroutine that writes down the vector
c
       subroutine ABC(AFUNC,M,nunit)
       implicit real*8 (a-h,o-z)
        dimension AFUNC(M)
        write(nunit) AFUNC
        return
       end
c--------------------------------------------------------------
