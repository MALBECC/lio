      subroutine exch(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,Md,
     >  ncontd,nshelld,Nucd,cd,ad,RMM,NCO,NCOb,M17,B1)
c
      implicit real*8 (a-h,o-z)
      logical NORM,SVD,dens1,integ,OPEN
      INCLUDE 'param'
      parameter (pi=3.14159265358979312D0,pi2=6.28318530717958623D0)
c
c input
      integer iconst,igrid,igrid2
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng),Iz(ntq)
      dimension r(nt,3),nshell(0:3),nshelld(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension xi(3),aux(ng)
      dimension AFUNC(ngd),RMM(*),ds(ntq)
c
c output
c fit for exchange correlation density and potential if different
      dimension B1(ngd,3)
c B1(i,1) exchange correlation density
c B1(i,2) exchange correlation potential
c
c
      common /Ll/ Ll(3)
      common /fit/ Nang,dens1,integ,Iexch,igrid,igrid2
      common /intg1/ e(50,3),wang(50),Nr(0:54)
      common /Sys/ SVD,iconst
      common /radii/ Rm(0:54)
      common /Nc/ Ndens
c---------------------
*
      NCOa = NCO
*
c
      index=0
c
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
c
      MMd=Md*(Md+1)/2
      M18=M17+MMd
c vectors of MO beta
      M18b=M18+M*NCO
c weights (open shell case)
      if (OPEN) then
      M19=M18b+M*NCOb
      else
      M19=M18b
      endif
c
c--------------------
c
c
c ------------------------------------------------
c
      do 43 l=1,3
   43    Ll(l)=l*(l-1)/2
c
c
c initialization
        do 1 k=1,Md
         B1(k,1)=0.D0
         B1(k,2)=0.D0
    1       B1(k,3)=0.0D0
c
c
c
      sq3=1.D0
      if (NORM) then
      sq3=sqrt(3.D0)
      endif
c
c-------------------------------------------------------------
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
c multiply also by radial part
c
       r1=Rm(Iz(na))*(1.D0+x)/(1.D0-x)
c
      do 15 i=1,Nang
c
c
c
       xi(1)=r(na,1)+r1*e(i,1)
       xi(2)=r(na,2)+r1*e(i,2)
       xi(3)=r(na,3)+r1*e(i,3)
c
        do 21 k=1,natom
         ds(k)=(xi(1)-r(k,1))**2+(xi(2)-r(k,2))**2+(xi(3)-r(k,3))**2
   21     continue
c
c
        if (Iexch.le.3) then
c local density functionals, no gradients needed
         if (OPEN) then
c
         call DNSOP(DENSA,DENSB,aux,Xi,ds,NORM,Nuc,ncont,nshell,
     >     a,c,r,M,M18,NCO,NCOb,RMM)
         call potop(Iexch,DENSA,DENSB,yiex,yiec,y2a,y2b)
         dxi=DENSA+DENSB
         DENS=dxi
c
         else
         call DNS(DENS,aux,Xi,ds,NORM,Nuc,ncont,nshell,
     >     a,c,r,M,M18,NCO,RMM)
         dxi=DENS
         call pot(Iexch,dxi,yiex,yiec,y2i)
         endif
c
        else
c non local density functionals, gradients and 2nd derivatives needed
         if (OPEN) then
         call DNSGOP(DENSA,DENSB,aux,aDx,bDx,aDy,bDy,aDz,bDz,aDxx,bDxx,
     >        aDyy,bDyy,aDzz,bDzz,aDxy,bDxy,aDyz,bDyz,aDxz,bDxz,
     >        Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,M,M18,NCOa,NCOb,RMM)
         DENS=DENSA+DENSB
         dxia = DENSA
         dxib = DENSB
         dxi = DENS
         call potgop(Iexch,dxia,dxib,adx,bdx,ady,bdy,adz,bdz,adxx,
     >        bdxx,adyy,bdyy,adzz,bdzz,adxy,bdxy,adyz,bdyz,adxz,bdxz,
     >        yiex,yiec,y2a,y2b)

         else
         call DNSG(DENS,aux,Dx,Dy,Dz,Dxx,Dyy,Dzz,Dxy,Dyz,Dxz,
     >           Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,M,M18,NCO,RMM)
         dxi=DENS
         call potg(Iexch,dxi,dx,dy,dz,dxx,dyy,dzz,dxy,dyz,dxz,yiex,
     >             yiec,y2i)
         endif
        endif
c
        yi = yiex + yiec
c---------------------------------------------------------
c
        if (dxi.eq.0.0D0) then
         tmp=0.0D0
         yi=0.0D0
         y2i=0.0D0
         y2a=0.0D0
         y2b=0.0D0
        endif
c
c calculates fitting functions  at xi ------
c
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
c -------
         t=pi2/ad(k,nk)*cd(k,nk)
         uu=ad(k,nk)*dd
         t0=FUNCT(0,uu)
         t1=FUNCT(1,uu)
         t2=FUNCT(2,uu)
         tt=(t0-t1)/(2.D0*ad(k,nk))
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
c--------
         l12=Ll(l1)+l2
         kk=k+l12-1
   39      AFUNC(kk)=AFUNC(kk)+term
c
c-------------------------------------------------------
c
c -- THIS PART IS THE NORMAL EQUATION METHOD -----
c
c
        if (OPEN) then
        DO 119 j=1,Md
         tt=AFUNC(j)*RMM(M19+index)
         B1(j,1)=B1(j,1)+yi*tt
         B1(j,2)=B1(j,2)+y2a*tt
         B1(j,3)=B1(j,3)+y2b*tt
  119     CONTINUE
c
          else
c
        DO 11 j=1,Md
         tt=AFUNC(j)*RMM(M19+index)
         B1(j,1)=B1(j,1)+yi*tt
         B1(j,2)=B1(j,2)+y2i*tt
   11      CONTINUE
c
        endif
         index=index+1
c
   15    CONTINUE
c
   16   continue
   12   continue
c-------------------------------------------------------
c ESSL OPTION
#ifdef essl
      CALL DPPS(RMM(M17),Md,B1(1,1),1)
      CALL DPPS(RMM(M17),Md,B1(1,2),1)
c
      if (OPEN) then
       CALL DPPS(RMM(M17),Md,B1(1,3),1)
      endif
c
#endif
c
c LINPACK OPTION
c
#ifdef pack
      call dppsl(RMM(M17),Md,B1(1,1))
      call dppsl(RMM(M17),Md,B1(1,2))
c
      if (OPEN) then
      call dppsl(RMM(M17),Md,B1(1,3))
      endif
c
#endif
c --------------------------------------------

      return
  500 format(2F9.5)
      END
c-------------------------------------------------------------
