c-------------------------------------------------------------------
c Integrals subroutine -Second part
c Gradients only
c 2 e integrals, 2 index : density fitting functions
c All of them are calculated
c using the Obara-Saika recursive method.
c
c
c loop over all basis functions
c now the basis is supposed to be ordered according to the type,
c all s, then all p, then all d, .....
c inside each type, are ordered in shells
c px,py,pz , dx2,dxy,dyy,dzx,dzy,dzz, .....
c
c ns ... marker for end of s
c np ... marker for end of p
c nd ... marker for end of d
c
c r(Nuc(i),j) j component of position of nucleus i , j=1,3
c Input :  density basis 
c Output: forces on nuclei 
c G matrix should be inverted, 
c later on, for evaluating  Coulomb terms
c-----------------------------------------------------------------
       module subm_int2G; contains
       subroutine int2G(f)

       use liotemp   , only: FUNCT
       use garcha_mod, only: RMM, ll, natom, M, Md, NORM, pi5, r
     >                     , af, nshelld, ad, nucd, ncontd, d, cd
!
       implicit none
! aux . things
       real*8  :: Q(3), f(natom,3)
       real*8  :: ti, tj, t0, t1, t2, t10, t11, t12, t12b
       real*8  :: t13, t13a, t13b, t14, t14a, t14b, t15, t15a, t15b
       real*8  :: t16, t16b, t17, t17b, t18, t18b, t20, t21, t22
       real*8  :: t23, t24, t25, t26, t27, t27a, t30, t31, t32, t33
       real*8  :: t40, t41, t42, t43

       real*8  :: z2, zij, z2a, zc, zc2
       real*8  :: f1, f2, fs, fp, fd
       real*8  :: cc, cc1, cc2, cci, ccj
       real*8  :: ccoef, factor, roz, roz2, sq3, alf, u

       real*8  :: sp, spj, sp1j, s1pk, s2pk
       real*8  :: s0s, s1s, s2s, s3s, s4s, s5s
       real*8  :: ps, pp, pd, pis, pip, pid, pjs, pjp, pjd
       real*8  :: pi0s, pi0p, pi0d, pj0s, pj0p, pj0d
       real*8  :: p1p, pi2s, pi2p, pi3s, pi4s, pj2s, pj2p, pj3s
       real*8  :: ds, dp, dpl, dsd
       real*8  :: d0p, d0pl, d1d, d1p, d1s, d2s, d3s, dd, df

       integer :: MM, MMd, Md2, M1, M3, M5, M7
       integer :: i, ii, j, jj, ni, nj, nsd, npd, ndd
       integer :: l, lk, lij, l1, l2, l3, l4, l5, l12, l34



c
c------------------------------------------
c

      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
c      do 50 i=1,natom
c      do 50 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 50   continue
c
      do 181 l=1,3
 181   Ll(l)=l*(l-1)/2
c
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2

c
c--- 2 index electron repulsion for density basis set
c
cloop  (s|s) gradients
c
c
      do 200 i=1,nsd
      do 200 j=1,i
c
      cc=-0.5D0*af(i)*af(j)
       f1=1.D0
      if (i.ne.j) then
       f1=2.D0
      endif
c
      dd=d(Nucd(i),Nucd(j))
c
      do 200 ni=1,ncontd(i)
      do 200 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
c
      u=alf*dd
      s1s=t2*FUNCT(1,u)
c
c
      ccoef=2.D0*cd(i,ni)*cd(j,nj)
      cc2=ccoef*cc*f1
c l2: different p in the p shell ( x,y,z respectively)
c
      do 205 l1=1,3
        t1=Q(l1)-r(Nucd(i),l1)
        t2=Q(l1)-r(Nucd(j),l1)
        ps=t1*s1s
        sp=t2*s1s
c
        f(Nucd(i),l1)=f(Nucd(i),l1)+ad(i,ni)*cc2*ps
        f(Nucd(j),l1)=f(Nucd(j),l1)+ad(j,nj)*cc2*sp
c
 205   continue
 200   continue
c------------------------------------------------------------
c
c (p|s) case  and gradients
c
      do 300 i=nsd+1,nsd+npd,3
      do 300 j=1,nsd
c
      dd=d(Nucd(i),Nucd(j))
c
      do 300 ni=1,ncontd(i)
      do 300 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      z2a=2.D0*ad(i,ni)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
c
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      t10=s1s/z2
      t20=(s0s-tj*s1s)/z2a

      ccoef=cd(i,ni)*cd(j,nj)
c
      do 305 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       ps=t1*s2s
       ii=i+l1-1
       cc=-ccoef*af(j)*af(ii)
c
       do 305 l2=1,3
c
       t1=Q(l2)-r(Nucd(j),l2)
       t2=Q(l2)-r(Nucd(i),l2)
       pp=t1*ps
       ds=t2*ps
c
       if (l1.eq.l2) then
       pp =pp +t10
       f(Nucd(i),l2)=f(Nucd(i),l2)-cc*s0s
       ds =ds +t20
       endif
c
       f(Nucd(j),l2)=f(Nucd(j),l2)+cc*pp*2.D0*ad(j,nj)
       f(Nucd(i),l2)=f(Nucd(i),l2)+cc*ds*2.D0*ad(i,ni)
 305  continue
c
 300  continue
c-------------------------------------------------------------------
c (p|p) case gradients
      do 400 i=nsd+1,nsd+npd,3
      do 400 j=nsd+1,i,3
c
      dd=d(Nucd(i),Nucd(j))
c
c
      do 400 ni=1,ncontd(i)
      do 400 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
c
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      t25=s2s/z2
c
      ccoef=cd(i,ni)*cd(j,nj)
c
      do 405 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pi0s=t1*s1s
       pis=t1*s2s
       t20=pis/z2
       pi2s=t1*s3s
       t15=(pi0s-ti*pis)/(2.D0*ad(j,nj))
       t17=pis/z2
c
       lij=3
       if(i.eq.j) then
        lij=l1
       endif
c
       ii=i+l1-1
       do 405 l2=1,lij
c
       t2=Q(l2)-r(Nucd(j),l2)
       spj=t2*s1s
       sp1j=t2*s2s
       t20=(spj-tj*sp1j)/(2.D0*ad(i,ni))
       t21=sp1j/z2
c
       p1p=t2*pi2s
c
       if (l1.eq.l2) then
        p1p=p1p+t25
       endif
c
       jj=j+l2-1
       f1=1.D0
      if (ii.ne.jj) then
       f1=2.D0
      endif
c
       cc=-0.5D0*af(ii)*af(jj)*f1*ccoef
c--------
       ii=i+l1-1
       jj=j+l2-1
c
c index of p
c
       do 405 l3=1,3
c
       t0=Q(l3)-r(Nucd(j),l3)
       t1=Q(l3)-r(Nucd(i),l3)
       dp=t1*p1p
       pd=t0*p1p
c
       if (l1.eq.l3) then
        dp=dp+t20
        pd=pd+t21
       f(Nucd(i),l3)=f(Nucd(i),l3)-cc*spj
       endif
c
       if (l2.eq.l3) then
        dp=dp+t17
        pd=pd+t15
       f(Nucd(j),l3)=f(Nucd(j),l3)-cc*pi0s
       endif
c
       f(Nucd(i),l3)=f(Nucd(i),l3)+cc*dp*2.D0*ad(i,ni)
       f(Nucd(j),l3)=f(Nucd(j),l3)+cc*pd*2.D0*ad(j,nj)
c
 405  continue
c
 400  continue
c
c-------------------------------------------------------------------
c (d|s)  gradients
      do 500 i=nsd+npd+1,Md,6
      do 500 j=1,nsd
c
      dd=d(Nucd(i),Nucd(j))
c
      do 500 ni=1,ncontd(i)
      do 500 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      zc=2.D0*ad(i,ni)
      roz=ad(j,nj)/zij
      alf=roz*ad(i,ni)
      t0=ad(i,ni)*ad(j,nj)
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
c
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
c
      t10=(s1s-alf*s2s/ad(i,ni))/zc
      ccoef=cd(i,ni)*cd(j,nj)
c
      do 505 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pi0s=t1*s1s
       pis=t1*s2s
       t12=pis/z2
       t22=(pi0s-roz*pis)/zc
       pi2s=t1*s3s
c
       do 505 l2=1,l1
c
       t1=Q(l2)-r(Nucd(i),l2)
       pj0s=t1*s1s
       pjs=t1*s2s
       t11=pjs/z2
       t21=(pj0s-roz*pjs)/zc
c
       ds=t1*pi2s
c
       f1=1.D0
       if (l1.eq.l2) then
        f1=sq3
        ds=ds+t10
       endif
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
c
       cc=-af(ii)*af(j)*ccoef/f1
       cc1=cc*2.D0
       cci=cc1*ad(i,ni)
       ccj=cc1*ad(j,nj)

c
c------------------------
c gradient part
       do 505 l3=1,3
c
       t0=Q(l3)-r(Nucd(j),l3)
       t1=Q(l3)-r(Nucd(i),l3)
       dp=t0*ds
       fs=t1*ds
c
       if (l1.eq.l3) then
        dp=dp+t11
        fs=fs+t21
        f(Nucd(i),l3)=f(Nucd(i),l3)-cc*pj0s
       endif
c
       if (l2.eq.l3) then
        dp=dp+t12
        fs=fs+t22
       f(Nucd(i),l3)=f(Nucd(i),l3)-cc*pi0s
       endif
c
       f(Nucd(i),l3)=f(Nucd(i),l3)+cci*fs
       f(Nucd(j),l3)=f(Nucd(j),l3)+ccj*dp
c
 505  continue
c
 500  continue
c
c-------------------------------------------------------------------
c (d|p)  gradients
      do 600 i=nsd+npd+1,Md,6
      do 600 j=nsd+1,nsd+npd,3
c
      dd=d(Nucd(i),Nucd(j))
c
      do 600 ni=1,ncontd(i)
      do 600 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
      zc=2.D0*ad(i,ni)
      zc2=2.D0*ad(j,nj)
      roz=ad(j,nj)/zij
      roz2=ad(i,ni)/zij
c
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
c
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      s4s=t2*FUNCT(4,u)
c
      t13a=s1s/z2
      t13=s2s/z2
      t10=(s0s-roz*s1s)/zc
      t11=(s1s-roz*s2s)/zc
      t12=(s2s-roz*s3s)/zc
c
      ccoef=cd(i,ni)*cd(j,nj)
c
      do 605 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pis=t1*s2s
       pi2s=t1*s3s
       pi3s=t1*s4s
       t14=pi2s/z2
c
       do 605 l2=1,l1
c
       t1=Q(l2)-r(Nucd(i),l2)
       pjs=t1*s2s
       pj2s=t1*s3s
       t15=pj2s/z2
c
       ds=t1*pis
       d1s=t1*pi2s
       d2s=t1*pi3s
c
       f1=1.D0
       if (l1.eq.l2) then
        ds=ds+t10
        d1s=d1s+t11
        d2s=d2s+t12
        f1=sq3
       endif
       t18=(ds-roz2*d1s)/zc2
       t23=d1s/z2
c
       do 605 l3=1,3
c
       t0=Q(l3)-r(Nucd(j),l3)
       dp=t0*d2s
       pip=t0*pi2s
       pi0p=t0*pis
       pjp=t0*pj2s
       pj0p=t0*pjs
c
       if (l1.eq.l3) then
        dp=dp+t15
        pip=pip+t13
        pi0p=pi0p+t13a
       endif
c
       if (l2.eq.l3) then
        dp=dp+t14
        pjp=pjp+t13
        pj0p=pj0p+t13a
       endif
c
       t16=pip/z2
       t17=pjp/z2
       t21=(pj0p-roz*pjp)/zc
       t22=(pi0p-roz*pip)/zc
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
       cc=-af(ii)*af(jj)*ccoef/f1
       cc1=cc*2.D0
       cci=cc1*ad(i,ni)
       ccj=cc1*ad(j,nj)

c
c------------------------
c gradients
c
       do 605 l4=1,3
c
       t0=Q(l4)-r(Nucd(j),l4)
       t1=Q(l4)-r(Nucd(i),l4)
       dsd=t0*dp
       fp=t1*dp
c
       if (l1.eq.l4) then
        dsd=dsd+t17
        fp=fp+t21
        f(Nucd(i),l4)=f(Nucd(i),l4)-cc*pj0p
       endif
c
       if (l2.eq.l4) then
        dsd=dsd+t16
        fp=fp+t22
        f(Nucd(i),l4)=f(Nucd(i),l4)-cc*pi0p
       endif
c
       if (l3.eq.l4) then
        f(Nucd(j),l4)=f(Nucd(j),l4)-cc*ds
        dsd=dsd+t18
        fp=fp+t23
       endif
c
       f(Nucd(i),l4)=f(Nucd(i),l4)+cci*fp
       f(Nucd(j),l4)=f(Nucd(j),l4)+ccj*dsd
 605  continue
c
 600  continue
c
c--------------------------------------------
c (d|d) gradients
      do 700 i=nsd+npd+1,Md,6
      do 700 j=nsd+npd+1,i,6
c
      dd=d(Nucd(i),Nucd(j))
c
      do 700 ni=1,ncontd(i)
      do 700 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
      zc=2.D0*ad(i,ni)
      zc2=2.D0*ad(j,nj)
      roz=ad(j,nj)/zij
      roz2=ad(i,ni)/zij
c
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
c
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      s4s=t2*FUNCT(4,u)
      s5s=t2*FUNCT(5,u)
c
      t13a=s1s/z2
      t13=s2s/z2
      t13b=s3s/z2
      t10=(s0s-roz*s1s)/zc
      t11=(s1s-roz*s2s)/zc
      t12=(s2s-roz*s3s)/zc
      t12b=(s3s-roz*s4s)/zc
c
      ccoef=cd(i,ni)*cd(j,nj)
c
      do 705 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pi0s=t1*s1s
       pis=t1*s2s
       pi2s=t1*s3s
       pi3s=t1*s4s
       pi4s=t1*s5s
       t14a=pis/z2
       t14b=pi3s/z2
       t14=pi2s/z2
       t25=(pi0s-roz2*pis)/zc2
       t26=(pis-roz2*pi2s)/zc2
c
       do 705 l2=1,l1
c
       t1=Q(l2)-r(Nucd(i),l2)
       pj0s=t1*s1s
       pjs=t1*s2s
       pj2s=t1*s3s
       pj3s=t1*s4s
       t15a=pjs/z2
       t15b=pj3s/z2
       t15=pj2s/z2
       t23=(pj0s-roz2*pjs)/zc2
       t24=(pjs-roz2*pj2s)/zc2
c
       ds=t1*pis
       d1s=t1*pi2s
       d2s=t1*pi3s
       d3s=t1*pi4s
c
       f1=1.D0
       if (l1.eq.l2) then
        ds=ds+t10
        d1s=d1s+t11
        d2s=d2s+t12
        d3s=d3s+t12b
        f1=sq3
       endif
       t18b=(d1s-roz2*d2s)/zc2
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
       do 705 l3=1,lij
c
       t0=Q(l3)-r(Nucd(j),l3)
       s1pk=t0*s2s
       s2pk=t0*s3s
       t27a=s1pk/z2
       t27=s2pk/z2
       d0p=t0*d1s
       dp=t0*d2s
       d1p=t0*d3s
c
       pi2p=t0*pi3s
       pip=t0*pi2s
       pi0p=t0*pis
       pj2p=t0*pj3s
       pjp=t0*pj2s
       pj0p=t0*pjs
c
       if (l1.eq.l3) then
        d0p=d0p+t15a
        dp=dp+t15
        d1p=d1p+t15b
        pi2p=pi2p+t13b
        pip=pip+t13
        pi0p=pi0p+t13a
       endif
c
       if (l2.eq.l3) then
        d0p=d0p+t14a
        dp=dp+t14
        d1p=d1p+t14b
        pj2p=pj2p+t13b
        pjp=pjp+t13
        pj0p=pj0p+t13a
       endif
c
       t16b=pi2p/z2
       t17b=pj2p/z2
       t33=dp/z2
       t43=(d0p-roz2*dp)/zc2
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
c
       do 705 l4=1,lk
c
       t0=Q(l4)-r(Nucd(j),l4)
       d1d=t0*d1p
       d0pl=t0*d1s
       dpl=t0*d2s
       pi0d=t0*pip
       pid=t0*pi2p
       pj0d=t0*pjp
       pjd=t0*pj2p
c
       if (l1.eq.l4) then
        d1d=d1d+t17b
        d0pl=d0pl+t15a
        dpl=dpl+t15
        pi0d=pi0d+t27a
        pid=pid+t27
       endif
c
       if (l2.eq.l4) then
        d1d=d1d+t16b
        d0pl=d0pl+t14a
        dpl=dpl+t14
        pj0d=pj0d+t27a
        pjd=pjd+t27
       endif
c
       f2=1.D0
       if (l3.eq.l4) then
        d1d=d1d+t18b
        pi0d=pi0d+t25
        pid=pid+t26
        pj0d=pj0d+t23
        pjd=pjd+t24
       f2=sq3
       endif
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
        factor=1.D0
       if (ii.ne.jj) then
        factor=2.D0
       endif
c
       t30=(pj0d-roz*pjd)/zc
       t31=(pi0d-roz*pid)/zc
       t32=dpl/z2
       t40=pjd/z2
       t41=pid/z2
       t42=(d0pl-roz2*dpl)/zc2
       cc=-0.5D0*factor*af(ii)*af(jj)*ccoef/(f1*f2)
       cc1=cc*2.D0
       cci=cc1*ad(i,ni)
       ccj=cc1*ad(j,nj)
c
c------------------------

c gradients
      do 705 l5=1,3
       t1=Q(l5)-r(Nucd(i),l5)
       t2=Q(l5)-r(Nucd(j),l5)
       df=t2*d1d
       fd=t1*d1d
c
       if (l1.eq.l5) then
        df=df+t40
        fd=fd+t30
        f(Nucd(i),l5)=f(Nucd(i),l5)-cc*pj0d
       endif
c
       if (l2.eq.l5) then
        df=df+t41
        fd=fd+t31
        f(Nucd(i),l5)=f(Nucd(i),l5)-cc*pi0d
       endif
c
       if (l3.eq.l5) then
        df=df+t42
        fd=fd+t32
        f(Nucd(j),l5)=f(Nucd(j),l5)-cc*d0pl
       endif
c
       if (l4.eq.l5) then
        df=df+t43
        fd=fd+t33
        f(Nucd(j),l5)=f(Nucd(j),l5)-cc*d0p
       endif
c
       f(Nucd(i),l5)=f(Nucd(i),l5)+cci*fd
       f(Nucd(j),l5)=f(Nucd(j),l5)+ccj*df
c
 705  continue
c
 700  continue
c
c-------------------------------------------------------------------
c
c RMM must be passed as argument if this debugging part is used
c
c debugged the s and p parts, using constant coefficients c(i)
c because the total 2 electron energy doesn't depend on these
c coefficients, but individual terms do.
c this is because d(E2)/dR = d(E2)/d C dC/dR and d(E2)/dC is 0
c by construction (variational Dunlap's fit)
c
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, also F later
      M5=M3+MM
c now G
      M7=M5+MM
c
c     do 177 i=1,nsd
c     do 178 j=1,i
c      k=i+((Md2-j)*(j-1))/2
c 178  cont=cont+af(i)*af(j)*RMM(M7+k-1)
c     do 179 j=i+1,nsd
c      k=j+((Md2-i)*(i-1))/2
c 179  cont=cont+af(i)*af(j)*RMM(M7+k-1)
c 177  continue
c
c     do 277 i=nsd+1,nsd+npd
c     do 278 j=nsd+1,i
c      k=i+((Md2-j)*(j-1))/2
c 278  cont=cont+af(i)*af(j)*RMM(M7+k-1)
c     do 279 j=i+1,nsd+npd
c      k=j+((Md2-i)*(i-1))/2
c 279  cont=cont+af(i)*af(j)*RMM(M7+k-1)
c 277  continue
c
c     do 377 i=nsd+npd+1,Md
c     do 378 j=nsd+npd+1,i
c      k=i+((Md2-j)*(j-1))/2
c 378  cont=cont+af(i)*af(j)*RMM(M7+k-1)
c     do 379 j=i+1,Md
c      k=j+((Md2-i)*(i-1))/2
c 379  cont=cont+af(i)*af(j)*RMM(M7+k-1)
c 377  continue
c
c      do 477 i=nsd+1,nsd+npd
c      do 477 j=1,nsd
c      k=i+((Md2-j)*(j-1))/2
c477    cont=cont+2.D0*af(i)*af(j)*RMM(M7+k-1)
c
c      do 577 i=nsd+npd+1,Md
c      do 577 j=1,nsd
c      k=i+((Md2-j)*(j-1))/2
c577    cont=cont+2.D0*af(i)*af(j)*RMM(M7+k-1)
c
c      do 677 i=nsd+npd+1,Md
c      do 677 j=nsd+1,nsd+npd
c      k=i+((Md2-j)*(j-1))/2
c677    cont=cont+2.D0*af(i)*af(j)*RMM(M7+k-1)
c
c       cont=-0.5D0*cont
c     write(*,*) 'cont =',cont
c     do i=1,natom
c      write(*,*) i,f(i,1),f(i,2),f(i,3)
c     enddo
c-------------------------------------------------------------------
      return
      end subroutine
      end module subm_int2G
