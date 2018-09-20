! Integrals subroutine -Second part
! Gradients only
! 2 e integrals, 2 index : density fitting functions
! All of them are calculated
! using the Obara-Saika recursive method.
       module subm_int2G; contains
       subroutine int2G(f)

       use liotemp   , only: FUNCT
       use garcha_mod, only: RMM, ll, natom, M, Md, NORM, pi5, r
     >                     , af, nshelld, ad, nucd, ncontd, d, cd

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

      SQ3 = 1.0D0
      if (NORM) SQ3 = sqrt(3.D0)

      do l=1,3
       Ll(l)=l*(l-1)/2
      enddo

      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2

      ! (s|s)
      do i=1,nsd
      do j=1,i
      cc=-0.5D0*af(i)*af(j)
       f1=1.D0
      if (i.ne.j) then
       f1=2.D0
      endif

      dd=d(Nucd(i),Nucd(j))

      do ni=1,ncontd(i)
      do nj=1,ncontd(j)
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
      u=alf*dd
      s1s=t2*FUNCT(1,u)

      ccoef=2.D0*cd(i,ni)*cd(j,nj)
      cc2=ccoef*cc*f1
      do l1=1,3
        t1=Q(l1)-r(Nucd(i),l1)
        t2=Q(l1)-r(Nucd(j),l1)
        ps=t1*s1s
        sp=t2*s1s
        f(Nucd(i),l1)=f(Nucd(i),l1)+ad(i,ni)*cc2*ps
        f(Nucd(j),l1)=f(Nucd(j),l1)+ad(j,nj)*cc2*sp
      enddo
      enddo
      enddo
      enddo
      enddo

      ! (p|s) case  and gradients
      do i=nsd+1,nsd+npd,3
      do j=1,nsd
      dd=d(Nucd(i),Nucd(j))
      do ni=1,ncontd(i)
      do nj=1,ncontd(j)
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      z2a=2.D0*ad(i,ni)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      t10=s1s/z2
      t20=(s0s-tj*s1s)/z2a

      ccoef=cd(i,ni)*cd(j,nj)
      do l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       ps=t1*s2s
       ii=i+l1-1
       cc=-ccoef*af(j)*af(ii)
       do l2=1,3
       t1=Q(l2)-r(Nucd(j),l2)
       t2=Q(l2)-r(Nucd(i),l2)
       pp=t1*ps
       ds=t2*ps
       if (l1.eq.l2) then
       pp =pp +t10
       f(Nucd(i),l2)=f(Nucd(i),l2)-cc*s0s
       ds =ds +t20
       endif
       f(Nucd(j),l2)=f(Nucd(j),l2)+cc*pp*2.D0*ad(j,nj)
       f(Nucd(i),l2)=f(Nucd(i),l2)+cc*ds*2.D0*ad(i,ni)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

      ! (p|p) case gradients
      do i=nsd+1,nsd+npd,3
      do j=nsd+1,i,3
      dd=d(Nucd(i),Nucd(j))
      do ni=1,ncontd(i)
      do nj=1,ncontd(j)
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      t25=s2s/z2
      ccoef=cd(i,ni)*cd(j,nj)
      do l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pi0s=t1*s1s
       pis=t1*s2s
       t20=pis/z2
       pi2s=t1*s3s
       t15=(pi0s-ti*pis)/(2.D0*ad(j,nj))
       t17=pis/z2
       lij=3
       if(i.eq.j) then
        lij=l1
       endif
       ii=i+l1-1
       do l2=1,lij
       t2=Q(l2)-r(Nucd(j),l2)
       spj=t2*s1s
       sp1j=t2*s2s
       t20=(spj-tj*sp1j)/(2.D0*ad(i,ni))
       t21=sp1j/z2
       p1p=t2*pi2s
       if (l1.eq.l2) then
        p1p=p1p+t25
       endif
       jj=j+l2-1
       f1=1.D0
      if (ii.ne.jj) then
       f1=2.D0
      endif
       cc=-0.5D0*af(ii)*af(jj)*f1*ccoef

       ii=i+l1-1
       jj=j+l2-1

       do l3=1,3
       t0=Q(l3)-r(Nucd(j),l3)
       t1=Q(l3)-r(Nucd(i),l3)
       dp=t1*p1p
       pd=t0*p1p
       if (l1.eq.l3) then
        dp=dp+t20
        pd=pd+t21
       f(Nucd(i),l3)=f(Nucd(i),l3)-cc*spj
       endif
       if (l2.eq.l3) then
        dp=dp+t17
        pd=pd+t15
       f(Nucd(j),l3)=f(Nucd(j),l3)-cc*pi0s
       endif
       f(Nucd(i),l3)=f(Nucd(i),l3)+cc*dp*2.D0*ad(i,ni)
       f(Nucd(j),l3)=f(Nucd(j),l3)+cc*pd*2.D0*ad(j,nj)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

      ! (d|s)  gradients
      do i=nsd+npd+1,Md,6
      do j=1,nsd
      dd=d(Nucd(i),Nucd(j))
      do ni=1,ncontd(i)
      do nj=1,ncontd(j)
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      zc=2.D0*ad(i,ni)
      roz=ad(j,nj)/zij
      alf=roz*ad(i,ni)
      t0=ad(i,ni)*ad(j,nj)
      t1=sqrt(zij)*t0
      t2=pi5/t1
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      t10=(s1s-alf*s2s/ad(i,ni))/zc
      ccoef=cd(i,ni)*cd(j,nj)
      do l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pi0s=t1*s1s
       pis=t1*s2s
       t12=pis/z2
       t22=(pi0s-roz*pis)/zc
       pi2s=t1*s3s
       do l2=1,l1
       t1=Q(l2)-r(Nucd(i),l2)
       pj0s=t1*s1s
       pjs=t1*s2s
       t11=pjs/z2
       t21=(pj0s-roz*pjs)/zc
       ds=t1*pi2s
       f1=1.D0
       if (l1.eq.l2) then
        f1=sq3
        ds=ds+t10
       endif
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       cc=-af(ii)*af(j)*ccoef/f1
       cc1=cc*2.D0
       cci=cc1*ad(i,ni)
       ccj=cc1*ad(j,nj)

       do l3=1,3
       t0=Q(l3)-r(Nucd(j),l3)
       t1=Q(l3)-r(Nucd(i),l3)
       dp=t0*ds
       fs=t1*ds
       if (l1.eq.l3) then
        dp=dp+t11
        fs=fs+t21
        f(Nucd(i),l3)=f(Nucd(i),l3)-cc*pj0s
       endif
       if (l2.eq.l3) then
        dp=dp+t12
        fs=fs+t22
       f(Nucd(i),l3)=f(Nucd(i),l3)-cc*pi0s
       endif
       f(Nucd(i),l3)=f(Nucd(i),l3)+cci*fs
       f(Nucd(j),l3)=f(Nucd(j),l3)+ccj*dp
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

      ! (d|p)  gradients
      do i=nsd+npd+1,Md,6
      do j=nsd+1,nsd+npd,3
      dd=d(Nucd(i),Nucd(j))
      do ni=1,ncontd(i)
      do nj=1,ncontd(j)
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
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      s4s=t2*FUNCT(4,u)
      t13a=s1s/z2
      t13=s2s/z2
      t10=(s0s-roz*s1s)/zc
      t11=(s1s-roz*s2s)/zc
      t12=(s2s-roz*s3s)/zc
      ccoef=cd(i,ni)*cd(j,nj)
      do l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pis=t1*s2s
       pi2s=t1*s3s
       pi3s=t1*s4s
       t14=pi2s/z2
       do l2=1,l1
       t1=Q(l2)-r(Nucd(i),l2)
       pjs=t1*s2s
       pj2s=t1*s3s
       t15=pj2s/z2
       ds=t1*pis
       d1s=t1*pi2s
       d2s=t1*pi3s
       f1=1.D0
       if (l1.eq.l2) then
        ds=ds+t10
        d1s=d1s+t11
        d2s=d2s+t12
        f1=sq3
       endif
       t18=(ds-roz2*d1s)/zc2
       t23=d1s/z2
       do l3=1,3
       t0=Q(l3)-r(Nucd(j),l3)
       dp=t0*d2s
       pip=t0*pi2s
       pi0p=t0*pis
       pjp=t0*pj2s
       pj0p=t0*pjs
       if (l1.eq.l3) then
        dp=dp+t15
        pip=pip+t13
        pi0p=pi0p+t13a
       endif
       if (l2.eq.l3) then
        dp=dp+t14
        pjp=pjp+t13
        pj0p=pj0p+t13a
       endif
       t16=pip/z2
       t17=pjp/z2
       t21=(pj0p-roz*pjp)/zc
       t22=(pi0p-roz*pip)/zc
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
       cc=-af(ii)*af(jj)*ccoef/f1
       cc1=cc*2.D0
       cci=cc1*ad(i,ni)
       ccj=cc1*ad(j,nj)

       do l4=1,3
       t0=Q(l4)-r(Nucd(j),l4)
       t1=Q(l4)-r(Nucd(i),l4)
       dsd=t0*dp
       fp=t1*dp
       if (l1.eq.l4) then
        dsd=dsd+t17
        fp=fp+t21
        f(Nucd(i),l4)=f(Nucd(i),l4)-cc*pj0p
       endif
       if (l2.eq.l4) then
        dsd=dsd+t16
        fp=fp+t22
        f(Nucd(i),l4)=f(Nucd(i),l4)-cc*pi0p
       endif
       if (l3.eq.l4) then
        f(Nucd(j),l4)=f(Nucd(j),l4)-cc*ds
        dsd=dsd+t18
        fp=fp+t23
       endif
       f(Nucd(i),l4)=f(Nucd(i),l4)+cci*fp
       f(Nucd(j),l4)=f(Nucd(j),l4)+ccj*dsd
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

      ! (d|d) gradients
      do i=nsd+npd+1,Md,6
      do j=nsd+npd+1,i,6
      dd=d(Nucd(i),Nucd(j))
      do ni=1,ncontd(i)
      do nj=1,ncontd(j)
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
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      s4s=t2*FUNCT(4,u)
      s5s=t2*FUNCT(5,u)
      t13a=s1s/z2
      t13=s2s/z2
      t13b=s3s/z2
      t10=(s0s-roz*s1s)/zc
      t11=(s1s-roz*s2s)/zc
      t12=(s2s-roz*s3s)/zc
      t12b=(s3s-roz*s4s)/zc
      ccoef=cd(i,ni)*cd(j,nj)
      do l1=1,3
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
       do l2=1,l1
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
       ds=t1*pis
       d1s=t1*pi2s
       d2s=t1*pi3s
       d3s=t1*pi4s
       f1=1.D0
       if (l1.eq.l2) then
        ds=ds+t10
        d1s=d1s+t11
        d2s=d2s+t12
        d3s=d3s+t12b
        f1=sq3
       endif
       t18b=(d1s-roz2*d2s)/zc2
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
       do l3=1,lij
       t0=Q(l3)-r(Nucd(j),l3)
       s1pk=t0*s2s
       s2pk=t0*s3s
       t27a=s1pk/z2
       t27=s2pk/z2
       d0p=t0*d1s
       dp=t0*d2s
       d1p=t0*d3s
       pi2p=t0*pi3s
       pip=t0*pi2s
       pi0p=t0*pis
       pj2p=t0*pj3s
       pjp=t0*pj2s
       pj0p=t0*pjs
       if (l1.eq.l3) then
        d0p=d0p+t15a
        dp=dp+t15
        d1p=d1p+t15b
        pi2p=pi2p+t13b
        pip=pip+t13
        pi0p=pi0p+t13a
       endif
       if (l2.eq.l3) then
        d0p=d0p+t14a
        dp=dp+t14
        d1p=d1p+t14b
        pj2p=pj2p+t13b
        pjp=pjp+t13
        pj0p=pj0p+t13a
       endif
       t16b=pi2p/z2
       t17b=pj2p/z2
       t33=dp/z2
       t43=(d0p-roz2*dp)/zc2
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif

       do l4=1,lk
       t0=Q(l4)-r(Nucd(j),l4)
       d1d=t0*d1p
       d0pl=t0*d1s
       dpl=t0*d2s
       pi0d=t0*pip
       pid=t0*pi2p
       pj0d=t0*pjp
       pjd=t0*pj2p
       if (l1.eq.l4) then
        d1d=d1d+t17b
        d0pl=d0pl+t15a
        dpl=dpl+t15
        pi0d=pi0d+t27a
        pid=pid+t27
       endif
       if (l2.eq.l4) then
        d1d=d1d+t16b
        d0pl=d0pl+t14a
        dpl=dpl+t14
        pj0d=pj0d+t27a
        pjd=pjd+t27
       endif
       f2=1.D0
       if (l3.eq.l4) then
        d1d=d1d+t18b
        pi0d=pi0d+t25
        pid=pid+t26
        pj0d=pj0d+t23
        pjd=pjd+t24
       f2=sq3
       endif
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
        factor=1.D0
       if (ii.ne.jj) then
        factor=2.D0
       endif
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
      do l5=1,3
       t1=Q(l5)-r(Nucd(i),l5)
       t2=Q(l5)-r(Nucd(j),l5)
       df=t2*d1d
       fd=t1*d1d
       if (l1.eq.l5) then
        df=df+t40
        fd=fd+t30
        f(Nucd(i),l5)=f(Nucd(i),l5)-cc*pj0d
       endif
       if (l2.eq.l5) then
        df=df+t41
        fd=fd+t31
        f(Nucd(i),l5)=f(Nucd(i),l5)-cc*pi0d
       endif
       if (l3.eq.l5) then
        df=df+t42
        fd=fd+t32
        f(Nucd(j),l5)=f(Nucd(j),l5)-cc*d0pl
       endif
       if (l4.eq.l5) then
        df=df+t43
        fd=fd+t33
        f(Nucd(j),l5)=f(Nucd(j),l5)-cc*d0p
       endif
       f(Nucd(i),l5)=f(Nucd(i),l5)+cci*fd
       f(Nucd(j),l5)=f(Nucd(j),l5)+ccj*df
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

      return
      end subroutine
      end module subm_int2G
