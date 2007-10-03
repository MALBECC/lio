c-------------------------------------------------------------------
c Calculation of electrical potential and
c electostatical field (its gradient) for an arbitrary point of
c the space
c 16 Feb. 1994 - Dario Estrin
c 1 e integrals
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
c Input : basis function information
c-------------------------------------------------------------------
      subroutine efield(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,RMM,xi,V,Ef)
c
      implicit real*8 (a-h,o-z)
      INCLUDE 'param'
      logical NORM
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0)
      parameter(rmax=30.0D0)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(natom)
      dimension r(nt,3),nshell(0:3)
      dimension RMM(*),xi(3),Ll(3),Ef(3),x0x(3)
c
      dimension x1x(3),x2x(3),x3x(3),x4x(3),x5x(3)
      dimension dn(3),dn1(3),dn2(3),dn3(3),dn4(3),dn5(3)
      dimension dn6(3),dn7(3),dn8(3),dn9(3),dn10(3)
      dimension dn2b(3),dn4b(3),dn5b(3),dn7b(3),dn8b(3),dn9b(3)
      COMMON /TABLE/ STR(880,0:21)
c auxiliar quantities
c
      dimension Q(3),d(ntq,ntq)
c distance between pairs of centers
c
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2
c
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      M2=2*M
c
      M1=1
c now Pnew
      M3=M1+MM
c now S, F also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
      M11=M9+MMd
c
      Ef(1)=0.0D0
      Ef(2)=0.0D0
      Ef(3)=0.0D0
c
       V=0.0D0
       do 11 n=1,natom
       tx=xi(1)-r(n,1)
       ty=xi(2)-r(n,2)
       tz=xi(3)-r(n,3)
       dd1=tx**2 + ty**2 + tz**2
       dd=sqrt(dd1)
       tv=Iz(n)/dd
       V=V-tv
c
       te=tv/dd1
       Ef(1)=Ef(1)+ te*tx
       Ef(2)=Ef(2)+ te*ty
       Ef(3)=Ef(3)+ te*tz
c
  11   continue
c
      do 50 i=1,natom
      do 50 j=1,natom
       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
 50   continue
c
c first loop (s|s) case -------------------------------------------
c
      do 200 i=1,ns
      do 200 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
c
      do 201 ni=1,ncont(i)
      do 201 nj=1,ncont(j)
c
c (0|0) calculation
      zij=a(i,ni)+a(j,nj)
      z2=2.0D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      kk=i+(M2-j)*(j-1)/2
      rexp=alf*dd
      if (rexp.gt.rmax) goto 201
c
      t1=a(i,ni)/zij
      t2=a(j,nj)/zij
      Q(1)=t1*r(Nuc(i),1)+t2*r(Nuc(j),1)
      Q(2)=t1*r(Nuc(i),2)+t2*r(Nuc(j),2)
      Q(3)=t1*r(Nuc(i),3)+t2*r(Nuc(j),3)
      ccoef=c(i,ni)*c(j,nj)
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=ss*2.D0*sqrt(zij/pi)
c
       q1=Q(1)-xi(1)
       q2=Q(2)-xi(2)
       q3=Q(3)-xi(3)
       u=q1**2+q2**2+q3**2
       u=u*zij
       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
c gradients
       temp=z2*s1s
       x0x(1)=temp*q1
       x0x(2)=temp*q2
       x0x(3)=temp*q3
c
      cc=ccoef*RMM(kk)
      term=cc*s0s
      V=V+term
c
c gradients
      Ef(1)=Ef(1)+cc*x0x(1)
      Ef(2)=Ef(2)+cc*x0x(2)
      Ef(3)=Ef(3)+cc*x0x(3)
 202  continue
 201  continue
 200  continue
c
c------------------------------------------------------------------
c
c second loop  (p|s) case
c
c
      do 300 i=ns+1,ns+np,3
      do 300 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
c
      do 300 ni=1,ncont(i)
      do 300 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.0D0*zij

      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.gt.rmax) goto 301
c
      t1=a(i,ni)/zij
      t2=a(j,nj)/zij
      Q(1)=t1*r(Nuc(i),1)+t2*r(Nuc(j),1)
      Q(2)=t1*r(Nuc(i),2)+t2*r(Nuc(j),2)
      Q(3)=t1*r(Nuc(i),3)+t2*r(Nuc(j),3)
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
      ccoef=c(i,ni)*c(j,nj)
c
       q1=Q(1)-xi(1)
       q2=Q(2)-xi(2)
       q3=Q(3)-xi(3)
       u=q1**2+q2**2+q3**2
       u=u*zij
c
       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
c
        temp=z2*s1s
       x0x(1)=temp*q1
       x0x(2)=temp*q2
       x0x(3)=temp*q3
        temp=z2*s2s
       x1x(1)=temp*q1
       x1x(2)=temp*q2
       x1x(3)=temp*q3
c
c
c l2: different p in the p shell ( x,y,z respectively)
c
      do 305 l2=1,3
        t1=Q(l2)-r(Nuc(i),l2)
        t2=Q(l2)-xi(l2)
        ii=i+l2-1
c ii index , taking into account different components of the shell
c
        kk=ii+((M2-j)*(j-1))/2
c
        tna=t1*s0s-(Q(l2)-xi(l2))*s1s
c gradients
       dn(1)=t1*x0x(1)-t2*x1x(1)
       dn(2)=t1*x0x(2)-t2*x1x(2)
       dn(3)=t1*x0x(3)-t2*x1x(3)
c
       dn(l2)=dn(l2)+s1s
c


        kk=ii+((M2-j)*(j-1))/2
        te=ccoef*RMM(kk)
        term= te*tna
c
        V=V+term
        Ef(1)=Ef(1)+te*dn(1)
        Ef(2)=Ef(2)+te*dn(2)
        Ef(3)=Ef(3)+te*dn(3)
c
 305    continue
c
 301   continue

 300  continue       
c-------------------------------------------------------------------
c 
c (p|p) case
c
      do 400 i=ns+1,ns+np,3
      do 400 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 400 ni=1,ncont(i)
      do 400 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.gt.rmax) goto 401
c
      t1=a(i,ni)/zij
      t2=a(j,nj)/zij
      Q(1)=t1*r(Nuc(i),1)+t2*r(Nuc(j),1)
      Q(2)=t1*r(Nuc(i),2)+t2*r(Nuc(j),2)
      Q(3)=t1*r(Nuc(i),3)+t2*r(Nuc(j),3)
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ccoef=c(i,ni)*c(j,nj)
       temp=2.D0*sqrt(zij/pi)*ss
c
       q1=Q(1)-xi(1)
       q2=Q(2)-xi(2)
       q3=Q(3)-xi(3)
       u=q1**2+q2**2+q3**2
       u=u*zij

       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
       s3s=temp*FUNCT(3,u)
c
c
        temp=z2*s1s
       x0x(1)=temp*q1
       x0x(2)=temp*q2
       x0x(3)=temp*q3
        temp=z2*s2s
       x1x(1)=temp*q1
       x1x(2)=temp*q2
       x1x(3)=temp*q3
        temp=z2*s3s
       x2x(1)=temp*q1
       x2x(2)=temp*q2
       x2x(3)=temp*q3
c
      t26=(x0x(1)-x1x(1))/z2
      t27=(x0x(2)-x1x(2))/z2
      t28=(x0x(3)-x1x(3))/z2
c
c
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-xi(l1)
       p0s=t1*s0s-t2*s1s
       p1s=t1*s1s-t2*s2s
c
c dn(u) (pi|Au|s)
      dn(1)=t1*x0x(1)-t2*x1x(1)
      dn(2)=t1*x0x(2)-t2*x1x(2)
      dn(3)=t1*x0x(3)-t2*x1x(3)
      dn(l1)=dn(l1)+s1s
c
      dn1(1)=t1*x1x(1)-t2*x2x(1)
      dn1(2)=t1*x1x(2)-t2*x2x(2)
      dn1(3)=t1*x1x(3)-t2*x2x(3)
      dn1(l1)=dn1(l1)+s2s
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
       do 406 l2=1,lij
c
       t1=Q(l2)-r(Nuc(j),l2)
       t2=Q(l2)-xi(l2)
       tna=t1*p0s-t2*p1s
c
       dn2(1)=t1*dn(1)-t2*dn1(1)
       dn2(2)=t1*dn(2)-t2*dn1(2)
       dn2(3)=t1*dn(3)-t2*dn1(3)
       dn2(l2)=dn2(l2)+p1s
c
       if (l1.eq.l2) then
        tna=tna+(s0s-s1s)/z2
       dn2(1)=dn2(1)+t26
       dn2(2)=dn2(2)+t27
       dn2(3)=dn2(3)+t28
       endif
c
c
       ii=i+l1-1
       jj=j+l2-1
       kk=ii+((M2-jj)*(jj-1))/2
       te=ccoef*RMM(kk)
       term=tna*te
        V=V+term
c
       Ef(1)=Ef(1)+te*dn2(1)
       Ef(2)=Ef(2)+te*dn2(2)
       Ef(3)=Ef(3)+te*dn2(3)
 406  continue
 401  continue
 400  continue
c-------------------------------------------------------------------
c (d|s) case
c
      do 500 i=ns+np+1,M,6
      do 500 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
c
      do 500 ni=1,ncont(i)
      do 500 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.gt.rmax) goto 501
c
      t1=a(i,ni)/zij
      t2=a(j,nj)/zij
      Q(1)=t1*r(Nuc(i),1)+t2*r(Nuc(j),1)
      Q(2)=t1*r(Nuc(i),2)+t2*r(Nuc(j),2)
      Q(3)=t1*r(Nuc(i),3)+t2*r(Nuc(j),3)
c
      ccoef=c(i,ni)*c(j,nj)
       ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c

       q1=Q(1)-xi(1)
       q2=Q(2)-xi(2)
       q3=Q(3)-xi(3)
       u=q1**2+q2**2+q3**2
       u=u*zij
c
       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
       s3s=temp*FUNCT(3,u)
c
        temp=z2*s1s
       x0x(1)=temp*q1
       x0x(2)=temp*q2
       x0x(3)=temp*q3
        temp=z2*s2s
       x1x(1)=temp*q1
       x1x(2)=temp*q2
       x1x(3)=temp*q3
        temp=z2*s3s
       x2x(1)=temp*q1
       x2x(2)=temp*q2
       x2x(3)=temp*q3
c
      t7=(s0s-s1s)/z2
      t8=(s1s-s2s)/z2
      t26=(x0x(1)-x1x(1))/z2
      t27=(x0x(2)-x1x(2))/z2
      t28=(x0x(3)-x1x(3))/z2
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-xi(l1)
       p0s=t1*s0s-t2*s1s
       p1s=t1*s1s-t2*s2s
c
c dn(u) (pi|Au|s)
      dn(1)=t1*x0x(1)-t2*x1x(1)
      dn(2)=t1*x0x(2)-t2*x1x(2)
      dn(3)=t1*x0x(3)-t2*x1x(3)
      dn(l1)=dn(l1)+s1s
c
      dn1(1)=t1*x1x(1)-t2*x2x(1)
      dn1(2)=t1*x1x(2)-t2*x2x(2)
      dn1(3)=t1*x1x(3)-t2*x2x(3)
      dn1(l1)=dn1(l1)+s2s
c
c
      do 506 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-xi(l2)
       tna=t1*p0s-t2*p1s
c
       dn2(1)=t1*dn(1)-t2*dn1(1)
       dn2(2)=t1*dn(2)-t2*dn1(2)
       dn2(3)=t1*dn(3)-t2*dn1(3)
       dn2(l2)=dn2(l2)+p1s
c
       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+(s0s-s1s)/z2
        dn2(1)=dn2(1)+t26
        dn2(2)=dn2(2)+t27
        dn2(3)=dn2(3)+t28
        f1=sq3
       endif
c
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       ii=i+l12-1
c
       kk=ii+((M2-j)*(j-1))/2
       cc=ccoef/f1
       te=cc*RMM(kk)
       term=te*tna
        V=V+term
c
        Ef(1)=Ef(1)+te*dn2(1)
        Ef(2)=Ef(2)+te*dn2(2)
        Ef(3)=Ef(3)+te*dn2(3)
 506  continue
c
 501  continue
 500  continue
c-----------------------------------------------------------------
c
c (d|p) case
      do 600 i=ns+np+1,M,6
      do 600 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 600 ni=1,ncont(i)
      do 600 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.gt.rmax) goto 601
c
      t1=a(i,ni)/zij
      t2=a(j,nj)/zij
      Q(1)=t1*r(Nuc(i),1)+t2*r(Nuc(j),1)
      Q(2)=t1*r(Nuc(i),2)+t2*r(Nuc(j),2)
      Q(3)=t1*r(Nuc(i),3)+t2*r(Nuc(j),3)
c
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
      ccoef=c(i,ni)*c(j,nj)
c
       q1=Q(1)-xi(1)
       q2=Q(2)-xi(2)
       q3=Q(3)-xi(3)
       u=q1**2+q2**2+q3**2
       u=u*zij
c
       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
       s3s=temp*FUNCT(3,u)
       s4s=temp*FUNCT(4,u)
c
        temp=z2*s1s
       x0x(1)=temp*q1
       x0x(2)=temp*q2
       x0x(3)=temp*q3
        temp=z2*s2s
       x1x(1)=temp*q1
       x1x(2)=temp*q2
       x1x(3)=temp*q3
        temp=z2*s3s
       x2x(1)=temp*q1
       x2x(2)=temp*q2
       x2x(3)=temp*q3
        temp=z2*s4s
       x3x(1)=temp*q1
       x3x(2)=temp*q2
       x3x(3)=temp*q3
c
      t7=(s0s-s1s)/z2
      t8=(s1s-s2s)/z2
      t9=(s2s-s3s)/z2
      t26=(x0x(1)-x1x(1))/z2
      t27=(x0x(2)-x1x(2))/z2
      t28=(x0x(3)-x1x(3))/z2
      t29=(x1x(1)-x2x(1))/z2
      t30=(x1x(2)-x2x(2))/z2
      t31=(x1x(3)-x2x(3))/z2
c
c
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-xi(l1)
       p0s=t1*s0s-t2*s1s
       p1s=t1*s1s-t2*s2s
       p2s=t1*s2s-t2*s3s
c
c dn(u) (pi|Au|s)
      dn(1)=t1*x0x(1)-t2*x1x(1)
      dn(2)=t1*x0x(2)-t2*x1x(2)
      dn(3)=t1*x0x(3)-t2*x1x(3)
      dn(l1)=dn(l1)+s1s
c
      dn1(1)=t1*x1x(1)-t2*x2x(1)
      dn1(2)=t1*x1x(2)-t2*x2x(2)
      dn1(3)=t1*x1x(3)-t2*x2x(3)
      dn1(l1)=dn1(l1)+s2s
c
       t51=(dn(1)-dn1(1))/z2
       t52=(dn(2)-dn1(2))/z2
       t53=(dn(3)-dn1(3))/z2
c
      dn2(1)=t1*x2x(1)-t2*x3x(1)
      dn2(2)=t1*x2x(2)-t2*x3x(2)
      dn2(3)=t1*x2x(3)-t2*x3x(3)
      dn2(l1)=dn2(l1)+s3s
c
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-xi(l2)
       pj0s=t1*s0s-t2*s1s
       pj1s=t1*s1s-t2*s2s
c
c dn3 (dij || s) order 0
       dn3(1)=t1*dn(1)-t2*dn1(1)
       dn3(2)=t1*dn(2)-t2*dn1(2)
       dn3(3)=t1*dn(3)-t2*dn1(3)
       dn3(l2)=dn3(l2)+p1s
c
c dn4 (dij || s) order 1
       dn4(1)=t1*dn1(1)-t2*dn2(1)
       dn4(2)=t1*dn1(2)-t2*dn2(2)
       dn4(3)=t1*dn1(3)-t2*dn2(3)
       dn4(l2)=dn4(l2)+p2s
c dn6 and dn7 used for (pj | s) order 0 and 1
c
      dn6(1)=t1*x0x(1)-t2*x1x(1)
      dn6(2)=t1*x0x(2)-t2*x1x(2)
      dn6(3)=t1*x0x(3)-t2*x1x(3)
      dn6(l2)=dn6(l2)+s1s
c
      dn7(1)=t1*x1x(1)-t2*x2x(1)
      dn7(2)=t1*x1x(2)-t2*x2x(2)
      dn7(3)=t1*x1x(3)-t2*x2x(3)
      dn7(l2)=dn7(l2)+s2s

       t54=(dn6(1)-dn7(1))/z2
       t55=(dn6(2)-dn7(2))/z2
       t56=(dn6(3)-dn7(3))/z2

       f1=1.D0
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s-s1s)/z2
        d1s=d1s+(s1s-s2s)/z2
        dn3(1)=dn3(1)+t26
        dn3(2)=dn3(2)+t27
        dn3(3)=dn3(3)+t28
        dn4(1)=dn4(1)+t29
        dn4(2)=dn4(2)+t30
        dn4(3)=dn4(3)+t31
c
       endif
c
      do 606 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-xi(l3)
       tna=t1*d0s-t2*d1s
c
       dn5(1)=t1*dn3(1)-t2*dn4(1)
       dn5(2)=t1*dn3(2)-t2*dn4(2)
       dn5(3)=t1*dn3(3)-t2*dn4(3)
       dn5(l3)=dn5(l3)+d1s
c
       if (l1.eq.l3) then
        tna=tna+(pj0s-pj1s)/z2
       dn5(1)=dn5(1)+t54
       dn5(2)=dn5(2)+t55
       dn5(3)=dn5(3)+t56
       endif
c
       if (l2.eq.l3) then
        tna=tna+(p0s-p1s)/z2
       dn5(1)=dn5(1)+t51
       dn5(2)=dn5(2)+t52
       dn5(3)=dn5(3)+t53
       endif
c
c
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
      kk=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
       te=cc*RMM(kk)
       term=te*tna
        V=V+term
c
       Ef(1)=Ef(1)+te*dn5(1)
       Ef(2)=Ef(2)+te*dn5(2)
       Ef(3)=Ef(3)+te*dn5(3)
 606  continue
 601  continue
 600  continue
c------------------------------------------------------------------
c
c (d|d) case
c
      do 700 i=ns+np+1,M,6
      do 700 j=ns+np+1,i,6
c
      dd=d(Nuc(i),Nuc(j))
c
      do 700 ni=1,ncont(i)
      do 700 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.gt.rmax) goto 701
c
      t1=a(i,ni)/zij
      t2=a(j,nj)/zij
      Q(1)=t1*r(Nuc(i),1)+t2*r(Nuc(j),1)
      Q(2)=t1*r(Nuc(i),2)+t2*r(Nuc(j),2)
      Q(3)=t1*r(Nuc(i),3)+t2*r(Nuc(j),3)
c
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
      ccoef=c(i,ni)*c(j,nj)
c
       q1=Q(1)-xi(1)
       q2=Q(2)-xi(2)
       q3=Q(3)-xi(3)
       u=q1**2+q2**2+q3**2
       u=u*zij
c

       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
       s3s=temp*FUNCT(3,u)
       s4s=temp*FUNCT(4,u)
       s5s=temp*FUNCT(5,u)
c
        temp=z2*s1s
       x0x(1)=temp*q1
       x0x(2)=temp*q2
       x0x(3)=temp*q3
        temp=z2*s2s
       x1x(1)=temp*q1
       x1x(2)=temp*q2
       x1x(3)=temp*q3
        temp=z2*s3s
       x2x(1)=temp*q1
       x2x(2)=temp*q2
       x2x(3)=temp*q3
        temp=z2*s4s
       x3x(1)=temp*q1
       x3x(2)=temp*q2
       x3x(3)=temp*q3
        temp=z2*s5s
       x4x(1)=temp*q1
       x4x(2)=temp*q2
       x4x(3)=temp*q3
c
      t50=(s0s-s1s)/z2
      t51=(s1s-s2s)/z2
      t52=(s2s-s3s)/z2
      t53=(s3s-s4s)/z2
c
      t26=(x0x(1)-x1x(1))/z2
      t27=(x0x(2)-x1x(2))/z2
      t28=(x0x(3)-x1x(3))/z2
      t29=(x1x(1)-x2x(1))/z2
      t30=(x1x(2)-x2x(2))/z2
      t31=(x1x(3)-x2x(3))/z2
      t32=(x2x(1)-x3x(1))/z2
      t33=(x2x(2)-x3x(2))/z2
      t34=(x2x(3)-x3x(3))/z2
c
c
      do 706 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-xi(l1)
       p0s=t1*s0s-t2*s1s
       p1s=t1*s1s-t2*s2s
       p2s=t1*s2s-t2*s3s
       p3s=t1*s3s-t2*s4s
c
cdn(u) (pi|Au|s)
      dn(1)=t1*x0x(1)-t2*x1x(1)
      dn(2)=t1*x0x(2)-t2*x1x(2)
      dn(3)=t1*x0x(3)-t2*x1x(3)
      dn(l1)=dn(l1)+s1s
c
      dn1(1)=t1*x1x(1)-t2*x2x(1)
      dn1(2)=t1*x1x(2)-t2*x2x(2)
      dn1(3)=t1*x1x(3)-t2*x2x(3)
      dn1(l1)=dn1(l1)+s2s
c
       t81=(dn(1)-dn1(1))/z2
       t82=(dn(2)-dn1(2))/z2
       t83=(dn(3)-dn1(3))/z2
c
      dn2(1)=t1*x2x(1)-t2*x3x(1)
      dn2(2)=t1*x2x(2)-t2*x3x(2)
      dn2(3)=t1*x2x(3)-t2*x3x(3)
      dn2(l1)=dn2(l1)+s3s
c
      dn2b(1)=t1*x3x(1)-t2*x4x(1)
      dn2b(2)=t1*x3x(2)-t2*x4x(2)
      dn2b(3)=t1*x3x(3)-t2*x4x(3)
      dn2b(l1)=dn2b(l1)+s4s
c
       t81b=(dn1(1)-dn2(1))/z2
       t82b=(dn1(2)-dn2(2))/z2
       t83b=(dn1(3)-dn2(3))/z2
c
      do 706 l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-xi(l2)
       pj0s=t1*s0s-t2*s1s
       pj1s=t1*s1s-t2*s2s
       pj2s=t1*s2s-t2*s3s
c
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
c
c dn3 (dij || s) order 0
       dn3(1)=t1*dn(1)-t2*dn1(1)
       dn3(2)=t1*dn(2)-t2*dn1(2)
       dn3(3)=t1*dn(3)-t2*dn1(3)
       dn3(l2)=dn3(l2)+p1s
c
c dn4 (dij || s) order 1
       dn4(1)=t1*dn1(1)-t2*dn2(1)
       dn4(2)=t1*dn1(2)-t2*dn2(2)
       dn4(3)=t1*dn1(3)-t2*dn2(3)
       dn4(l2)=dn4(l2)+p2s
c
c dn4b (dij || s) order 2
       dn4b(1)=t1*dn2(1)-t2*dn2b(1)
       dn4b(2)=t1*dn2(2)-t2*dn2b(2)
       dn4b(3)=t1*dn2(3)-t2*dn2b(3)
       dn4b(l2)=dn4b(l2)+p3s
c dn6 and dn7 used for (pj | s) order 0 and 1
c
      dn6(1)=t1*x0x(1)-t2*x1x(1)
      dn6(2)=t1*x0x(2)-t2*x1x(2)
      dn6(3)=t1*x0x(3)-t2*x1x(3)
      dn6(l2)=dn6(l2)+s1s
c
      dn7(1)=t1*x1x(1)-t2*x2x(1)
      dn7(2)=t1*x1x(2)-t2*x2x(2)
      dn7(3)=t1*x1x(3)-t2*x2x(3)
      dn7(l2)=dn7(l2)+s2s
c
      dn7b(1)=t1*x2x(1)-t2*x3x(1)
      dn7b(2)=t1*x2x(2)-t2*x3x(2)
      dn7b(3)=t1*x2x(3)-t2*x3x(3)
      dn7b(l2)=dn7b(l2)+s3s
c
       t84=(dn6(1)-dn7(1))/z2
       t85=(dn6(2)-dn7(2))/z2
       t86=(dn6(3)-dn7(3))/z2
c
       t84b=(dn7(1)-dn7b(1))/z2
       t85b=(dn7(2)-dn7b(2))/z2
       t86b=(dn7(3)-dn7b(3))/z2
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s-s1s)/z2
        d1s=d1s+(s1s-s2s)/z2
        d2s=d2s+(s2s-s3s)/z2
        dn3(1)=dn3(1)+t26
        dn3(2)=dn3(2)+t27
        dn3(3)=dn3(3)+t28
        dn4(1)=dn4(1)+t29
        dn4(2)=dn4(2)+t30
        dn4(3)=dn4(3)+t31
        dn4b(1)=dn4b(1)+t32
        dn4b(2)=dn4b(2)+t33
        dn4b(3)=dn4b(3)+t34
       endif
c
       t96=(dn3(1)-dn4(1))/z2
       t97=(dn3(2)-dn4(2))/z2
       t98=(dn3(3)-dn4(3))/z2
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 706 l3=1,lij
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-xi(l3)
c
       d0p=t1*d0s-t2*d1s
       d1p=t1*d1s-t2*d2s
c
       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s
c
c dn8 and dn8b (pi||pk) ,   dn9 and dn9b (pj||pk)
       dn8(1)=t1*dn(1)-t2*dn1(1)
       dn8(2)=t1*dn(2)-t2*dn1(2)
       dn8(3)=t1*dn(3)-t2*dn1(3)
       dn8(l3)=dn8(l3)+p1s
       dn8b(1)=t1*dn1(1)-t2*dn2(1)
       dn8b(2)=t1*dn1(2)-t2*dn2(2)
       dn8b(3)=t1*dn1(3)-t2*dn2(3)
       dn8b(l3)=dn8b(l3)+p2s
c
       dn9(1)=t1*dn6(1)-t2*dn7(1)
       dn9(2)=t1*dn6(2)-t2*dn7(2)
       dn9(3)=t1*dn6(3)-t2*dn7(3)
       dn9(l3)=dn9(l3)+pj1s
       dn9b(1)=t1*dn7(1)-t2*dn7b(1)
       dn9b(2)=t1*dn7(2)-t2*dn7b(2)
       dn9b(3)=t1*dn7(3)-t2*dn7b(3)
       dn9b(l3)=dn9b(l3)+pj2s
c
c dn5 (dij || pk) dn5b (dij ||pk) order 1
       dn5(1)=t1*dn3(1)-t2*dn4(1)
       dn5(2)=t1*dn3(2)-t2*dn4(2)
       dn5(3)=t1*dn3(3)-t2*dn4(3)
       dn5(l3)=dn5(l3)+d1s
c
       dn5b(1)=t1*dn4(1)-t2*dn4b(1)
       dn5b(2)=t1*dn4(2)-t2*dn4b(2)
       dn5b(3)=t1*dn4(3)-t2*dn4b(3)
       dn5b(l3)=dn5b(l3)+d2s
c
       if (l1.eq.l3) then
        d0p=d0p+(pj0s-pj1s)/z2
        d1p=d1p+(pj1s-pj2s)/z2
        pi0p=pi0p+(s0s-s1s)/z2
        pi1p=pi1p+(s1s-s2s)/z2
       dn5(1)=dn5(1)+t84
       dn5(2)=dn5(2)+t85
       dn5(3)=dn5(3)+t86
       dn5b(1)=dn5b(1)+t84b
       dn5b(2)=dn5b(2)+t85b
       dn5b(3)=dn5b(3)+t86b
       dn8(1)=dn8(1)+t26
       dn8(2)=dn8(2)+t27
       dn8(3)=dn8(3)+t28
       dn8b(1)=dn8b(1)+t29
       dn8b(2)=dn8b(2)+t30
       dn8b(3)=dn8b(3)+t31
c
       endif
c
       if (l2.eq.l3) then
        d0p=d0p+(p0s-p1s)/z2
        d1p=d1p+(p1s-p2s)/z2
        pj0p=pj0p+(s0s-s1s)/z2
        pj1p=pj1p+(s1s-s2s)/z2
       dn5(1)=dn5(1)+t81
       dn5(2)=dn5(2)+t82
       dn5(3)=dn5(3)+t83
       dn5b(1)=dn5b(1)+t81b
       dn5b(2)=dn5b(2)+t82b
       dn5b(3)=dn5b(3)+t83b
       dn9(1)=dn9(1)+t26
       dn9(2)=dn9(2)+t27
       dn9(3)=dn9(3)+t28
       dn9b(1)=dn9b(1)+t29
       dn9b(2)=dn9b(2)+t30
       dn9b(3)=dn9b(3)+t31
c
       endif
c
        t90=(dn9(1)-dn9b(1))/z2
        t91=(dn9(2)-dn9b(2))/z2
        t92=(dn9(3)-dn9b(3))/z2
        t93=(dn8(1)-dn8b(1))/z2
        t94=(dn8(2)-dn8b(2))/z2
        t95=(dn8(3)-dn8b(3))/z2
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
       do 706 l4=1,lk
c
       f2=1.D0
       t1=Q(l4)-R(Nuc(j),l4)
       t2=Q(l4)-xi(l4)
       tna=t1*d0p-t2*d1p
c
c dn10 : (dij || dkl) nuclear derivative needed
c
       dn10(1)=t1*dn5(1)-t2*dn5b(1)
       dn10(2)=t1*dn5(2)-t2*dn5b(2)
       dn10(3)=t1*dn5(3)-t2*dn5b(3)
       dn10(l4)=dn10(l4)+d1p
c
       if (l4.eq.l1) then
        tna=tna+(pj0p-pj1p)/z2
        dn10(1)=dn10(1)+t90
        dn10(2)=dn10(2)+t91
        dn10(3)=dn10(3)+t92
c
       endif
c
       if (l4.eq.l2) then
        tna=tna+(pi0p-pi1p)/z2
        dn10(1)=dn10(1)+t93
        dn10(2)=dn10(2)+t94
        dn10(3)=dn10(3)+t95

       endif
c
       if (l4.eq.l3) then
        f2=sq3
        tna=tna+(d0s-d1s)/z2
        dn10(1)=dn10(1)+t96
        dn10(2)=dn10(2)+t97
        dn10(3)=dn10(3)+t98
       endif
c
       cc=ccoef/(f1*f2)
       term=cc*tna
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       kk=ii+((M2-jj)*(jj-1))/2
       te=cc*RMM(kk)
        V=V+term*RMM(kk)
c
        Ef(1)=Ef(1)+dn10(1)*te
        Ef(2)=Ef(2)+dn10(2)*te
        Ef(3)=Ef(3)+dn10(3)*te
c
c
 706  continue
 701  continue
 700  continue
c
c     write(*,*) 'Electrostatic potential =',V
c     write(*,990) xi(1),xi(2),xi(3)
c     write(*,*) Ef(1),Ef(2),Ef(3)
c
 990  format ('At point ',3(F13.4))
      return
      end
c-------------------------------------------------------------------


 
