!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module subm_int1G; contains
       subroutine int1G(ff)
!------------------------------------------------------------------------------!
c GRADIENT VERSION
c calculates 1 e part and all gradients, to be used with MD.f
c
c Integrals subroutine 
c 1 e integrals
c using the Obara-Saika recursive method.
c
c 
c loop over all basis functions
c now the basis is supposed to be ordered according to the type,
c all s, then all p, then all d, .....
c inside each type, are ordered in shells
c
c ns ... marker for end of s
c np ... marker for end of p
c nd ... marker for end of d
c
c r(Nuc(i),j) j component of position of nucleus i , j=1,3
c Input : basis function information
c Output: F matrix, and S matrix and forces on nuclei
c all gradients, up to d functions
c debugged ( or supposed to) 28-7-92
c Dario Estrin
!------------------------------------------------------------------------------!
       use liotemp   , only: FUNCT
!       use garcha_mod
       use garcha_mod, only: RMM, Nuc, a, c, d, r, Iz, ncont
     >                     , nshell, pi, pi32, NORM, natom, M, Md 
     >                     , ll, ntq
       implicit none

c-----auxiliar quantities
       real*8, intent(inout) :: ff(natom,3)

       integer :: natomold, igpu
       integer :: n, i, j, k, ii, jj, ni, nj
       integer :: l1, l2, l3, l4, l12, l34
       integer :: MM, MMd, ns, np, nd
       integer :: M1, M2, M3, M5, M7, M9, M11

       real*8  :: En, ovlap, alf, alf2, alf3, alf4
       real*8  :: Q(3), term, temp, sq3, cc, ccoef
       real*8  :: f1, f2, tn, tna, u, z2, zij
       real*8  :: ss, ps, dd, p0s, p1s, p2s, p3s
       real*8  :: pi0p, pi1p, piks, pikpk, pipk, pis
       real*8  :: pj0s, pj1s, pj2s, pj0p, pj1p, pjkpk
       real*8  :: pjks, pjpk, pjs, pks, sks
       real*8  :: dijs, dijpk, dijks, dijkpk
       real*8  :: d0s, d0p, d1p, d1s, d2s
       real*8  :: t0, t1, t2

       integer :: l, lij, lk, l5
       real*8  :: temp0, pp, pd, ds, dp, df, fs, fp, fd, q1, q2, q3
       real*8  :: p4s, d0pl, d1pl, d2p, d3s, dijkpl, dijpl
       real*8  :: dkd, dkf, dkp, dks, dsd, fkd, fkp, fks
       real*8  :: dN1s, dNs, dNp, dNd, dNf, fNd, fNp
       real*8  :: pNd, pN1p, pkp, pkd, pjkdkl, pjdkl, pj3s, pj2p
       real*8  :: pj1d, pj0d, piNs, pikdkl, pidkl, pi2p, pi1d, pi0d
       real*8  :: spk, spj, sNpi, skpk, skpj, skpi, s2p, s1p, s0p
       real*8  :: pNp, fNs
       real*8  :: tt, tx, ty, te, ti, tj, tn1a
       real*8  :: t3, t4, t5, t7, t8, t9
       real*8  :: t10, t11, t12, t13, t14, t15, t16, t17, t18, t19
       real*8  :: t20, t21, t22, t23, t24, t25, t26, t27, t28, t29
       real*8  :: t30, t31, t32, t33, t34, t35, t36, t37, t38, t39
       real*8  :: t40, t41
       real*8  :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59
       real*8  :: t60, t61, t62, t63, t64, t65, t66, t67, t68, t69
       real*8  :: t70, t71, t72, t73, t74
       real*8  :: t81, t82, t83, t84, t85, t86
       real*8  :: t81b, t82b, t83b, t84b, t85b, t86b
       real*8  :: t90, t91, t92, t93, t94, t95, t96, t97, t98

       real*8  :: s0s(ntq), s1s(ntq), s2s(ntq), s3s(ntq)
       real*8  :: s4s(ntq), s5s(ntq), s6s(ntq)
       real*8  :: x0x(ntq,3), x1x(ntq,3), x2x(ntq,3)
       real*8  :: x3x(ntq,3), x4x(ntq,3)
       real*8  :: dn1(3), dn2(3), dn3(3), dn4(3), dn5(3)
       real*8  :: dn6(3), dn7(3), dn8(3), dn9(3), dn(3)
       real*8  :: dn10(3)
       real*8  :: dn2b(3), dn4b(3), dn5b(3), dn7b(3), dn8b(3), dn9b(3)
c distance between pairs of centers
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c
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
      do 181 l=1,3
 181   Ll(l)=l*(l-1)/2
c
c Pointers
c first P
      M1=1
c now  S
      M3=M1+MM
c now F also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c---- Overlap ,Kinetic energy and Nuclear Attraction
c      matrix elements evaluation
c Overlap matrix will be kept, kinetic energy and nuclear attraction
c matrix elements no,
c they're stored in Fock matrix and in the Energy directly
c in order to reduce the memory requirements
c
c
c      do 50 i=1,natom
c      do 50 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 50   continue
c
c Nuclear Repulsion part ------------------------------------------
      En=0.D0
      do 51 i=1,natom
      do 51 j=1,i-1
 51    En=En+Iz(i)*Iz(j)/sqrt(d(i,j))
c
      do 52 i=1,natom
c
       ff(i,1)=0.D0
       ff(i,2)=0.D0
       ff(i,3)=0.D0
      do 53 j=1,i-1
       tt=Iz(i)*Iz(j)/d(i,j)**1.5D0

       ff(i,1)=ff(i,1)-tt*(r(i,1)-r(j,1))
       ff(i,2)=ff(i,2)-tt*(r(i,2)-r(j,2))
 53    ff(i,3)=ff(i,3)-tt*(r(i,3)-r(j,3))
c
      do 54 j=i+1,natom
       tt=Iz(i)*Iz(j)/d(i,j)**1.5D0

       ff(i,1)=ff(i,1)-tt*(r(i,1)-r(j,1))
       ff(i,2)=ff(i,2)-tt*(r(i,2)-r(j,2))
 54    ff(i,3)=ff(i,3)-tt*(r(i,3)-r(j,3))
c
 52   continue

      call aint_query_gpu_level(igpu)
      ! doing nuclear attraction part on GPU - KE part still is
      ! done here
      if (igpu.gt.3) then
        natomold = natom
        natom = 0
      endif
c
c first loop (s|s) case -------------------------------------------
c
      do 200 i=1,ns
      do 200 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
c
      do 200 ni=1,ncont(i)
      do 200 nj=1,ncont(j)
c
c (0|0) calculation
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=alf*2.D0
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      ccoef=c(i,ni)*c(j,nj)
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ovlap=ss
      sks=alf*(3.D0-alf2*dd)*ovlap
      tn=sks
c
      k=i+((M2-j)*(j-1))/2
c 
c loop over nuclei, nuclear attraction matrix elements
c tna: accumulates nuc. attraction over all nuclei
c
       tna=0.D0
       temp0=2.D0*sqrt(zij/pi)*ss
c
      do 202 n=1,natom
       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2
       u=u*zij
       temp=-Iz(n)*temp0
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       temp=z2*s1s(n)
       x0x(n,1)=temp*q1
       x0x(n,2)=temp*q2
       x0x(n,3)=temp*q3
 202   tna=tna+s0s(n)
c
c l2: different p in the p shell GRADIENT PART ----------
c
      t3=alf2*ss
      te=RMM(k)*ccoef
      ty=te*2.D0
      t5=ty*a(j,nj)
      t4=ty*a(i,ni)
      
c
      do 205 l2=1,3
        t1=Q(l2)-r(Nuc(i),l2)
        pis=t1*ss
        piks=t1*sks+alf2*pis
        tx=r(Nuc(i),l2)-r(Nuc(j),l2)
        skpi=piks+ tx*(sks+t3)
c
        ff(Nuc(i),l2)=ff(Nuc(i),l2)+t4*piks
        ff(Nuc(j),l2)=ff(Nuc(j),l2)+t5*skpi
c
c loop over nuclei, specific part
      do 203 n=1,natom
       piNs=t1*s0s(n)-(Q(l2)-r(n,l2))*s1s(n)
       sNpi=piNs+tx*s0s(n)
       ff(Nuc(i),l2)=ff(Nuc(i),l2)+t4*piNs
       ff(Nuc(j),l2)=ff(Nuc(j),l2)+t5*sNpi
c
        ff(n,l2)=ff(n,l2)+te*x0x(n,l2)
 203  continue

 205    continue
c --------------------------------------------------------
 200  continue
c
c------------------------------------------------------------------
c (p|s) case  and gradients
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
      z2=2.D0*zij
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      alf=ti*a(j,nj)
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
c
c loop over nuclei, part common for all shell
      temp0=2.D0*sqrt(zij/pi)*ss
c
      do 302 n=1,natom
       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2
       u=u*zij
       temp=-temp0*Iz(n)
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
        temp=z2*s1s(n)
       x0x(n,1)=temp*q1
       x0x(n,2)=temp*q2
       x0x(n,3)=temp*q3
        temp=z2*s2s(n)
       x1x(n,1)=temp*q1
       x1x(n,2)=temp*q2
       x1x(n,3)=temp*q3
 302  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
       t10=ss/z2
       t15=sks/z2
       t20=t15-alf*ss/a(i,ni)
c
      do 305 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       ps=ss*t1
       pks=sks*t1+alf2*ps
c
        ii=i+l1-1
c ii index , taking into account different components of the shell
c
        k=ii+((M2-j)*(j-1))/2
c
        te=RMM(k)*ccoef
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
      do 307 l2=1,3
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(Nuc(j),l2)
       ds=t1*ps
       dks=t1*pks
       pp=t2*ps
       pkp=t2*pks
c
       if (l1.eq.l2) then
        ff(Nuc(i),l2)=ff(Nuc(i),l2)-te*sks
        ds=ds+t10
        dks=dks+t20
        pp=pp+t10
        pkp=pkp+t15
       endif
c
       dks=dks+alf2*ds
       pkp=pkp+alf2*pp
c
        ff(Nuc(i),l2)=ff(Nuc(i),l2)+t4*dks
        ff(Nuc(j),l2)=ff(Nuc(j),l2)+t5*pkp
c
 307  continue
 305  continue
cc nuclear attraction part
c
      do 303 n=1,natom
c
      t50=(s0s(n)-s1s(n))/z2
c
      do 306 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-R(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
c
       dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
       dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
       dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
c  
       dn(l1)=dn(l1)+s1s(n)
c
        ii=i+l1-1
c ii index , taking into account different components of the shell
c
        k=ii+((M2-j)*(j-1))/2
c
        te=RMM(k)*ccoef
c
      do 308 l2=1,3
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       dNs=t1*p0s-t2*p1s
c
       if (l1.eq.l2) then
        dNs=dNs+t50
        ff(Nuc(i),l2)=ff(Nuc(i),l2)-te*s0s(n)
       endif
c
      tx=r(Nuc(i),l2)-r(Nuc(j),l2)
      pNp=dNs+tx*p0s
c
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
      ff(Nuc(i),l2)=ff(Nuc(i),l2)+t4*dNs
      ff(Nuc(j),l2)=ff(Nuc(j),l2)+t5*pNp
        ff(n,l2)=ff(n,l2)+te*dn(l2)
 308  continue
 306  continue
c
 303  continue
c end nuclear attr. part ----------
 300  continue
c-----------------------------------------------------------------
c (p|p) case and gradients
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
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      alf=ti*a(j,nj)
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
      temp0=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 402 n=1,natom
       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2
       u=u*zij
       temp=-temp0*Iz(n)
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
       s3s(n)=temp*FUNCT(3,u)
        temp=z2*s1s(n)
       x0x(n,1)=temp*q1
       x0x(n,2)=temp*q2
       x0x(n,3)=temp*q3
        temp=z2*s2s(n)
       x1x(n,1)=temp*q1
       x1x(n,2)=temp*q2
       x1x(n,3)=temp*q3
        temp=z2*s3s(n)
       x2x(n,1)=temp*q1
       x2x(n,2)=temp*q2
       x2x(n,3)=temp*q3
 402  continue
c
c
      t10=ss/z2
      t20=sks/z2
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 405 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       piks=sks*t1+alf2*pis
       t11=pis/z2
       t12=piks/z2
       t16=t12-alf/a(j,nj)*pis
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 405 l2=1,lij
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(Nuc(j),l2)
       spj=t2*ss
       skpj=sks*t2+alf2*spj
       t13=spj/z2
       t15=skpj/z2
       t14=t15-alf*spj/a(i,ni)
c
       pp=t2*pis
       pkp=t2*piks
c
       if (l1.eq.l2) then
        pp=pp+t10
        pkp=pkp+t20
       endif
c
       pkp=pkp+alf2*pp
c
      ii=i+l1-1
      jj=j+l2-1
      k=ii+((M2-jj)*(jj-1))/2
c
        te=RMM(k)*ccoef
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
      do 405 l3=1,3
c
       t1=Q(l3)-r(Nuc(i),l3)
       t2=Q(l3)-r(Nuc(j),l3)
       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
       dp=t1*pp
       dkp=t1*pkp
       pkd=t2*pkp
c
       if (l1.eq.l3) then
        dp=dp+t13
        dkp=dkp+t14
        pkd=pkd+t15
        ff(Nuc(i),l3)=ff(Nuc(i),l3)-te*skpj
       endif
c
       if (l2.eq.l3) then
        dp=dp+t11
        dkp=dkp+t12
        pkd=pkd+t16
        ff(Nuc(j),l3)=ff(Nuc(j),l3)-te*piks
       endif
c
       pd=dp+tx*pp
       dkp=dkp+alf2*dp
       pkd=pkd+alf2*pd
c
        ff(Nuc(i),l3)=ff(Nuc(i),l3)+t4*dkp
        ff(Nuc(j),l3)=ff(Nuc(j),l3)+t5*pkd
c
 405  continue
c
c Nuclear attraction part ----------
      do 403 n=1,natom
c
      t15=(s0s(n)-s1s(n))/z2
      t25=(s1s(n)-s2s(n))/z2
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2
c
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       t30=(p0s-p1s)/z2
       p2s=t1*s2s(n)-t2*s3s(n)
c
c dn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)
c
      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 406 l2=1,lij
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       t3=Q(l2)-r(Nuc(j),l2)
       s0p=t3*s0s(n)-t2*s1s(n)
       s1p=t3*s1s(n)-t2*s2s(n)
       t29=(s0p-s1p)/z2
       pNp=t3*p0s-t2*p1s
       pN1p=t3*p1s-t2*p2s
c
       dn2(1)=t3*dn(1)-t2*dn1(1)
       dn2(2)=t3*dn(2)-t2*dn1(2)
       dn2(3)=t3*dn(3)-t2*dn1(3)
       dn2(l2)=dn2(l2)+p1s
c
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
c
       if (l1.eq.l2) then
        pNp=pNp+t15
        pN1p=pN1p+t25
        dn2(1)=dn2(1)+t26
        dn2(2)=dn2(2)+t27
        dn2(3)=dn2(3)+t28
       endif
c
c
      ii=i+l1-1
      jj=j+l2-1
      k=ii+((M2-jj)*(jj-1))/2
c
        te=RMM(k)*ccoef
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
c gradient part
      do 406 l3=1,3
       t1=Q(l3)-r(Nuc(i),l3)
       t2=Q(l3)-r(n,l3)
       dNp=t1*pNp-t2*pN1p
c
       if (l1.eq.l3) then
        dNp=dNp+t29
        ff(Nuc(i),l3)=ff(Nuc(i),l3)-te*s0p
       endif
c
       if (l2.eq.l3) then
        dNp=dNp+t30
       ff(Nuc(j),l3)=ff(Nuc(j),l3)-te*p0s
       endif
c
       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
       pNd=dNp+tx*pNp
c
        ff(Nuc(i),l3)=ff(Nuc(i),l3)+t4*dNp
        ff(Nuc(j),l3)=ff(Nuc(j),l3)+t5*pNd
        ff(n,l3)=ff(n,l3)+te*dn2(l3)
 406  continue
c
 403  continue
c
 400  continue
c
c-------------------------------------------------------------------
c (d|s) case and gradients
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
      alf3=alf/a(i,ni)
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
      temp0=2.D0*sqrt(zij/pi)*ss
c
      t10=ss/z2
      t11=sks/z2-alf3*ss
c
c loop over nuclei, part common for all shell
      do 502 n=1,natom
       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2
       u=u*zij
       temp=-temp0*Iz(n)
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
       s3s(n)=temp*FUNCT(3,u)
        temp=z2*s1s(n)
       x0x(n,1)=temp*q1
       x0x(n,2)=temp*q2
       x0x(n,3)=temp*q3
        temp=z2*s2s(n)
       x1x(n,1)=temp*q1
       x1x(n,2)=temp*q2
       x1x(n,3)=temp*q3
        temp=z2*s3s(n)
       x2x(n,1)=temp*q1
       x2x(n,2)=temp*q2
       x2x(n,3)=temp*q3
 502  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 505 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       ps=ss*t1
       pks=sks*t1+alf2*ps
       t12=ps/z2
       t13=pks/z2
       t17=t13-alf3*ps
c
      do 505 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       pjks=t1*sks+alf2*pjs
       t14=pjs/z2
       t15=pjks/z2
       t16=t15-alf3*pjs
       ovlap=t1*ps
       tn=t1*pks
c
       f1=1.D0
       if (l1.eq.l2) then
        ovlap=ovlap+t10
        tn=tn+t11
        f1=sq3
       endif
c
       tn=tn+alf2*ovlap
c
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       ii=i+l12-1
c
       cc=ccoef/f1
       term=cc*tn
       k=ii+((M2-j)*(j-1))/2
c now gradients
c
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
      do 505 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(Nuc(i),l3)
       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
       dp=t1*ovlap
       dkp=t1*tn
       fks=t2*tn
c
       if (l1.eq.l3) then
        dp=dp+t14
        dkp=dkp+t15
        fks=fks+t16
        ff(Nuc(i),l1)=ff(Nuc(i),l1)-te*pjks
       endif
c
       if (l2.eq.l3) then
        dp=dp+t12
        dkp=dkp+t13
        fks=fks+t17
        ff(Nuc(i),l2)=ff(Nuc(i),l2)-te*pks
       endif
c
       fs=dp-tx*ovlap
       dkp=dkp+alf2*dp
       fks=fks+alf2*fs
c
        ff(Nuc(i),l3)=ff(Nuc(i),l3)+t4*fks
        ff(Nuc(j),l3)=ff(Nuc(j),l3)+t5*dkp
c
 505  continue
c
cc nuclear attraction part 
c
      do 503 n=1,natom
c
      t7=(s0s(n)-s1s(n))/z2
      t8=(s1s(n)-s2s(n))/z2
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-R(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       p2s=t1*s2s(n)-t2*s3s(n)
       t30=(p0s-p1s)/z2
c
c dn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)
c
      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)
c
      do 506 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       tna=t1*p0s-t2*p1s 
       tn1a=t1*p1s-t2*p2s
c
       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
       t29=(pj0s-pj1s)/z2
c
       dn2(1)=t1*dn(1)-t2*dn1(1)
       dn2(2)=t1*dn(2)-t2*dn1(2)
       dn2(3)=t1*dn(3)-t2*dn1(3)
       dn2(l2)=dn2(l2)+p1s
c
       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+t7
        tn1a=tn1a+t8
        f1=sq3
        dn2(1)=dn2(1)+t26
        dn2(2)=dn2(2)+t27
        dn2(3)=dn2(3)+t28
       endif
c
       dNs=tna
       dN1s=tn1a
c
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       ii=i+l12-1
c
       k=ii+((M2-j)*(j-1))/2
       cc=ccoef/f1
       term=cc*dNs
c
c now gradients
c
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
        do 506 l3=1,3
        t1=Q(l3)-r(Nuc(j),l3)
        t2=Q(l3)-r(n,l3)
        tx=r(Nuc(i),l3)-r(Nuc(j),l3)
c
        dNp=t1*dNs-t2*dN1s
c
       if (l1.eq.l3) then
        dNp=dNp+t29
        ff(Nuc(i),l3)=ff(Nuc(i),l3)-te*pj0s
       endif
c
       if (l2.eq.l3) then
        dNp=dNp+t30
       ff(Nuc(i),l3)=ff(Nuc(i),l3)-te*p0s
       endif
c
        fNs=dNp-tx*dNs
c
        ff(Nuc(i),l3)=ff(Nuc(i),l3)+t4*fNs
        ff(Nuc(j),l3)=ff(Nuc(j),l3)+t5*dNp
        ff(n,l3)=ff(n,l3)+te*dn2(l3)
 506  continue
 503  continue
c end nuclear attr. part ----------
 500  continue
c-----------------------------------------------------------------
c
c (d|p) case 
c
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
      alf3=alf/a(i,ni)
      alf4=alf/a(j,nj)
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
      temp0=2.D0*sqrt(zij/pi)*ss
c
      t10=ss/z2
      t30=sks/z2
      t11=t30-alf3*ss
c loop over nuclei, part common for all shell
      do 602 n=1,natom
       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2
       u=u*zij
       temp=-temp0*Iz(n)
c
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
       s3s(n)=temp*FUNCT(3,u)
       s4s(n)=temp*FUNCT(4,u)
        temp=z2*s1s(n)
       x0x(n,1)=temp*q1
       x0x(n,2)=temp*q2
       x0x(n,3)=temp*q3
        temp=z2*s2s(n)
       x1x(n,1)=temp*q1
       x1x(n,2)=temp*q2
       x1x(n,3)=temp*q3
        temp=z2*s3s(n)
       x2x(n,1)=temp*q1
       x2x(n,2)=temp*q2
       x2x(n,3)=temp*q3
        temp=z2*s4s(n)
       x3x(n,1)=temp*q1
       x3x(n,2)=temp*q2
       x3x(n,3)=temp*q3
 602  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 605 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       piks=sks*t1+alf2*pis
c
       t14=pis/z2
       t15=piks/z2
      do 605 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       pjks=sks*t1+alf2*pjs
c
       t12=pjs/z2
       t13=pjks/z2
       dijs=t1*pis
       dijks=t1*piks
       f1=1.D0
c
       if (l1.eq.l2) then
        f1=sq3
        dijs=dijs+t10
        dijks=dijks+t11
       endif
c
       dijks=dijks+alf2*dijs
c
       t22=dijs/z2
       t23=dijks/z2
       t24=t23-alf4*dijs
c
      do 605 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       ovlap=t1*dijs 
       tn=t1*dijks
c
       pipk=t1*pis
       pikpk=t1*piks
       pjpk=t1*pjs
       pjkpk=t1*pjks
c
       if (l1.eq.l3) then
        pipk=pipk+t10
        pikpk=pikpk+t30
        ovlap=ovlap+t12
        tn=tn+t13
       endif
c
       if (l2.eq.l3) then
        pjpk=pjpk+t10
        pjkpk=pjkpk+t30
        ovlap=ovlap+t14
        tn=tn+t15
       endif
c
       pjkpk=pjkpk+alf2*pjpk
       pikpk=pikpk+alf2*pipk
c
       t16=pjpk/z2
       t17=pjkpk/z2
       t18=t17-alf3*pjpk
       t19=pipk/z2
       t20=pikpk/z2
       t21=t20-alf3*pipk
c
       tn=tn+alf2*ovlap
       
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
       cc=ccoef/f1
       term=cc*tn
       k=ii+((M2-jj)*(jj-1))/2
c gradients
c
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
      do 605 l4=1,3
       t1=Q(l4)-r(Nuc(j),l4)
       t2=Q(l4)-r(Nuc(i),l4)
       tx=r(Nuc(i),l4)-r(Nuc(j),l4)
c
       fkp=t2*tn
       dsd=t1*ovlap
       dkd=t1*tn
c
       if (l1.eq.l4) then
        dsd=dsd+t16
        dkd=dkd+t17
        fkp=fkp+t18
        ff(Nuc(i),l4)=ff(Nuc(i),l4)-te*pjkpk
       endif
c
       if (l2.eq.l4) then
        dsd=dsd+t19
        dkd=dkd+t20
        fkp=fkp+t21
        ff(Nuc(i),l4)=ff(Nuc(i),l4)-te*pikpk
       endif
c
       if (l3.eq.l4) then
        dsd=dsd+t22
        dkd=dkd+t24
        fkp=fkp+t23
        ff(Nuc(j),l4)=ff(Nuc(j),l4)-te*dijks
       endif
c
       fp=dsd-tx*ovlap
       dkd=dkd+alf2*dsd
       fkp=fkp+alf2*fp
c
        ff(Nuc(i),l4)=ff(Nuc(i),l4)+t4*fkp
        ff(Nuc(j),l4)=ff(Nuc(j),l4)+t5*dkd
c
 605  continue
c
c Nuclear attraction part ----------
      do 603 n=1,natom
c
      t7=(s0s(n)-s1s(n))/z2
      t8=(s1s(n)-s2s(n))/z2
      t9=(s2s(n)-s3s(n))/z2
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2
      t29=(x1x(n,1)-x2x(n,1))/z2
      t30=(x1x(n,2)-x2x(n,2))/z2
      t31=(x1x(n,3)-x2x(n,3))/z2
c
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       p2s=t1*s2s(n)-t2*s3s(n)
       p3s=t1*s3s(n)-t2*s4s(n)
c
c dn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)
c
      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)
c
       t51=(dn(1)-dn1(1))/z2
       t52=(dn(2)-dn1(2))/z2
       t53=(dn(3)-dn1(3))/z2
c
      dn2(1)=t1*x2x(n,1)-t2*x3x(n,1)
      dn2(2)=t1*x2x(n,2)-t2*x3x(n,2)
      dn2(3)=t1*x2x(n,3)-t2*x3x(n,3)
      dn2(l1)=dn2(l1)+s3s(n)
c
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
       pj2s=t1*s2s(n)-t2*s3s(n)
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
      dn6(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn6(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn6(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn6(l2)=dn6(l2)+s1s(n)
c
      dn7(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn7(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn7(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn7(l2)=dn7(l2)+s2s(n)

       t54=(dn6(1)-dn7(1))/z2
       t55=(dn6(2)-dn7(2))/z2
       t56=(dn6(3)-dn7(3))/z2
       
       f1=1.D0
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+t7
        d1s=d1s+t8
        d2s=d2s+t9
        dn3(1)=dn3(1)+t26
        dn3(2)=dn3(2)+t27
        dn3(3)=dn3(3)+t28
        dn4(1)=dn4(1)+t29
        dn4(2)=dn4(2)+t30
        dn4(3)=dn4(3)+t31
       endif
c
c
      do 606 l3=1,3
c
c dn5 (dij || Pk ) order 0, the one needed for derivatives
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(n,l3)
       tna=t1*d0s-t2*d1s
       tn1a=t1*d1s-t2*d2s
c
       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s
c
       dn5(1)=t1*dn3(1)-t2*dn4(1)
       dn5(2)=t1*dn3(2)-t2*dn4(2)
       dn5(3)=t1*dn3(3)-t2*dn4(3)
       dn5(l3)=dn5(l3)+d1s
c
       if (l1.eq.l3) then
        tna=tna+(pj0s-pj1s)/z2
        tn1a=tn1a+(pj1s-pj2s)/z2
       pi0p=pi0p+t7
       pi1p=pi1p+t8
       dn5(1)=dn5(1)+t54
       dn5(2)=dn5(2)+t55
       dn5(3)=dn5(3)+t56
       endif
c
       if (l2.eq.l3) then
        tna=tna+(p0s-p1s)/z2
        tn1a=tn1a+(p1s-p2s)/z2
       pj0p=pj0p+t7
       pj1p=pj1p+t8
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
      k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
       term=cc*tna
c
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
c gradient part
      do 606 l4=1,3
        t1=Q(l4)-r(Nuc(j),l4)
        t2=Q(l4)-r(n,l4)
        tx=r(Nuc(i),l4)-r(Nuc(j),l4)
c
        dNd=t1*tna-t2*tn1a
c
        if(l1.eq.l4) then
         ff(Nuc(i),l4)=ff(Nuc(i),l4)-te*pj0p
         dNd=dNd+(pj0p-pj1p)/z2
        endif
c
        if(l2.eq.l4) then
         ff(Nuc(i),l4)=ff(Nuc(i),l4)-te*pi0p
         dNd=dNd+(pi0p-pi1p)/z2
        endif
c
        if(l3.eq.l4) then
         ff(Nuc(j),l4)=ff(Nuc(j),l4)-te*d0s
         dNd=dNd+(d0s-d1s)/z2
        endif
c
        fNp=dNd-tx*tna
        ff(Nuc(i),l4)=ff(Nuc(i),l4)+t4*fNp
        ff(Nuc(j),l4)=ff(Nuc(j),l4)+t5*dNd
        ff(n,l4)=ff(n,l4)+te*dn5(l4)
c
        fNp=dNd-tx*tna
 606  continue
c
 603  continue
c
 600  continue
c
c-------------------------------------------------------------------
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
      alf3=alf/a(i,ni)
      alf4=alf/a(j,nj)
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
      temp0=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 702 n=1,natom
       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2
       u=u*zij
       temp=-temp0*Iz(n)
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
       s3s(n)=temp*FUNCT(3,u)
       s4s(n)=temp*FUNCT(4,u)
       s5s(n)=temp*FUNCT(5,u)
       s6s(n)=temp*FUNCT(6,u)
        temp=z2*s1s(n)
       x0x(n,1)=temp*q1
       x0x(n,2)=temp*q2
       x0x(n,3)=temp*q3
        temp=z2*s2s(n)
       x1x(n,1)=temp*q1
       x1x(n,2)=temp*q2
       x1x(n,3)=temp*q3
        temp=z2*s3s(n)
       x2x(n,1)=temp*q1
       x2x(n,2)=temp*q2
       x2x(n,3)=temp*q3
        temp=z2*s4s(n)
       x3x(n,1)=temp*q1
       x3x(n,2)=temp*q2
       x3x(n,3)=temp*q3
        temp=z2*s5s(n)
       x4x(n,1)=temp*q1
       x4x(n,2)=temp*q2
       x4x(n,3)=temp*q3
 702  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      t0=ss/z2
      t12=sks/z2
      t10=t12-alf3*ss
c
      do 705 l1=1,3
c
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       piks=sks*t1+alf2*pis
       t14=pis/z2
       t15=piks/z2
       t26=t15-alf4*pis
c
      do 705 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       pjks=sks*t1+alf2*pjs
       t11=pjs/z2
       t13=pjks/z2
       t25=t13-alf4*pjs
c
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       dijs=t1*pis
       dijks=t1*piks
c
       if (l1.eq.l2) then
        f1=sq3
        dijs=dijs+t0
        dijks=dijks+t10
       endif
c
       dijks=dijks+alf2*dijs
c
       t20=dijs/z2
       t21=dijks/z2-alf4*dijs
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c

       do 705 l3=1,lij
c
       t2=Q(l3)-r(Nuc(j),l3)
       spk=t2*ss
       skpk=t2*sks+alf2*spk
       t22=spk/z2
       t23=skpk/z2
c
       pipk=t2*pis
       pjpk=t2*pjs
       pikpk=t2*piks
       pjkpk=t2*pjks
       dijpk=t2*dijs
       dijkpk=t2*dijks
c
       if (l1.eq.l3) then
        pipk=pipk+t0
        dijpk=dijpk+t11
        pikpk=pikpk+t12
        dijkpk=dijkpk+t13
       endif
c
       if (l2.eq.l3) then
        pjpk=pjpk+t0
        dijpk=dijpk+t14
        pjkpk=pjkpk+t12
        dijkpk=dijkpk+t15
       endif
c
       pikpk=pikpk+alf2*pipk
       pjkpk=pjkpk+alf2*pjpk
       dijkpk=dijkpk+alf2*dijpk
c 
       t16=pjpk/z2
       t17=pjkpk/z2
       t18=pipk/z2
       t19=pikpk/z2
c
      t39=dijpk/z2
      t40=dijkpk/z2
      t41=t40-alf4*dijpk
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
       do 705 l4=1,lk
c
       f2=1.D0
       t1=Q(l4)-r(Nuc(j),l4)
       ovlap=t1*dijpk
       tn=t1*dijkpk
c
       pjdkl=t1*pjpk
       pjkdkl=t1*pjkpk
       pidkl=t1*pipk
       pikdkl=t1*pikpk
       dijpl=t1*dijs
       dijkpl=t1*dijks
c
       if (l1.eq.l4) then
        ovlap=ovlap+t16
        tn=tn+t17
        pidkl=pidkl+t22
        pikdkl=pikdkl+t23
        dijpl=dijpl+t11
        dijkpl=dijkpl+t13
       endif
c
       if (l2.eq.l4) then
        ovlap=ovlap+t18
        tn=tn+t19
        pjdkl=pjdkl+t22
        pjkdkl=pjkdkl+t23
        dijpl=dijpl+t14
        dijkpl=dijkpl+t15
       endif
c
       if (l3.eq.l4) then
        ovlap=ovlap+t20
        tn=tn+t21
        pjdkl=pjdkl+t11
        pjkdkl=pjkdkl+t25
        pidkl=pidkl+t14
        pikdkl=pikdkl+t26
        f2=sq3
       endif
c
       pikdkl=pikdkl+alf2*pidkl
       pjkdkl=pjkdkl+alf2*pjdkl
       dijkpl=dijkpl+alf2*dijpl
c
c l12 and l34 go from 1 to 6, spanning the d shell in
c the order xx, xy, yy, zx, zy, zz. The same order should be used
c in ordering the basis set, before calculating the matrix elements
c
       l12=Ll(l1)+l2
       l34=Ll(l3)+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
       tn= tn+alf2*ovlap
c
       cc=ccoef/(f1*f2)
       term=cc*tn
c
c
c gradients
c
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
      t30=pjdkl/z2
      t31=pjkdkl/z2
      t32=t31-alf3*pjdkl
      t33=pidkl/z2
      t34=pikdkl/z2
      t35=t34-alf3*pidkl
      t36=dijpl/z2
      t37=dijkpl/z2
      t38=t37-alf4*dijpl
c

      do 705 l5=1,3
      t1=Q(l5)-r(Nuc(i),l5)
      t2=Q(l5)-r(Nuc(j),l5)
      tx=r(Nuc(i),l5)-r(Nuc(j),l5)
c
      fd=t1*ovlap
      fkd=t1*tn
      dkf=t2*tn
c
      if (l1.eq.l5) then
       fd=fd+t30
       fkd=fkd+t32
       dkf=dkf+t31
       ff(Nuc(i),l5)=ff(Nuc(i),l5)-te*pjkdkl
      endif
c
      if (l2.eq.l5) then
       fd=fd+t33
       fkd=fkd+t35
       dkf=dkf+t34
       ff(Nuc(i),l5)=ff(Nuc(i),l5)-te*pikdkl
      endif
c
      if (l3.eq.l5) then
       fd=fd+t36
       fkd=fkd+t37
       dkf=dkf+t38
       ff(Nuc(j),l5)=ff(Nuc(j),l5)-te*dijkpl
      endif
c
      if (l4.eq.l5) then
       fd=fd+t39
       fkd=fkd+t40
       dkf=dkf+t41
       ff(Nuc(j),l5)=ff(Nuc(j),l5)-te*dijkpk
      endif
c
        df=fd+tx*ovlap
        fkd=fkd+alf2*fd
        dkf=dkf+alf2*df
c
        ff(Nuc(i),l5)=ff(Nuc(i),l5)+t4*fkd
        ff(Nuc(j),l5)=ff(Nuc(j),l5)+t5*dkf
c
 705  continue
c
c Loop over nuclei - Nuclear attraction part ---
      do 703 n=1,natom
c
      t50=(s0s(n)-s1s(n))/z2
      t51=(s1s(n)-s2s(n))/z2
      t52=(s2s(n)-s3s(n))/z2
      t53=(s3s(n)-s4s(n))/z2
c
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2
      t29=(x1x(n,1)-x2x(n,1))/z2
      t30=(x1x(n,2)-x2x(n,2))/z2
      t31=(x1x(n,3)-x2x(n,3))/z2
      t32=(x2x(n,1)-x3x(n,1))/z2
      t33=(x2x(n,2)-x3x(n,2))/z2
      t34=(x2x(n,3)-x3x(n,3))/z2
c
      do 706 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       p2s=t1*s2s(n)-t2*s3s(n)
       p3s=t1*s3s(n)-t2*s4s(n)
       p4s=t1*s4s(n)-t2*s5s(n)
       t54=(p0s-p1s)/z2
       t55=(p1s-p2s)/z2
       t56=(p2s-p3s)/z2
       t57=(p3s-p4s)/z2
c
cdn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)
c
      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)
c
       t81=(dn(1)-dn1(1))/z2
       t82=(dn(2)-dn1(2))/z2
       t83=(dn(3)-dn1(3))/z2
c
      dn2(1)=t1*x2x(n,1)-t2*x3x(n,1)
      dn2(2)=t1*x2x(n,2)-t2*x3x(n,2)
      dn2(3)=t1*x2x(n,3)-t2*x3x(n,3)
      dn2(l1)=dn2(l1)+s3s(n)
c
      dn2b(1)=t1*x3x(n,1)-t2*x4x(n,1)
      dn2b(2)=t1*x3x(n,2)-t2*x4x(n,2)
      dn2b(3)=t1*x3x(n,3)-t2*x4x(n,3)
      dn2b(l1)=dn2b(l1)+s4s(n)
c
       t81b=(dn1(1)-dn2(1))/z2
       t82b=(dn1(2)-dn2(2))/z2
       t83b=(dn1(3)-dn2(3))/z2

      do 706 l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
       pj2s=t1*s2s(n)-t2*s3s(n)
       pj3s=t1*s3s(n)-t2*s4s(n)
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
      dn6(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn6(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn6(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn6(l2)=dn6(l2)+s1s(n)
c
      dn7(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn7(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn7(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn7(l2)=dn7(l2)+s2s(n)
c
      dn7b(1)=t1*x2x(n,1)-t2*x3x(n,1)
      dn7b(2)=t1*x2x(n,2)-t2*x3x(n,2)
      dn7b(3)=t1*x2x(n,3)-t2*x3x(n,3)
      dn7b(l2)=dn7b(l2)+s3s(n)
c

       t84=(dn6(1)-dn7(1))/z2
       t85=(dn6(2)-dn7(2))/z2
       t86=(dn6(3)-dn7(3))/z2
c
       t84b=(dn7(1)-dn7b(1))/z2
       t85b=(dn7(2)-dn7b(2))/z2
       t86b=(dn7(3)-dn7b(3))/z2
c
       t58=(pj0s-pj1s)/z2
       t59=(pj1s-pj2s)/z2
       t60=(pj2s-pj3s)/z2
       
c
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
       d3s=t1*p3s-t2*p4s
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+t50
        d1s=d1s+t51
        d2s=d2s+t52
        d3s=d3s+t53
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
       t61=(d0s-d1s)/z2
       t62=(d1s-d2s)/z2
       t63=(d2s-d3s)/z2
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
       t2=Q(l3)-r(n,l3)
c
       s0p=t1*s0s(n)-t2*s1s(n)
       s1p=t1*s1s(n)-t2*s2s(n)
       s2p=t1*s2s(n)-t2*s3s(n)
       t70=(s0p-s1p)/z2
       t71=(s1p-s2p)/z2
c
       d0p=t1*d0s-t2*d1s
       d1p=t1*d1s-t2*d2s
       d2p=t1*d2s-t2*d3s
c
       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pi2p=t1*p2s-t2*p3s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s
       pj2p=t1*pj2s-t2*pj3s
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
        d0p=d0p+t58
        d1p=d1p+t59
        d2p=d2p+t60
        pi0p=pi0p+t50
        pi1p=pi1p+t51
        pi2p=pi2p+t52
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
       endif
c
       if (l2.eq.l3) then
        d0p=d0p+t54
        d1p=d1p+t55
        d2p=d2p+t56
        pj0p=pj0p+t50
        pj1p=pj1p+t51
        pj2p=pj2p+t52
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
       endif
c
        t64=(d0p-d1p)/z2
        t65=(d1p-d2p)/z2
        t66=(pi0p-pi1p)/z2
        t67=(pi1p-pi2p)/z2
        t68=(pj0p-pj1p)/z2
        t69=(pj1p-pj2p)/z2
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
c
       do 706 l4=1,lk
c
       f2=1.D0
       t1=Q(l4)-R(Nuc(j),l4)
       t2=Q(l4)-R(n,l4)
       tna=t1*d0p-t2*d1p
       tn1a=t1*d1p-t2*d2p
c
c dn10 : (dij || dkl) nuclear derivative needed
c
       dn10(1)=t1*dn5(1)-t2*dn5b(1)
       dn10(2)=t1*dn5(2)-t2*dn5b(2)
       dn10(3)=t1*dn5(3)-t2*dn5b(3)
       dn10(l4)=dn10(l4)+d1p
c
       d0pl=t1*d0s-t2*d1s
       d1pl=t1*d1s-t2*d2s
       pj0d=t1*pj0p-t2*pj1p
       pj1d=t1*pj1p-t2*pj2p
       pi0d=t1*pi0p-t2*pi1p
       pi1d=t1*pi1p-t2*pi2p
c
       if (l4.eq.l1) then
        tna=tna+t68
        tn1a=tn1a+t69
        d0pl=d0pl+t58
        d1pl=d1pl+t59
        pi0d=pi0d+t70
        pi1d=pi1d+t71
        dn10(1)=dn10(1)+t90
        dn10(2)=dn10(2)+t91
        dn10(3)=dn10(3)+t92
       endif
c
       if (l4.eq.l2) then
        tna=tna+t66
        tn1a=tn1a+t67
        d0pl=d0pl+t54
        d1pl=d1pl+t55
        pj0d=pj0d+t70
        pj1d=pj1d+t71
        dn10(1)=dn10(1)+t93
        dn10(2)=dn10(2)+t94
        dn10(3)=dn10(3)+t95
       endif
c
       if (l4.eq.l3) then
        f2=sq3
        tna=tna+t61
        tn1a=tn1a+t62
        pj0d=pj0d+t58
        pj1d=pj1d+t59
        pi0d=pi0d+t54
        pi1d=pi1d+t55
        dn10(1)=dn10(1)+t96
        dn10(2)=dn10(2)+t97
        dn10(3)=dn10(3)+t98
       endif
c
       t72=(pj0d-pj1d)/z2
       t73=(pi0d-pi1d)/z2
       t74=(d0pl-d1pl)/z2
       cc=ccoef/(f1*f2)
       term=cc*tna
c
       l12=Ll(l1)+l2
       l34=Ll(l3)+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
c
c gradients
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
      do 706 l5=1,3
       t1=Q(l5)-r(Nuc(i),l5)
       t2=Q(l5)-r(n,l5)
       tx=r(Nuc(i),l5)-r(Nuc(j),l5)
c
       fNd=t1*tna-t2*tn1a
c
       if (l1.eq.l5) then
        fNd=fNd+t72
        ff(Nuc(i),l5)=ff(Nuc(i),l5)-te*pj0d
       endif
c
       if (l2.eq.l5) then
        fNd=fNd+t73
        ff(Nuc(i),l5)=ff(Nuc(i),l5)-te*pi0d
       endif
c
       if (l3.eq.l5) then
        fNd=fNd+t74
        ff(Nuc(j),l5)=ff(Nuc(j),l5)-te*d0pl
       endif
c
       if (l4.eq.l5) then
        fNd=fNd+t64
        ff(Nuc(j),l5)=ff(Nuc(j),l5)-te*d0p
       endif
c
        dNf=fNd+tx*tna
        ff(Nuc(i),l5)=ff(Nuc(i),l5)+t4*fNd
        ff(Nuc(j),l5)=ff(Nuc(j),l5)+t5*dNf
        ff(n,l5)=ff(n,l5)+te*dn10(l5)
 706  continue
c
 703  continue
c end nuclear attraction part --------
c
 700  continue

      if (igpu.gt.3) natom = natomold
c
c test of penalty function ------------------------------------
c     f1=d(1,2)-2.89D0
c     f2=4.0D0*0.0D0*f1
c     ff(1,1)=ff(1,1)+f2*(r(1,1)-r(2,1))
c     ff(1,2)=ff(1,2)+f2*(r(1,2)-r(2,2))
c     ff(1,3)=ff(1,3)+f2*(r(1,3)-r(2,3))
c
c     ff(2,1)=ff(2,1)-f2*(r(1,1)-r(2,1))
c     ff(2,2)=ff(2,2)-f2*(r(1,2)-r(2,2))
c     ff(2,3)=ff(2,3)-f2*(r(1,3)-r(2,3))
c
c----------------------------------------------------------------
c     write(*,*) 'int1G'
c     do i=1,natom
c      write(*,*) i,ff(i,1),ff(i,2),ff(i,3)
c     enddo
      return
      end subroutine
      end module subm_int1G
c-------------------------------------------------------------------


 
