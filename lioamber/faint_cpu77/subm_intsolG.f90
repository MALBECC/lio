
! GRADIENT VERSION
! 1 e integrals
! using the Obara-Saika recursive method.
! debugged ( or supposed to) 28-7-92
! Dario Estrin

      module subm_intsolG; contains
      subroutine intsolG(f,ff)

      use liotemp   , only: FUNCT
      use garcha_mod, only: RMM, ll, a, c, d, r, pc, Iz, Nuc, Ncont, &
                           nshell, watermod, rmax, pi, pi32, NORM, &
                           natom, ntatom, M, Md, nsol

      implicit none
      double precision  :: f(natom,3), fcampo(natom,3), ff(ntatom,3), rr(3)

      double precision  :: dn(3),  dn1(3), dn2(3), dn3(3), dn4(3), dn5(3)
      double precision  :: dn6(3), dn7(3), dn8(3), dn9(3), dn10(3)
      double precision  :: dn2b(3), dn4b(3), dn5b(3), dn7b(3), dn8b(3)
      double precision  :: dn9b(3), dn11(3),dn12(3)
      double precision  :: Q(3)

      logical, dimension(:),   allocatable :: qmmmdo
      double precision,  dimension(:),   allocatable :: s0s, s1s, s2s, s3s
      double precision,  dimension(:),   allocatable :: s4s, s5s, s6s
      double precision,  dimension(:,:), allocatable :: x0x, x1x, x2x
      double precision,  dimension(:,:), allocatable :: x3x, x4x, x5x

! Implicits:
      integer :: i, j, k, n, ni, nj, ii, jj, j1, j2
      integer :: ncerca, nvecin, ns, np, nd
      integer :: l, lk, lij, l1, l2, l3, l4, l5, l12, l34
      integer :: MM, MMd
      integer :: M1, M2, M3, M5, M7, M9, M11

      double precision  :: f1, f2, te, et, q1, q2, q3, rexp, sq3
      double precision  :: term, temp, temp0, zij, z2, u
      double precision  :: distint, distx, disty, distz
      double precision  :: tx, ty, tz, tx1, ty1, tz1, ti, tj, tna, tn1a
      double precision  :: ss, t1, dd2, dd, alf, cc, ccoef
      double precision  :: dNs, dNp, dNd, dNf, dN1s
      double precision  :: fNs, fNp, fNd, piNs, sNpi
      double precision  :: pNp, pNd, pN1p
      double precision  :: s0p, s1p, s2p
      double precision  :: p0s, p1s, p2s, p3s, p4s
      double precision  :: pi0p, pi0d, pi1p, pi1d, pi2p
      double precision  :: pj0s, pj0p, pj0d, pj1s, pj1p, pj1d, pj2s, pj2p, pj3s
      double precision  :: d0s, d0p, d1s, d1p, d2s, d3s, d2p, d0pl, d1pl
      double precision  :: t2, t3, t4, t5, t7, t8, t9, t15
      double precision  :: t25, t26, t27, t28, t29, t30, t31, t32, t33, t34
      double precision  :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59
      double precision  :: t60, t61, t62, t63, t64, t65, t66, t67, t68, t69
      double precision  :: t70, t71, t72, t73, t74, t81, t81b, t82, t82b
      double precision  :: t83, t83b, t84, t84b, t85, t85b, t86, t86b
      double precision  :: t90, t91, t92, t93, t94, t95, t96, t97, t98

      allocate (qmmmdo(natom+nsol))
      allocate(s0s(ntatom),s1s(ntatom),s2s(ntatom),s3s(ntatom), &
      s4s(ntatom),s5s(ntatom),s6s(ntatom))

      allocate(x0x(ntatom,3),x1x(ntatom,3),x2x(ntatom,3),x3x(ntatom,3), &
      x4x(ntatom,3),x5x(ntatom,3))


      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif

      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      M2=2*M

      do l=1,3
      Ll(l)=l*(l-1)/2
     enddo

       qmmmdo=.true.
        if(watermod.eq.1) nvecin=2 !tre sites
        if(watermod.eq.2) nvecin=3 !four sites
        if(watermod.eq.1.or.watermod.eq.2) then
        qmmmdo=.false.

       do i=natom+1,natom+nsol
          ncerca=0
          if(pc(i).gt.1.-6D0) qmmmdo(i)=.true. ! eliminate the sites without charge

        do j=natom+1,natom+nsol


          distx=(r(i,1)-r(j,1))**2
          disty=(r(i,2)-r(j,2))**2
          distz=(r(i,3)-r(j,3))**2
          distint=distx+disty+distz
          if (distint.lt.8.45D0) then
            ncerca=ncerca+1
          endif

         enddo
          if(ncerca.le.nvecin) qmmmdo(i)=.false. ! eliminate the broken waters
         enddo
          endif

       fcampo=0
       do j1=1,natom
       do j2=natom+1,natom+nsol

       tx1=r(j1,1)-r(j2,1)
       ty1=r(j1,2)-r(j2,2)
       tz1=r(j1,3)-r(j2,3)
       dd2=tx1**2+ty1**2+tz1**2
        if(qmmmdo(j2)) then
       dd2=sqrt(dd2)
       et = -Iz(j1)*pc(j2)/dd2**3
       f(j1,1)=f(j1,1) + et*tx1
       f(j1,2)=f(j1,2) + et*ty1
       f(j1,3)=f(j1,3) + et*tz1

       ff(j2,1)=ff(j2,1) - et*tx1
       ff(j2,2)=ff(j2,2) - et*ty1
       ff(j2,3)=ff(j2,3) - et*tz1

       endif
    enddo
enddo

! first ioop (s|s) case
      do i=1,ns
      do j=1,i
     dd=d(Nuc(i),Nuc(j))

      do ni=1,ncont(i)
      do nj=1,ncont(j)
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      ccoef=c(i,ni)*c(j,nj)
      rexp=alf*dd

      if (rexp.lt.rmax) then

      ss=pi32*exp(-rexp)/(zij*sqrt(zij))

      k=i+((M2-j)*(j-1))/2
       tna=0.D0
       temp0=2.D0*sqrt(zij/pi)*ss

      do n=natom+1,natom+nsol
       if (qmmmdo(n)) then

       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=(q1**2+q2**2+q3**2)


       u=u*zij
       temp=-pc(n)*temp0
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       temp=z2*s1s(n)
       x0x(n,1)=temp*q1
       x0x(n,2)=temp*q2
       x0x(n,3)=temp*q3
       tna=tna+s0s(n)
       endif
    enddo
      te=RMM(k)*ccoef
      ty=te*2.D0
      t5=ty*a(j,nj)
      t4=ty*a(i,ni)

      do l2=1,3
        t1=Q(l2)-r(Nuc(i),l2)
        tx=r(Nuc(i),l2)-r(Nuc(j),l2)

      do n=natom+1,natom+nsol
       if(qmmmdo(n)) then

       piNs=t1*s0s(n)-(Q(l2)-r(n,l2))*s1s(n)
       sNpi=piNs+tx*s0s(n)
       f(Nuc(i),l2)=f(Nuc(i),l2)+t4*piNs
       f(Nuc(j),l2)=f(Nuc(j),l2)+t5*sNpi

        ff(n,l2)=ff(n,l2)+te*x0x(n,l2)
        endif
enddo
enddo
        endif
     enddo
     enddo
  enddo
  enddo

! (p|s) case  and gradients

      do i=ns+1,ns+np,3
      do j=1,ns

      dd=d(Nuc(i),Nuc(j))

      do ni=1,ncont(i)
      do nj=1,ncont(j)

      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      alf=ti*a(j,nj)
      rexp=alf*dd
       if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss

      do n=natom+1,natom+nsol
       if (qmmmdo(n)) then

       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2



       u=u*zij
       temp=-temp0*pc(n)
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
       endif
 enddo

      ccoef=c(i,ni)*c(j,nj)
      do n=natom+1,natom+nsol
       if(qmmmdo(n)) then
      t50=(s0s(n)-s1s(n))/z2

      do l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)

       dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
       dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
       dn(3)=t1*x0x(n,3)-t2*x1x(n,3)

       dn(l1)=dn(l1)+s1s(n)

        ii=i+l1-1
        k=ii+((M2-j)*(j-1))/2
        te=RMM(k)*ccoef
      do l2=1,3
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       dNs=t1*p0s-t2*p1s
       if (l1.eq.l2) then
        dNs=dNs+t50
        f(Nuc(i),l2)=f(Nuc(i),l2)-te*s0s(n)
       endif

      tx=r(Nuc(i),l2)-r(Nuc(j),l2)
      pNp=dNs+tx*p0s

        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
      f(Nuc(i),l2)=f(Nuc(i),l2)+t4*dNs
      f(Nuc(j),l2)=f(Nuc(j),l2)+t5*pNp

         ff(n,l2)=ff(n,l2)+te*dn(l2)
      enddo
      enddo
       endif
    enddo
           endif

    enddo
    enddo
enddo
enddo

! (p|p) case and gradients
      do i=ns+1,ns+np,3
      do j=ns+1,i,3

      dd=d(Nuc(i),Nuc(j))

      do ni=1,ncont(i)
      do nj=1,ncont(j)

      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      alf=ti*a(j,nj)
       rexp=alf*dd
       if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss

      do n=natom+1,natom+nsol

       if (qmmmdo(n)) then

       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2


       u=u*zij
       temp=-temp0*pc(n)
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

        endif
     enddo

      ccoef=c(i,ni)*c(j,nj)
      do n=natom+1,natom+nsol
       if(qmmmdo(n)) then

      t15=(s0s(n)-s1s(n))/z2
      t25=(s1s(n)-s2s(n))/z2
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2

      do l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       t30=(p0s-p1s)/z2
       p2s=t1*s2s(n)-t2*s3s(n)

       !dn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)

      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)

      lij=3
      if (i.eq.j) then
       lij=l1
      endif

      do l2=1,lij
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       t3=Q(l2)-r(Nuc(j),l2)
       s0p=t3*s0s(n)-t2*s1s(n)
       s1p=t3*s1s(n)-t2*s2s(n)
       t29=(s0p-s1p)/z2
       pNp=t3*p0s-t2*p1s
       pN1p=t3*p1s-t2*p2s

       dn2(1)=t3*dn(1)-t2*dn1(1)
       dn2(2)=t3*dn(2)-t2*dn1(2)
       dn2(3)=t3*dn(3)-t2*dn1(3)
       dn2(l2)=dn2(l2)+p1s

       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s

       if (l1.eq.l2) then
        pNp=pNp+t15
       pN1p=pN1p+t25
        dn2(1)=dn2(1)+t26
        dn2(2)=dn2(2)+t27
        dn2(3)=dn2(3)+t28
       endif

      ii=i+l1-1
      jj=j+l2-1
      k=ii+((M2-jj)*(jj-1))/2

        te=RMM(k)*ccoef
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)


      do l3=1,3
       t1=Q(l3)-r(Nuc(i),l3)
       t2=Q(l3)-r(n,l3)
       dNp=t1*pNp-t2*pN1p

       if (l1.eq.l3) then
        dNp=dNp+t29
        f(Nuc(i),l3)=f(Nuc(i),l3)-te*s0p
       endif

       if (l2.eq.l3) then
        dNp=dNp+t30
       f(Nuc(j),l3)=f(Nuc(j),l3)-te*p0s
       endif

       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
       pNd=dNp+tx*pNp

        f(Nuc(i),l3)=f(Nuc(i),l3)+t4*dNp
        f(Nuc(j),l3)=f(Nuc(j),l3)+t5*pNd

         ff(n,l3)=ff(n,l3)+te*dn2(l3)

      enddo
   enddo
   enddo
       endif
    enddo
      endif
   enddo
   enddo
enddo
enddo

   ! (d|s) case and gradients

      do i=ns+np+1,M,6
      do j=1,ns

      dd=d(Nuc(i),Nuc(j))

      do ni=1,ncont(i)
      do nj=1,ncont(j)

      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss

      do n=natom+1,natom+nsol
       if (qmmmdo(n)) then

       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2


       u=u*zij
       temp=-temp0*pc(n)
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
       endif
    enddo

      ccoef=c(i,ni)*c(j,nj)

      do n=natom+1,natom+nsol
       if(qmmmdo(n)) then

      t7=(s0s(n)-s1s(n))/z2
      t8=(s1s(n)-s2s(n))/z2
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2

      do l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       p2s=t1*s2s(n)-t2*s3s(n)
       t30=(p0s-p1s)/z2
       ! dn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)

      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)

      do l2=1,l1

       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       tna=t1*p0s-t2*p1s
       tn1a=t1*p1s-t2*p2s

       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
       t29=(pj0s-pj1s)/z2

       dn2(1)=t1*dn(1)-t2*dn1(1)
       dn2(2)=t1*dn(2)-t2*dn1(2)
       dn2(3)=t1*dn(3)-t2*dn1(3)
       dn2(l2)=dn2(l2)+p1s

       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+t7
        tn1a=tn1a+t8
        f1=sq3
        dn2(1)=dn2(1)+t26
        dn2(2)=dn2(2)+t27
        dn2(3)=dn2(3)+t28
       endif

       dNs=tna
       dN1s=tn1a

       l12=l1*(l1-1)/2+l2
       ii=i+l12-1

       k=ii+((M2-j)*(j-1))/2
       cc=ccoef/f1
       term=cc*dNs
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)

        do l3=1,3
        t1=Q(l3)-r(Nuc(j),l3)
        t2=Q(l3)-r(n,l3)
        tx=r(Nuc(i),l3)-r(Nuc(j),l3)

        dNp=t1*dNs-t2*dN1s

       if (l1.eq.l3) then
        dNp=dNp+t29
        f(Nuc(i),l3)=f(Nuc(i),l3)-te*pj0s
       endif

       if (l2.eq.l3) then
        dNp=dNp+t30
       f(Nuc(i),l3)=f(Nuc(i),l3)-te*p0s
       endif

        fNs=dNp-tx*dNs

        f(Nuc(i),l3)=f(Nuc(i),l3)+t4*fNs
        f(Nuc(j),l3)=f(Nuc(j),l3)+t5*dNp

         ff(n,l3)=ff(n,l3)+te*dn2(l3)
      enddo
      enddo
   enddo
       endif
    enddo

       endif

    enddo
    enddo
    enddo
    enddo

! (d|p) case
      do i=ns+np+1,M,6
      do j=ns+1,ns+np,3

      dd=d(Nuc(i),Nuc(j))

      do ni=1,ncont(i)
      do nj=1,ncont(j)

      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss

      do n=natom+1,natom+nsol
       if (qmmmdo(n)) then
       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2


       u=u*zij
       temp=-temp0*pc(n)

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
       endif
    enddo
      ccoef=c(i,ni)*c(j,nj)
      do n=natom+1,natom+nsol
       if(qmmmdo(n)) then

      t7=(s0s(n)-s1s(n))/z2
      t8=(s1s(n)-s2s(n))/z2
      t9=(s2s(n)-s3s(n))/z2
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2
      t29=(x1x(n,1)-x2x(n,1))/z2
      t30=(x1x(n,2)-x2x(n,2))/z2
      t31=(x1x(n,3)-x2x(n,3))/z2

      do l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       p2s=t1*s2s(n)-t2*s3s(n)
       p3s=t1*s3s(n)-t2*s4s(n)

   ! dn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)

      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)

       t51=(dn(1)-dn1(1))/z2
       t52=(dn(2)-dn1(2))/z2
       t53=(dn(3)-dn1(3))/z2

      dn2(1)=t1*x2x(n,1)-t2*x3x(n,1)
      dn2(2)=t1*x2x(n,2)-t2*x3x(n,2)
      dn2(3)=t1*x2x(n,3)-t2*x3x(n,3)
      dn2(l1)=dn2(l1)+s3s(n)

      do l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
       pj2s=t1*s2s(n)-t2*s3s(n)

       ! dn3 (dij || s) order 0
       dn3(1)=t1*dn(1)-t2*dn1(1)
       dn3(2)=t1*dn(2)-t2*dn1(2)
       dn3(3)=t1*dn(3)-t2*dn1(3)
       dn3(l2)=dn3(l2)+p1s

       ! dn4 (dij || s) order 1
       dn4(1)=t1*dn1(1)-t2*dn2(1)
       dn4(2)=t1*dn1(2)-t2*dn2(2)
       dn4(3)=t1*dn1(3)-t2*dn2(3)
       dn4(l2)=dn4(l2)+p2s

       ! dn6 and dn7 used for (pj | s) order 0 and 1
      dn6(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn6(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn6(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn6(l2)=dn6(l2)+s1s(n)

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

      do l3=1,3
         ! dn5 (dij || Pk ) order 0, the one needed for derivatives
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(n,l3)
       tna=t1*d0s-t2*d1s
       tn1a=t1*d1s-t2*d2s

       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s

       dn5(1)=t1*dn3(1)-t2*dn4(1)
       dn5(2)=t1*dn3(2)-t2*dn4(2)
       dn5(3)=t1*dn3(3)-t2*dn4(3)
       dn5(l3)=dn5(l3)+d1s

       if (l1.eq.l3) then
        tna=tna+(pj0s-pj1s)/z2
        tn1a=tn1a+(pj1s-pj2s)/z2
       pi0p=pi0p+t7
       pi1p=pi1p+t8
       dn5(1)=dn5(1)+t54
       dn5(2)=dn5(2)+t55
       dn5(3)=dn5(3)+t56
       endif

       if (l2.eq.l3) then
        tna=tna+(p0s-p1s)/z2
        tn1a=tn1a+(p1s-p2s)/z2
       pj0p=pj0p+t7
       pj1p=pj1p+t8
       dn5(1)=dn5(1)+t51
       dn5(2)=dn5(2)+t52
       dn5(3)=dn5(3)+t53
       endif

       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1

      k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
       term=cc*tna

        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
      do l4=1,3
        t1=Q(l4)-r(Nuc(j),l4)
        t2=Q(l4)-r(n,l4)
        tx=r(Nuc(i),l4)-r(Nuc(j),l4)

        dNd=t1*tna-t2*tn1a

       if(l1.eq.l4) then
         f(Nuc(i),l4)=f(Nuc(i),l4)-te*pj0p
         dNd=dNd+(pj0p-pj1p)/z2
        endif

        if(l2.eq.l4) then
         f(Nuc(i),l4)=f(Nuc(i),l4)-te*pi0p
         dNd=dNd+(pi0p-pi1p)/z2
        endif

        if(l3.eq.l4) then
         f(Nuc(j),l4)=f(Nuc(j),l4)-te*d0s
         dNd=dNd+(d0s-d1s)/z2
        endif

        fNp=dNd-tx*tna
        f(Nuc(i),l4)=f(Nuc(i),l4)+t4*fNp
        f(Nuc(j),l4)=f(Nuc(j),l4)+t5*dNd

         ff(n,l4)=ff(n,l4)+te*dn5(l4)

        fNp=dNd-tx*tna
     enddo
     enddo
  enddo
  enddo
      endif
   enddo
      endif
   enddo
   enddo
enddo
enddo

! (d|d) case
      do i=ns+np+1,M,6
      do j=ns+np+1,i,6

      dd=d(Nuc(i),Nuc(j))

      do ni=1,ncont(i)
      do nj=1,ncont(j)

      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      rexp=alf*dd
      if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss

      ! loop over nuclei, part common for all shell
      do  n=natom+1,natom+nsol
       if (qmmmdo(n)) then

       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2

       u=u*zij
       temp=-temp0*pc(n)
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
       endif
    enddo

      ccoef=c(i,ni)*c(j,nj)
      do n=natom+1,natom+nsol
       if(qmmmdo(n)) then
      t50=(s0s(n)-s1s(n))/z2
      t51=(s1s(n)-s2s(n))/z2
      t52=(s2s(n)-s3s(n))/z2
      t53=(s3s(n)-s4s(n))/z2

      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2
      t29=(x1x(n,1)-x2x(n,1))/z2
      t30=(x1x(n,2)-x2x(n,2))/z2
      t31=(x1x(n,3)-x2x(n,3))/z2
      t32=(x2x(n,1)-x3x(n,1))/z2
      t33=(x2x(n,2)-x3x(n,2))/z2
      t34=(x2x(n,3)-x3x(n,3))/z2

      do l1=1,3
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
       ! dn(u) (pi|Au|s)
      dn(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn(l1)=dn(l1)+s1s(n)

      dn1(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn1(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn1(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn1(l1)=dn1(l1)+s2s(n)

       t81=(dn(1)-dn1(1))/z2
       t82=(dn(2)-dn1(2))/z2
       t83=(dn(3)-dn1(3))/z2

      dn2(1)=t1*x2x(n,1)-t2*x3x(n,1)
      dn2(2)=t1*x2x(n,2)-t2*x3x(n,2)
      dn2(3)=t1*x2x(n,3)-t2*x3x(n,3)
      dn2(l1)=dn2(l1)+s3s(n)

      dn2b(1)=t1*x3x(n,1)-t2*x4x(n,1)
      dn2b(2)=t1*x3x(n,2)-t2*x4x(n,2)
      dn2b(3)=t1*x3x(n,3)-t2*x4x(n,3)
      dn2b(l1)=dn2b(l1)+s4s(n)

       t81b=(dn1(1)-dn2(1))/z2
       t82b=(dn1(2)-dn2(2))/z2
       t83b=(dn1(3)-dn2(3))/z2

      do l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
       pj2s=t1*s2s(n)-t2*s3s(n)
       pj3s=t1*s3s(n)-t2*s4s(n)

       ! dn3 (dij || s) order 0
       dn3(1)=t1*dn(1)-t2*dn1(1)
       dn3(2)=t1*dn(2)-t2*dn1(2)
       dn3(3)=t1*dn(3)-t2*dn1(3)
       dn3(l2)=dn3(l2)+p1s
       ! dn4 (dij || s) order 1
       dn4(1)=t1*dn1(1)-t2*dn2(1)
       dn4(2)=t1*dn1(2)-t2*dn2(2)
       dn4(3)=t1*dn1(3)-t2*dn2(3)
       dn4(l2)=dn4(l2)+p2s
       ! dn4b (dij || s) order 2
       dn4b(1)=t1*dn2(1)-t2*dn2b(1)
       dn4b(2)=t1*dn2(2)-t2*dn2b(2)
       dn4b(3)=t1*dn2(3)-t2*dn2b(3)
       dn4b(l2)=dn4b(l2)+p3s
       ! dn6 and dn7 used for (pj | s) order 0 and 1

      dn6(1)=t1*x0x(n,1)-t2*x1x(n,1)
      dn6(2)=t1*x0x(n,2)-t2*x1x(n,2)
      dn6(3)=t1*x0x(n,3)-t2*x1x(n,3)
      dn6(l2)=dn6(l2)+s1s(n)

      dn7(1)=t1*x1x(n,1)-t2*x2x(n,1)
      dn7(2)=t1*x1x(n,2)-t2*x2x(n,2)
      dn7(3)=t1*x1x(n,3)-t2*x2x(n,3)
      dn7(l2)=dn7(l2)+s2s(n)

      dn7b(1)=t1*x2x(n,1)-t2*x3x(n,1)
      dn7b(2)=t1*x2x(n,2)-t2*x3x(n,2)
      dn7b(3)=t1*x2x(n,3)-t2*x3x(n,3)
      dn7b(l2)=dn7b(l2)+s3s(n)

       t84=(dn6(1)-dn7(1))/z2
       t85=(dn6(2)-dn7(2))/z2
       t86=(dn6(3)-dn7(3))/z2

       t84b=(dn7(1)-dn7b(1))/z2
       t85b=(dn7(2)-dn7b(2))/z2
       t86b=(dn7(3)-dn7b(3))/z2

       t58=(pj0s-pj1s)/z2
       t59=(pj1s-pj2s)/z2
       t60=(pj2s-pj3s)/z2

       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
       d3s=t1*p3s-t2*p4s

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

       t61=(d0s-d1s)/z2
       t62=(d1s-d2s)/z2
       t63=(d2s-d3s)/z2

       t96=(dn3(1)-dn4(1))/z2
       t97=(dn3(2)-dn4(2))/z2
       t98=(dn3(3)-dn4(3))/z2

      lij=3
      if (i.eq.j) then
       lij=l1
      endif

      do l3=1,lij

       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(n,l3)

       s0p=t1*s0s(n)-t2*s1s(n)
       s1p=t1*s1s(n)-t2*s2s(n)
       s2p=t1*s2s(n)-t2*s3s(n)
       t70=(s0p-s1p)/z2
       t71=(s1p-s2p)/z2

       d0p=t1*d0s-t2*d1s
       d1p=t1*d1s-t2*d2s
       d2p=t1*d2s-t2*d3s

       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pi2p=t1*p2s-t2*p3s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s
       pj2p=t1*pj2s-t2*pj3s

       ! dn8 and dn8b (pi||pk) ,   dn9 and dn9b (pj||pk)
       dn8(1)=t1*dn(1)-t2*dn1(1)
       dn8(2)=t1*dn(2)-t2*dn1(2)
       dn8(3)=t1*dn(3)-t2*dn1(3)
       dn8(l3)=dn8(l3)+p1s
       dn8b(1)=t1*dn1(1)-t2*dn2(1)
       dn8b(2)=t1*dn1(2)-t2*dn2(2)
       dn8b(3)=t1*dn1(3)-t2*dn2(3)
       dn8b(l3)=dn8b(l3)+p2s

       dn9(1)=t1*dn6(1)-t2*dn7(1)
       dn9(2)=t1*dn6(2)-t2*dn7(2)
       dn9(3)=t1*dn6(3)-t2*dn7(3)
       dn9(l3)=dn9(l3)+pj1s
       dn9b(1)=t1*dn7(1)-t2*dn7b(1)
       dn9b(2)=t1*dn7(2)-t2*dn7b(2)
       dn9b(3)=t1*dn7(3)-t2*dn7b(3)
       dn9b(l3)=dn9b(l3)+pj2s

       ! dn5 (dij || pk) dn5b (dij ||pk) order 1
       dn5(1)=t1*dn3(1)-t2*dn4(1)
       dn5(2)=t1*dn3(2)-t2*dn4(2)
       dn5(3)=t1*dn3(3)-t2*dn4(3)
       dn5(l3)=dn5(l3)+d1s

       dn5b(1)=t1*dn4(1)-t2*dn4b(1)
       dn5b(2)=t1*dn4(2)-t2*dn4b(2)
       dn5b(3)=t1*dn4(3)-t2*dn4b(3)
       dn5b(l3)=dn5b(l3)+d2s

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

        t64=(d0p-d1p)/z2
        t65=(d1p-d2p)/z2
        t66=(pi0p-pi1p)/z2
        t67=(pi1p-pi2p)/z2
        t68=(pj0p-pj1p)/z2
        t69=(pj1p-pj2p)/z2

        t90=(dn9(1)-dn9b(1))/z2
        t91=(dn9(2)-dn9b(2))/z2
        t92=(dn9(3)-dn9b(3))/z2
        t93=(dn8(1)-dn8b(1))/z2
        t94=(dn8(2)-dn8b(2))/z2
        t95=(dn8(3)-dn8b(3))/z2

      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
       do l4=1,lk

       f2=1.D0
       t1=Q(l4)-R(Nuc(j),l4)
       t2=Q(l4)-r(n,l4)
       tna=t1*d0p-t2*d1p
       tn1a=t1*d1p-t2*d2p

       ! dn10 : (dij || dkl) nuclear derivative needed
       dn10(1)=t1*dn5(1)-t2*dn5b(1)
       dn10(2)=t1*dn5(2)-t2*dn5b(2)
       dn10(3)=t1*dn5(3)-t2*dn5b(3)
       dn10(l4)=dn10(l4)+d1p

       d0pl=t1*d0s-t2*d1s
       d1pl=t1*d1s-t2*d2s
       pj0d=t1*pj0p-t2*pj1p
       pj1d=t1*pj1p-t2*pj2p
       pi0d=t1*pi0p-t2*pi1p
       pi1d=t1*pi1p-t2*pi2p

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

       t72=(pj0d-pj1d)/z2
       t73=(pi0d-pi1d)/z2
       t74=(d0pl-d1pl)/z2
       cc=ccoef/(f1*f2)
       term=cc*tna

       l12=Ll(l1)+l2
       l34=Ll(l3)+l4
       ii=i+l12-1
       jj=j+l34-1

       k=ii+((M2-jj)*(jj-1))/2

! gradients
        te=RMM(k)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
      do l5=1,3
       t1=Q(l5)-r(Nuc(i),l5)
       t2=Q(l5)-r(n,l5)
       tx=r(Nuc(i),l5)-r(Nuc(j),l5)

       fNd=t1*tna-t2*tn1a

       if (l1.eq.l5) then
        fNd=fNd+t72
        f(Nuc(i),l5)=f(Nuc(i),l5)-te*pj0d
       endif

       if (l2.eq.l5) then
        fNd=fNd+t73
        f(Nuc(i),l5)=f(Nuc(i),l5)-te*pi0d
       endif

       if (l3.eq.l5) then
        fNd=fNd+t74
        f(Nuc(j),l5)=f(Nuc(j),l5)-te*d0pl
       endif

       if (l4.eq.l5) then
        fNd=fNd+t64
        f(Nuc(j),l5)=f(Nuc(j),l5)-te*d0p
       endif

        dNf=fNd+tx*tna
        f(Nuc(i),l5)=f(Nuc(i),l5)+t4*fNd
        f(Nuc(j),l5)=f(Nuc(j),l5)+t5*dNf

         ff(n,l5)=ff(n,l5)+te*dn10(l5)
      enddo
      enddo
   enddo
      enddo
      enddo
       endif
    enddo
       endif
    enddo
    enddo
enddo
enddo
      deallocate (qmmmdo)
      deallocate(s0s,s1s,s2s,s3s, &
      s4s,s5s,s6s,x0x,x1x,x2x,x3x, &
      x4x,x5x)


      return
      end subroutine
      end module subm_intsolG
