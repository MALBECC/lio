c-------------------------------------------------------------------
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
c-------------------------------------------------------------------
      module subm_intsolG; contains
      subroutine intsolG(f,ff)

      use liotemp   , only: FUNCT
      use garcha_mod, only: RMM, ll, a, c, d, r, pc, Iz, Nuc, Ncont
     >                    , nshell, watermod, rmax, pi, pi32, NORM
     >                    , natom, ntatom, M, Md, nsol
c
      implicit none
      real*8  :: f(natom,3), fcampo(natom,3), ff(ntatom,3), rr(3)

c auxiliar quantities
      real*8  :: dn(3),  dn1(3), dn2(3), dn3(3), dn4(3), dn5(3)
      real*8  :: dn6(3), dn7(3), dn8(3), dn9(3), dn10(3)
      real*8  :: dn2b(3), dn4b(3), dn5b(3), dn7b(3), dn8b(3)
      real*8  :: dn9b(3), dn11(3),dn12(3)
      real*8  :: Q(3)

      logical, dimension(:),   allocatable :: qmmmdo
      real*8,  dimension(:),   allocatable :: s0s, s1s, s2s, s3s
      real*8,  dimension(:),   allocatable :: s4s, s5s, s6s
      real*8,  dimension(:,:), allocatable :: x0x, x1x, x2x
      real*8,  dimension(:,:), allocatable :: x3x, x4x, x5x

! Implicits:
      integer :: i, j, k, n, ni, nj, ii, jj, j1, j2
      integer :: ncerca, nvecin, ns, np, nd
      integer :: l, lk, lij, l1, l2, l3, l4, l5, l12, l34
      integer :: MM, MMd
      integer :: M1, M2, M3, M5, M7, M9, M11

      real*8  :: f1, f2, te, et, q1, q2, q3, rexp, sq3
      real*8  :: term, temp, temp0, zij, z2, u
      real*8  :: distint, distx, disty, distz
      real*8  :: tx, ty, tz, tx1, ty1, tz1, ti, tj, tna, tn1a
      real*8  :: ss, t1, dd2, dd, alf, cc, ccoef
      real*8  :: dNs, dNp, dNd, dNf, dN1s
      real*8  :: fNs, fNp, fNd, piNs, sNpi
      real*8  :: pNp, pNd, pN1p
      real*8  :: s0p, s1p, s2p
      real*8  :: p0s, p1s, p2s, p3s, p4s
      real*8  :: pi0p, pi0d, pi1p, pi1d, pi2p
      real*8  :: pj0s, pj0p, pj0d, pj1s, pj1p, pj1d, pj2s, pj2p, pj3s
      real*8  :: d0s, d0p, d1s, d1p, d2s, d3s, d2p, d0pl, d1pl
      real*8  :: t2, t3, t4, t5, t7, t8, t9, t15
      real*8  :: t25, t26, t27, t28, t29, t30, t31, t32, t33, t34
      real*8  :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59
      real*8  :: t60, t61, t62, t63, t64, t65, t66, t67, t68, t69
      real*8  :: t70, t71, t72, t73, t74, t81, t81b, t82, t82b
      real*8  :: t83, t83b, t84, t84b, t85, t85b, t86, t86b
      real*8  :: t90, t91, t92, t93, t94, t95, t96, t97, t98

c
c -----------------------------
c
c
c       write(*,*) rmax,rmax
      allocate (qmmmdo(natom+nsol))
c      allocate(d(natom,natom))
      allocate(s0s(ntatom),s1s(ntatom),s2s(ntatom),s3s(ntatom)
     > ,s4s(ntatom),s5s(ntatom),s6s(ntatom))
 
      allocate(x0x(ntatom,3),x1x(ntatom,3),x2x(ntatom,3),x3x(ntatom,3)
     >,x4x(ntatom,3),x5x(ntatom,3))


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

c SOLVENT-SOLUTE ELECTROSTATIC AND LJ  INTERACTIONS
c
c       write(*,*) 'qmmmcut=',qmmmcut
c       qmmmcut2=qmmmcut**2
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
       do 125 j1=1,natom
       do 125 j2=natom+1,natom+nsol
c
c
       tx1=r(j1,1)-r(j2,1)
       ty1=r(j1,2)-r(j2,2)
       tz1=r(j1,3)-r(j2,3)
       dd2=tx1**2+ty1**2+tz1**2
        if(qmmmdo(j2)) then
       dd2=sqrt(dd2)
c
c
c        Ese=Ese+Iz(j1)*pc(ni)/dd2
c
c
       et = -Iz(j1)*pc(j2)/dd2**3
c       et2= -pc(j2)/dd2**3

c        fcampo(j1,1) = fcampo(j1,1) + et2*tx1
c        fcampo(j1,2) = fcampo(j1,2) + et2*ty1
c        fcampo(j1,3) = fcampo(j1,3) + et2*tz1
c       write(99,*) Iz(j1),pc(j2),dd2
c
       f(j1,1)=f(j1,1) + et*tx1
       f(j1,2)=f(j1,2) + et*ty1
       f(j1,3)=f(j1,3) + et*tz1
c
       ff(j2,1)=ff(j2,1) - et*tx1
       ff(j2,2)=ff(j2,2) - et*ty1
       ff(j2,3)=ff(j2,3) - et*tz1
c

       endif
 125   continue
c        do j1=1, natom
c       write(6,*) fcampo(j1,1:3)
c        enddo

c      do 50 i=1,natom
c      do 50 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 50   continue

c      write(33,*) natom+nsol
c      do i=1,natom+nsol
c          if (i.lt.natom.or.abs(pc(i)).gt.1.D-5) then
c       write(33,*) 'C',r(i,1),r(i,2),r(i,3)
c         endif
c      enddo

c      write(34,*) natom+nsol
c      do i=1,natom+nsol
c       write(34,*) 'C',pc(I)
c      enddo


c
c first ioop (s|s) case -------------------------------------------
c
       call g2g_timer_start('ss')

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
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      ccoef=c(i,ni)*c(j,nj)
      rexp=alf*dd
c
      if (rexp.lt.rmax) then

      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
c
      k=i+((M2-j)*(j-1))/2
c 
c loop over nuclei, nuclear attraction matrix elements
c tna: accumulates nuc. attraction over all nuclei
c
       tna=0.D0
       temp0=2.D0*sqrt(zij/pi)*ss
c
      do 202 n=natom+1,natom+nsol
       if (qmmmdo(n)) then
c
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
202    continue
 
c
c l2: different p in the p shell GRADIENT PART ----------
c
      te=RMM(k)*ccoef
      ty=te*2.D0
      t5=ty*a(j,nj)
      t4=ty*a(i,ni)
      
c
      do 205 l2=1,3
        t1=Q(l2)-r(Nuc(i),l2)
        tx=r(Nuc(i),l2)-r(Nuc(j),l2)
c
c loop over classical nuclei
      do 203 n=natom+1,natom+nsol
       if(qmmmdo(n)) then
c
       piNs=t1*s0s(n)-(Q(l2)-r(n,l2))*s1s(n)
       sNpi=piNs+tx*s0s(n)
       f(Nuc(i),l2)=f(Nuc(i),l2)+t4*piNs
       f(Nuc(j),l2)=f(Nuc(j),l2)+t5*sNpi
c
        ff(n,l2)=ff(n,l2)+te*x0x(n,l2)
        endif
 203  continue

 205    continue
        endif
c --------------------------------------------------------
 200  continue
       call g2g_timer_stop('ss')

       call g2g_timer_start('ps')
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
      rexp=alf*dd
       if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
c
c loop over nuclei, part common for all shell
      temp0=2.D0*sqrt(zij/pi)*ss
c
      do 302 n=natom+1,natom+nsol
       if (qmmmdo(n)) then
c
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
 302  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c        te=RMM(k)*ccoef
c        ty=te*2.D0
c        t5=ty*a(j,nj)
c        t4=ty*a(i,ni)
c
cc nuclear attraction part
c
      do 303 n=natom+1,natom+nsol
       if(qmmmdo(n)) then
c
c
      t50=(s0s(n)-s1s(n))/z2
c
      do 306 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
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
        f(Nuc(i),l2)=f(Nuc(i),l2)-te*s0s(n)
       endif
c
      tx=r(Nuc(i),l2)-r(Nuc(j),l2)
      pNp=dNs+tx*p0s
c
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
      f(Nuc(i),l2)=f(Nuc(i),l2)+t4*dNs
      f(Nuc(j),l2)=f(Nuc(j),l2)+t5*pNp
c
         ff(n,l2)=ff(n,l2)+te*dn(l2)
c
 308  continue
 306  continue
       endif
c
 303  continue
c end nuclear attr. part ----------
       endif
 300  continue
       call g2g_timer_stop('ps')
       call g2g_timer_start('pp')
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
       rexp=alf*dd
       if (rexp.lt.rmax) then 
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 402 n=natom+1,natom+nsol
c
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
 402  continue
c
      ccoef=c(i,ni)*c(j,nj)
c Nuclear attraction part ----------
      do 403 n=natom+1,natom+nsol
       if(qmmmdo(n)) then

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
        f(Nuc(i),l3)=f(Nuc(i),l3)-te*s0p
       endif
c
       if (l2.eq.l3) then
        dNp=dNp+t30
       f(Nuc(j),l3)=f(Nuc(j),l3)-te*p0s
       endif
c
       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
       pNd=dNp+tx*pNp
c
        f(Nuc(i),l3)=f(Nuc(i),l3)+t4*dNp
        f(Nuc(j),l3)=f(Nuc(j),l3)+t5*pNd
c
         ff(n,l3)=ff(n,l3)+te*dn2(l3)

 406  continue
c
       endif
 403  continue
      endif
c
 400  continue
       call g2g_timer_stop('pp')
       call g2g_timer_start('ds')
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
      rexp=alf*dd
      if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss
c
c
c loop over nuclei, part common for all shell
      do 502 n=natom+1,natom+nsol
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
 502  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
cc nuclear attraction part 
c
      do 503 n=natom+1,natom+nsol
       if(qmmmdo(n)) then
c
      t7=(s0s(n)-s1s(n))/z2
      t8=(s1s(n)-s2s(n))/z2
      t26=(x0x(n,1)-x1x(n,1))/z2
      t27=(x0x(n,2)-x1x(n,2))/z2
      t28=(x0x(n,3)-x1x(n,3))/z2
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
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
        f(Nuc(i),l3)=f(Nuc(i),l3)-te*pj0s
       endif
c
       if (l2.eq.l3) then
        dNp=dNp+t30
       f(Nuc(i),l3)=f(Nuc(i),l3)-te*p0s
       endif
c
        fNs=dNp-tx*dNs
c
        f(Nuc(i),l3)=f(Nuc(i),l3)+t4*fNs
        f(Nuc(j),l3)=f(Nuc(j),l3)+t5*dNp
c
         ff(n,l3)=ff(n,l3)+te*dn2(l3)
 506  continue
       endif
 503  continue
c end nuclear attr. part ----------
       endif
 500  continue
       call g2g_timer_stop('ds')
       call g2g_timer_start('dp')
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
      rexp=alf*dd
      if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 602 n=natom+1,natom+nsol
       if (qmmmdo(n)) then
c

       q1=Q(1)-r(n,1)
       q2=Q(2)-r(n,2)
       q3=Q(3)-r(n,3)
       u=q1**2+q2**2+q3**2
       
       
       u=u*zij
       temp=-temp0*pc(n)
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
       endif
 602  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c Nuclear attraction part ----------
      do 603 n=natom+1,natom+nsol
       if(qmmmdo(n)) then
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
         f(Nuc(i),l4)=f(Nuc(i),l4)-te*pj0p
         dNd=dNd+(pj0p-pj1p)/z2
        endif
c
        if(l2.eq.l4) then
         f(Nuc(i),l4)=f(Nuc(i),l4)-te*pi0p
         dNd=dNd+(pi0p-pi1p)/z2
        endif
c
        if(l3.eq.l4) then
         f(Nuc(j),l4)=f(Nuc(j),l4)-te*d0s
         dNd=dNd+(d0s-d1s)/z2
        endif
c
        fNp=dNd-tx*tna
        f(Nuc(i),l4)=f(Nuc(i),l4)+t4*fNp
        f(Nuc(j),l4)=f(Nuc(j),l4)+t5*dNd
c
         ff(n,l4)=ff(n,l4)+te*dn5(l4)
c
        fNp=dNd-tx*tna
 606  continue
      endif
c
 603  continue
      endif
c
 600  continue
       call g2g_timer_stop('dp')
       call g2g_timer_start('dd')
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
      rexp=alf*dd
      if (rexp.lt.rmax) then
      ss=pi32*exp(-rexp)/(zij*sqrt(zij))
      temp0=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 702 n=natom+1,natom+nsol
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
 702  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c Loop over nuclei - Nuclear attraction part ---
      do 703 n=natom+1,natom+nsol
       if(qmmmdo(n)) then
c
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
       t2=Q(l4)-r(n,l4)
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
        f(Nuc(i),l5)=f(Nuc(i),l5)-te*pj0d
       endif
c
       if (l2.eq.l5) then
        fNd=fNd+t73
        f(Nuc(i),l5)=f(Nuc(i),l5)-te*pi0d
       endif
c
       if (l3.eq.l5) then
        fNd=fNd+t74
        f(Nuc(j),l5)=f(Nuc(j),l5)-te*d0pl
       endif
c
       if (l4.eq.l5) then
        fNd=fNd+t64
        f(Nuc(j),l5)=f(Nuc(j),l5)-te*d0p
       endif
c
        dNf=fNd+tx*tna
        f(Nuc(i),l5)=f(Nuc(i),l5)+t4*fNd
        f(Nuc(j),l5)=f(Nuc(j),l5)+t5*dNf
c
         ff(n,l5)=ff(n,l5)+te*dn10(l5)
 706  continue
       endif
c
 703  continue
c end nuclear attraction part --------
       endif
c
 700  continue
       call g2g_timer_stop('dd')
c       do 128 j1=1,natom
c128     write(79,*) f(j1,1),f(j1,2),f(j1,3)
c       do j1=natom+1,nsol+natom
c         write(*,*) 'nsol',natom,nsol
c       write(78,*) ff(j1,1),ff(j1,2),ff(j1,3)
c        enddo

c      stop
c
c----------------------------------------------------------------
    
      deallocate (qmmmdo)
      deallocate(s0s,s1s,s2s,s3s
     > ,s4s,s5s,s6s,x0x,x1x,x2x,x3x
     >,x4x,x5x)


      return
      end subroutine
      end module subm_intsolG
c-------------------------------------------------------------------


 
