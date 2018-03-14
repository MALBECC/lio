c-------------------------------------------------------------------
c GRADIENT VERSION
c calculates gradients of overlap , to be used with
c geometry optimization subroutines
c
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
c Output, derivative part coming from overlap
c debugged ( or supposed to) 29-7-92
c Dario Estrin
c-------------------------------------------------------------------
      module subm_intSG; contains
      subroutine intSG(ff)
      use garcha_mod, only: RMM, ll, a, c, d, r, nuc, ncont, nshell
     >                    , pi32, natom, M, Md, NORM
c
      implicit none
      real*8, intent(inout) :: ff(natom,3)
      real*8                :: Q(3)

      integer :: i, j, k, ii, jj, ni, nj
      integer :: l, lk, lij, l1, l2, l3, l4, l5, l12, l34
      integer :: MM, MMd, ns, np, nd
      integer :: M1, M2, M3, M5, M7, M9, M11, M13, M15, M17

      real*8  :: ovlap, fsp, sq3, alf, cc, ccoef, con
      real*8  :: zij, z2, fs, fd, f1, f2
      real*8  :: ti, tj, tx, ty, te, t0, t1, t2, t4, t5
      real*8  :: t10, t11, t12, t13, t14, t15, t16, t17
      real*8  :: ss, spi, spj, spk
      real*8  :: ps, pp, pd, pidkl, pipk, pis, pjdkl, pjpk, pjs
      real*8  :: ds, dp, dd, df, dsd, dijpk, dijpl, dijs

! Implicits:

c distance between pairs of centers
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
      con=0.D0
      do 181 l=1,3
 181   Ll(l)=l*(l-1)/2
c
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
c W ( eigenvalues ), also space used in least-squares
      M13=M11+MM
c aux ( vector for ESSl)
c and also energy weighted density matrix 
      M15=M13+M
c Least squares
      M17=M15+MM
c
c
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
c      ff(i,1)=0.D0
c      ff(i,2)=0.D0
c      ff(i,3)=0.D0
c 50   continue
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
      alf=a(i,ni)*a(j,nj)/zij
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      ccoef=c(i,ni)*c(j,nj)
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
c
      k=i+((M2-j)*(j-1))/2
c 
c l2: different p in the p shell GRADIENT PART ----------
c
      te=RMM(M15+k-1)*ccoef
      ty=te*2.D0
      t5=ty*a(j,nj)
      t4=ty*a(i,ni)
      con=con+te*ss
c
      do 205 l2=1,3
        t1=Q(l2)-r(Nuc(i),l2)
        pis=t1*ss
        tx=r(Nuc(i),l2)-r(Nuc(j),l2)
        spi=pis+ tx*ss
c
        ff(Nuc(i),l2)=ff(Nuc(i),l2)+t4*pis
        ff(Nuc(j),l2)=ff(Nuc(j),l2)+t5*spi
 205  continue
c
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
c
c
      ccoef=c(i,ni)*c(j,nj)
       t10=ss/z2
c
      do 305 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       ps=ss*t1
c
        ii=i+l1-1
c ii index , taking into account different components of the shell
c
        k=ii+((M2-j)*(j-1))/2
c
        te=RMM(M15+k-1)*ccoef
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
        con=con+te*ps
      do 307 l2=1,3
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(Nuc(j),l2)
       ds=t1*ps
       pp=t2*ps
c
       if (l1.eq.l2) then
        ff(Nuc(i),l2)=ff(Nuc(i),l2)-te*ss
        ds=ds+t10
        pp=pp+t10
       endif
c
c
        ff(Nuc(i),l2)=ff(Nuc(i),l2)+t4*ds
        ff(Nuc(j),l2)=ff(Nuc(j),l2)+t5*pp
c
 307  continue
 305  continue
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
c
c
      t10=ss/z2
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 405 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       t11=pis/z2
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
       t13=spj/z2
c
       pp=t2*pis
c
       if (l1.eq.l2) then
        pp=pp+t10
       endif
c
c
      ii=i+l1-1
      jj=j+l2-1
      k=ii+((M2-jj)*(jj-1))/2
c
        te=RMM(M15+k-1)*ccoef
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
        con=con+te*pp
      do 405 l3=1,3
c
       t1=Q(l3)-r(Nuc(i),l3)
       t2=Q(l3)-r(Nuc(j),l3)
       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
       dp=t1*pp
c
       if (l1.eq.l3) then
        dp=dp+t13
        ff(Nuc(i),l3)=ff(Nuc(i),l3)-te*spj
       endif
c
       if (l2.eq.l3) then
        dp=dp+t11
        ff(Nuc(j),l3)=ff(Nuc(j),l3)-te*pis
       endif
c
       pd=dp+tx*pp
c
        ff(Nuc(i),l3)=ff(Nuc(i),l3)+t4*dp
        ff(Nuc(j),l3)=ff(Nuc(j),l3)+t5*pd
c
 405  continue
c
c
 400  continue
c
c-------------------------------------------------------------------
c (d|s)  gradients
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
c
      t10=ss/z2
      ccoef=c(i,ni)*c(j,nj)
c
      do 505 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       t12=pis/z2
      do 505 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       t11=pjs/z2
c
       dijs=t1*pis
       f1=1.D0
c
       if (l1.eq.l2) then
        f1=sq3
        dijs=dijs+t10
       endif
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       k=ii+((M2-j)*(j-1))/2
       cc=ccoef/f1
c
        te=RMM(M15+k-1)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
        con=con+te*dijs
c gradient part
      do 505 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
       dp=t1*dijs 
c
       if (l1.eq.l3) then
        ff(Nuc(i),l3)=ff(Nuc(i),l3)-te*pjs
        dp=dp+t11
       endif
c
       if (l2.eq.l3) then
        ff(Nuc(i),l3)=ff(Nuc(i),l3)-te*pis
        dp=dp+t12
       endif
c
       fs=dp-tx*dijs
c
        ff(Nuc(i),l3)=ff(Nuc(i),l3)+t4*fs
        ff(Nuc(j),l3)=ff(Nuc(j),l3)+t5*dp
c
 505  continue
c
 500  continue
c-------------------------------------------------------------------
c
c (d|p)  gradients
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
c
      ccoef=c(i,ni)*c(j,nj)
c
      t0=ss/z2
c
      do 605 l1=1,3
c
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       t11=pis/z2
      do 605 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       t12=pjs/z2
c
       f1=1.D0
       dijs=t1*pis
c
       if (l1.eq.l2) then
        f1=sq3
        dijs=dijs+t0
       endif
c
        t13=dijs/z2
       do 605 l3=1,3
c
       t2=Q(l3)-r(Nuc(j),l3)
       pipk=t2*pis
       pjpk=t2*pjs
       dijpk=t2*dijs
c
       if (l1.eq.l3) then
        pipk=pipk+t0
        dijpk=dijpk+t12
       endif
c
       if (l2.eq.l3) then
        pjpk=pjpk+t0
        dijpk=dijpk+t11
       endif
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
       k=ii+((M2-jj)*(jj-1))/2
c
        cc=ccoef/f1
        te=RMM(M15+k-1)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
        con=con+te*dijpk
        t14=pipk/z2
        t15=pjpk/z2
c gradients
       do 605 l4=1,3
c
       t1=Q(l4)-r(Nuc(j),l4)
       tx=r(Nuc(i),l4)-r(Nuc(j),l4)
       dsd=t1*dijpk
c
       if (l1.eq.l4) then
        dsd=dsd+t15
        ff(Nuc(i),l4)=ff(Nuc(i),l4)-te*pjpk
       endif
c
       if (l2.eq.l4) then
        dsd=dsd+t14
       ff(Nuc(i),l4)=ff(Nuc(i),l4)-te*pipk
       endif
c
       if (l3.eq.l4) then
        dsd=dsd+t13
       ff(Nuc(j),l4)=ff(Nuc(j),l4)-te*dijs
       endif
c
       fsp=dsd-tx*dijpk
c
        ff(Nuc(i),l4)=ff(Nuc(i),l4)+t4*fsp
        ff(Nuc(j),l4)=ff(Nuc(j),l4)+t5*dsd
c
 605  continue
c
 600  continue
c-------------------------------------------------------
c (d|d)  gradients
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
c
      ccoef=c(i,ni)*c(j,nj)
c
      t0=ss/z2
c
      do 705 l1=1,3
c
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       t17=pis/z2
      do 705 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       t16=pjs/z2
c
       f1=1.D0
       dijs=t1*pis
c
       if (l1.eq.l2) then
        f1=sq3
        dijs=dijs+t0
       endif
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
       t15=spk/z2
       pipk=t2*pis
       pjpk=t2*pjs
       dijpk=t2*dijs
c
       if (l1.eq.l3) then
        pipk=pipk+t0
        dijpk=dijpk+t16
       endif
c
       if (l2.eq.l3) then
        pjpk=pjpk+t0
        dijpk=dijpk+t17
       endif
c
      t13=dijpk/z2
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
       do 705 l4=1,lk
c
       f2=1.D0
       t1=Q(l4)-r(Nuc(j),l4)
       ovlap=t1*dijpk
       pjdkl=t1*pjpk
       pidkl=t1*pipk
       dijpl=t1*dijs
c
       if (l1.eq.l4) then
        pidkl=pidkl+t15
        dijpl=dijpl+t16
        ovlap=ovlap+pjpk/z2
       endif
c
       if (l2.eq.l4) then
        pjdkl=pjdkl+t15
        dijpl=dijpl+t17
        ovlap=ovlap+pipk/z2
       endif
c
       if (l3.eq.l4) then
        pjdkl=pjdkl+t16
        pidkl=pidkl+t17
        ovlap=ovlap+dijs/z2
        f2=sq3
       endif
c
       t10=pjdkl/z2
       t11=pidkl/z2
       t12=dijpl/z2
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
c
       cc=ccoef/(f1*f2)
        te=RMM(M15+k-1)*cc
        ty=te*2.D0
        t5=ty*a(j,nj)
        t4=ty*a(i,ni)
c
       con=con+te*ovlap
c gradients
        do 705 l5=1,3
        t1=Q(l5)-r(Nuc(i),l5)
        tx=r(Nuc(i),l5)-r(Nuc(j),l5)
c
        fd=t1*ovlap
c
        if (l1.eq.l5) then
        fd=fd+t10
        ff(Nuc(i),l5)=ff(Nuc(i),l5)-te*pjdkl
        endif
c
        if (l2.eq.l5) then
        fd=fd+t11
        ff(Nuc(i),l5)=ff(Nuc(i),l5)-te*pidkl
        endif
c
        if (l3.eq.l5) then
        fd=fd+t12
        ff(Nuc(j),l5)=ff(Nuc(j),l5)-te*dijpl
        endif
c
        if (l4.eq.l5) then
        fd=fd+t13
        ff(Nuc(j),l5)=ff(Nuc(j),l5)-te*dijpk
        endif
c
        df=fd+tx*ovlap
c
        ff(Nuc(i),l5)=ff(Nuc(i),l5)+t4*fd
        ff(Nuc(j),l5)=ff(Nuc(j),l5)+t5*df
c
 705  continue
c
 700  continue
c-------------------------------------------------------------
c
c     write(*,*) 'con =',con
c     do i=1,natom
c      write(*,*) i,ff(i,1),ff(i,2),ff(i,3)
c     enddo
c
c Trick used in order to eliminate motion of center of mass
c due to the errors in exchange-correlation forces
c
c      fx=0.0
c      fy=0.0
c      fz=0.0
c      Pp=0.0
c      do i=1,natom
c       fx=fx+ff(i,1)
c       fy=fy+ff(i,2)
c       fz=fz+ff(i,3)
cc       Pp=Pp+Pm(i)
c      enddo
c
c     write(*,*) 'en intSG'
c     do i=1,natom
c      ff(i,1)=ff(i,1)-Pm(i)*fx/Pp
c      ff(i,2)=ff(i,2)-Pm(i)*fy/Pp
c      ff(i,3)=ff(i,3)-Pm(i)*fz/Pp
c      write(*,*) i,ff(i,1),ff(i,2),ff(i,3)
c     enddo
c     pause
c
      return
      end subroutine
      end module subm_intSG
