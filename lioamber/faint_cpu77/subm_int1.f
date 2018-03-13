!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      module subm_int1; contains
       subroutine int1(En)
!------------------------------------------------------------------------------!
c
c      Integrals subroutine 
c      1 e integrals
c      using the Obara-Saika recursive method.
c
c 
c      loop over all basis functions
c      now the basis is supposed to be ordered according to the type,
c      all s, then all p, then all d, .....
c      inside each type, are ordered in shells
c      px,py,pz , dx2,dxy,dyy,dzx,dzy,dzz, .....
c
c      ns ... marker for end of s
c      np ... marker for end of p
c      nd ... marker for end of d
c
c      r(Nuc(i),j) j component of position of nucleus i , j=1,3
c      Input : basis function information
c      Output: F matrix, and S matrix
c
!------------------------------------------------------------------------------!
       use liotemp   , only: FUNCT
       use garcha_mod, only: RMM, Smat, Nuc, a, c, d, r, Iz, ncont
     >                     , nshell, pi, pi32, NORM, natom, M, Md
       implicit none

c-----auxiliar quantities
       real*8, intent(inout) :: En

       integer :: natomold, igpu
       integer :: n, i, j, k, ii, jj, ni, nj, iin
       integer :: l1, l2, l3, l4, l12, l34
       integer :: MM, MMd, ns, np, nd
       integer :: M1, M2, M3, M5, M7, M9, M11

       real*8  :: E1, ovlap
       real*8  :: Q(3), term, temp, sq3, alf, alf2, cc, ccoef
       real*8  :: f1, f2, tn, tna, u, z2, zij
       real*8  :: ss, ps, dd, p0s, p1s, p2s, p3s
       real*8  :: pi0p, pi1p, piks, pikpk, pipk, pis
       real*8  :: pj0s, pj1s, pj2s, pj0p, pj1p, pjkpk
       real*8  :: pjks, pjpk, pjs, pks, sks
       real*8  :: dijs, dijpk, dijks, dijkpk
       real*8  :: d0s, d0p, d1p, d1s, d2s
       real*8  :: t0, t1, t2

       integer, allocatable, dimension(:) :: Iaux
       real*8 , allocatable, dimension(:) :: s0s, s1s, s2s, s3s, s4s
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       allocate(s0s(natom),s1s(natom),s2s(natom))
!       allocate(s3s(natom),s4s(natom),Iaux(natom))
       allocate(s3s(natom))
       allocate(s4s(natom))
!       allocate(Iaux(natom))
       if (.not.allocated(Smat)) allocate(Smat(M,M))

c-----distance between pairs of centers
c
c --- BSSE------------------------------------
c      sets to 0 atomic charges, but since they are used to
c      construct the grids, they are stored in an auxiliar array
c      if (BSSE) then
c      do i=1,natom
c       Iaux(i)=Iz(i)
c       Iz(i)=Iz(i)*ighost(i)
c      enddo
c      endif
c----------------------------------------------

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
c
c Pointers
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, F also uses the same position after S was used
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

       Smat=0.0d0
       do i=1,MM
         RMM(M5+i-1)=0.D0
         RMM(M11+i-1)=0.D0
       enddo
c
c      do 50 i=1,natom
c      do 50 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 50   continue
c
c Nuclear Repulsion part ------------------------------------------
      En=0.D0

       do i=1,natom
       do j=1,i-1
         En=En+Iz(i)*Iz(j)/sqrt(d(i,j))
       enddo
       enddo

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
      alf=a(i,ni)*a(j,nj)/zij
      alf2=alf*2.D0
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      ccoef=c(i,ni)*c(j,nj)
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ovlap=ss
      tn=alf*(3.D0-alf2*dd)*ovlap
c
      k=i+((M2-j)*(j-1))/2
      RMM(M5+k-1)=RMM(M5+k-1)+ccoef*ovlap
      Smat(i,j)=Smat(i,j)+ccoef*ovlap
      Smat(j,i)=Smat(j,i)+ccoef*ovlap
c 
c loop over nuclei, nuclear attraction matrix elements
c tna: accumulates nuc. attraction over all nuclei
c
       tna=0.D0

      do n=1,natom
       u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
       u=u*zij

       s0s(n)=Iz(n)*2.*sqrt(zij/pi)*ss*FUNCT(0,u)
       tna=tna-s0s(n)
      enddo
c
      term=ccoef*(tn+tna)
      RMM(M11+k-1)=RMM(M11+k-1)+ term
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
c
c loop over nuclei, part common for all shell
      do 302 n=1,natom
       u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
       u=u*zij
       temp=2.D0*sqrt(zij/pi)*ss

       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
 302  continue
c
      ccoef=c(i,ni)*c(j,nj)
c
c l2: different p in the p shell ( x,y,z respectively)
c
      do 305 l2=1,3
        t1=Q(l2)-r(Nuc(i),l2)
        ovlap=t1*ss
        tn=t1*sks+alf2*ovlap
        iin=i+l2-1
c ii index , taking into account different components of the shell
c
        k=iin+((M2-j)*(j-1))/2
        RMM(M5+k-1)=RMM(M5+k-1)+ovlap*ccoef
        Smat(iin,j)=Smat(iin,j)+ovlap*ccoef
        Smat(j,iin)=Smat(j,iin)+ovlap*ccoef
c
c loop over nuclei, specific part
       tna=0.D0
      do 303 n=1,natom
       term=t1*s0s(n)-(Q(l2)-r(n,l2))*s1s(n)
       tna=tna-Iz(n)*term
 303  continue

        term=ccoef*(tn+tna)
        RMM(M11+k-1)=RMM(M11+k-1)+term
 305    continue
c

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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
c
c loop over nuclei, part common for all shell
      do 402 n=1,natom
       u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
       u=u*zij
       temp=2.D0*sqrt(zij/pi)*ss

       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
 402  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 405 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       ps=ss*t1
       pks=sks*t1+alf2*ps
c
       do 405 l2=1,3
c
       t1=Q(l2)-r(Nuc(j),l2)
       ovlap=t1*ps
       tn=t1*pks+alf2*ovlap
c
       if (l1.eq.l2) then
        ovlap=ovlap+ss/z2
        tn=tn+(sks+alf2*ss)/z2
       endif
c
       iin=i+l1-1
       jj=j+l2-1
       term=tn*ccoef
c
       if (iin.ge.jj) then
         k=iin+((M2-jj)*(jj-1))/2
         RMM(M5+k-1)=RMM(M5+k-1)+ovlap*ccoef
         Smat(iin,jj)=Smat(iin,jj)+ovlap*ccoef
         Smat(jj,iin)=Smat(jj,iin)+ovlap*ccoef
         RMM(M11+k-1)=RMM(M11+k-1)+term
       endif
 405  continue
c
c loop over nuclei ( specific part)
      do 403 n=1,natom
c
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
c
       do 406 l2=1,3
c
       t1=Q(l2)-r(Nuc(j),l2)
       t2=Q(l2)-r(n,l2)
       tna=t1*p0s-t2*p1s
c
       if (l1.eq.l2) then
        tna=tna+(s0s(n)-s1s(n))/z2
       endif
c
c
       iin=i+l1-1
       jj=j+l2-1
       if (iin.ge.jj) then
         k=iin+((M2-jj)*(jj-1))/2
         term=-tna*ccoef*Iz(n)
         RMM(M11+k-1)=RMM(M11+k-1)+term
       endif
 406  continue
 403   continue
c ---------------
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
c
c loop over nuclei, part common for all shell
      do 502 n=1,natom
       u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
       u=u*zij
       temp=2.D0*sqrt(zij/pi)*ss
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
 502  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 505 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       ps=ss*t1
       pks=sks*t1+alf2*ps
c
      do 505 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       ovlap=t1*ps
       tn=t1*pks
c
       f1=1.D0
       if (l1.eq.l2) then
        ovlap=ovlap+ss/z2
        tn=tn+sks/z2-alf2*ss/(2.*a(i,ni))
        f1=sq3
       endif
c
       tn=tn+alf2*ovlap
c
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       iin=i+l12-1
c
       cc=ccoef/f1
       term=cc*tn
       k=iin+((M2-j)*(j-1))/2
       RMM(M5+k-1)=RMM(M5+k-1)+ovlap*cc
       Smat(iin,j)=Smat(iin,j)+ovlap*cc
       Smat(j,iin)=Smat(j,iin)+ovlap*cc
       RMM(M11+k-1)=RMM(M11+k-1)+term
 505  continue
cc nuclear attraction part 
c
      do 503 n=1,natom
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-R(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
c
      do 506 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       tna=t1*p0s-t2*p1s 
c
       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+(s0s(n)-s1s(n))/z2
        f1=sq3
       endif
c
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       iin=i+l12-1
c
       k=iin+((M2-j)*(j-1))/2
       cc=ccoef/f1
       term=-cc*tna*Iz(n)
       RMM(M11+k-1)=RMM(M11+k-1)+term
 506  continue
c
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
c
c loop over nuclei, part common for all shell
      do 602 n=1,natom
       u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
       u=u*zij
       temp=2.D0*sqrt(zij/pi)*ss
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
       s3s(n)=temp*FUNCT(3,u)
 602  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 605 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       piks=sks*t1+alf2*pis
      do 605 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       pjks=sks*t1+alf2*pjs
c
       dijs=t1*pis
       dijks=t1*piks
       f1=1.D0
c
       if (l1.eq.l2) then
        f1=sq3
        dijs=dijs+ss/z2
        dijks=dijks+sks/z2-alf2*ss/(2.D0*a(i,ni))
       endif
c
       dijks=dijks+alf2*dijs
c
      do 605 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       ovlap=t1*dijs 
       tn=t1*dijks
c
       if (l1.eq.l3) then
        ovlap=ovlap+pjs/z2
        tn=tn+pjks/z2
       endif
c
       if (l2.eq.l3) then
        ovlap=ovlap+pis/z2
        tn=tn+piks/z2
       endif
c
       tn=tn+alf2*ovlap
       
c
       l12=l1*(l1-1)/2+l2
       iin=i+l12-1
       jj=j+l3-1
c
       cc=ccoef/f1
       term=cc*tn
       k=iin+((M2-jj)*(jj-1))/2
       RMM(M5+k-1)=RMM(M5+k-1)+cc*ovlap
       Smat(iin,jj)=Smat(iin,jj)+cc*ovlap
       Smat(jj,iin)=Smat(jj,iin)+cc*ovlap
       RMM(M11+k-1)=RMM(M11+k-1)+term
 605  continue
c
c Nuclear attraction part ----------
      do 603 n=1,natom
c
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       p2s=t1*s2s(n)-t2*s3s(n)
c
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
c
       f1=1.D0
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(n)-s1s(n))/z2
        d1s=d1s+(s1s(n)-s2s(n))/z2
       endif
c
c
      do 606 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(n,l3)
       tna=t1*d0s-t2*d1s
c
       if (l1.eq.l3) then
        tna=tna+(pj0s-pj1s)/z2
       endif
c
       if (l2.eq.l3) then
        tna=tna+(p0s-p1s)/z2
       endif
c
c
c
       l12=l1*(l1-1)/2+l2
       iin=i+l12-1
       jj=j+l3-1
c
      k=iin+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
       term=-cc*tna*Iz(n)
       RMM(M11+k-1)=RMM(M11+k-1)+term
 606  continue
c
 603  continue
c
 600  continue
c
c-------------------------------------------------------------------
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sks=alf*(3.D0-alf2*dd)*ss
c
c loop over nuclei, part common for all shell
      do 702 n=1,natom
       u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
       u=u*zij
       temp=2.D0*sqrt(zij/pi)*ss
       s0s(n)=temp*FUNCT(0,u)
       s1s(n)=temp*FUNCT(1,u)
       s2s(n)=temp*FUNCT(2,u)
       s3s(n)=temp*FUNCT(3,u)
       s4s(n)=temp*FUNCT(4,u)
 702  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      t0=ss/z2
c
      do 705 l1=1,3
c
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       piks=sks*t1+alf2*pis
      do 705 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       pjks=sks*t1+alf2*pjs
c
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       dijs=t1*pis
       dijks=t1*piks
c
       if (l1.eq.l2) then
        f1=sq3
        dijs=dijs+t0
        dijks=dijks+sks/z2-alf2*ss/(2.D0*a(i,ni))
       endif
c
       dijks=dijks+alf2*dijs
c
       do 705 l3=1,3
c
       t2=Q(l3)-r(Nuc(j),l3)
       pipk=t2*pis
       pjpk=t2*pjs
       pikpk=t2*piks
       pjkpk=t2*pjks
       dijpk=t2*dijs
       dijkpk=t2*dijks
c
       if (l1.eq.l3) then
        pipk=pipk+t0
        dijpk=dijpk+pjs/z2
        pikpk=pikpk+sks/z2
        dijkpk=dijkpk+pjks/z2
       endif
c
       if (l2.eq.l3) then
        pjpk=pjpk+t0
        dijpk=dijpk+pis/z2
        pjkpk=pjkpk+sks/z2
        dijkpk=dijkpk+piks/z2
       endif
c
       pikpk=pikpk+alf2*pipk
       pjkpk=pjkpk+alf2*pjpk
       dijkpk=dijkpk+alf2*dijpk
c 
       do 705 l4=1,l3
c
       f2=1.D0
       t1=Q(l4)-r(Nuc(j),l4)
       ovlap=t1*dijpk
       tn=t1*dijkpk
c
       if (l1.eq.l4) then
        ovlap=ovlap+pjpk/z2
        tn=tn+pjkpk/z2
       endif
c
       if (l2.eq.l4) then
        ovlap=ovlap+pipk/z2
        tn=tn+pikpk/z2
       endif
c
       if (l3.eq.l4) then
        ovlap=ovlap+dijs/z2
        tn=tn+dijks/z2-alf2*dijs/(2.D0*a(j,nj))
        f2=sq3
       endif
c
c l12 and l34 go from 1 to 6, spanning the d shell in
c the order xx, xy, yy, zx, zy, zz. The same order should be used
c in ordering the basis set, before calculating the matrix elements
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       iin=i+l12-1
       jj=j+l34-1
c
       if (iin.ge.jj) then
       k=iin+((M2-jj)*(jj-1))/2
       tn= tn+alf2*ovlap
c
       cc=ccoef/(f1*f2)
       term=cc*tn
       RMM(M5+k-1)=RMM(M5+k-1)+ovlap*cc
       Smat(iin,jj)=Smat(iin,jj)+ovlap*cc
       Smat(jj,iin)=Smat(jj,iin)+ovlap*cc
       RMM(M11+k-1)=RMM(M11+k-1)+term
       endif
c
 705  continue
c
c Loop over nuclei - Nuclear attraction part ---
      do 703 n=1,natom
c
      do 706 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(n,l1)
       p0s=t1*s0s(n)-t2*s1s(n)
       p1s=t1*s1s(n)-t2*s2s(n)
       p2s=t1*s2s(n)-t2*s3s(n)
       p3s=t1*s3s(n)-t2*s4s(n)
c
      do 706 l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(n,l2)
       pj0s=t1*s0s(n)-t2*s1s(n)
       pj1s=t1*s1s(n)-t2*s2s(n)
       pj2s=t1*s2s(n)-t2*s3s(n)
c
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
c
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(n)-s1s(n))/z2
        d1s=d1s+(s1s(n)-s2s(n))/z2
        d2s=d2s+(s2s(n)-s3s(n))/z2
       endif
c
c
      do 706 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(n,l3)
c
       d0p=t1*d0s-t2*d1s
       d1p=t1*d1s-t2*d2s
c
       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s
c
       if (l1.eq.l3) then
        d0p=d0p+(pj0s-pj1s)/z2
        d1p=d1p+(pj1s-pj2s)/z2
        pi0p=pi0p+(s0s(n)-s1s(n))/z2
        pi1p=pi1p+(s1s(n)-s2s(n))/z2
       endif
c
       if (l2.eq.l3) then
        d0p=d0p+(p0s-p1s)/z2
        d1p=d1p+(p1s-p2s)/z2
        pj0p=pj0p+(s0s(n)-s1s(n))/z2
        pj1p=pj1p+(s1s(n)-s2s(n))/z2
       endif
c
c
       do 706 l4=1,l3
c
       f2=1.D0
       t1=Q(l4)-R(Nuc(j),l4)
       t2=Q(l4)-R(n,l4)
       tna=t1*d0p-t2*d1p
c
       if (l4.eq.l1) then
        tna=tna+(pj0p-pj1p)/z2
       endif
c
       if (l4.eq.l2) then
        tna=tna+(pi0p-pi1p)/z2
       endif
c
       if (l4.eq.l3) then
        f2=sq3
        tna=tna+(d0s-d1s)/z2
       endif
c
       cc=ccoef/(f1*f2)
       term=-cc*Iz(n)*tna
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       iin=i+l12-1
       jj=j+l34-1
c
       if (iin.ge.jj) then
       k=iin+((M2-jj)*(jj-1))/2
       RMM(M11+k-1)=RMM(M11+k-1)+term
       endif
c
 706  continue
c
 703  continue
c end nuclear attraction part --------
c
 700  continue
c
c
c--- Debugging and tests -----------------------------------------------
c
c     do 1 k=1,M*(M+1)/2
c      write(*,*) k,RMM(M5+k-1),RMM(M11+k-1)
c gradient debuggings
c
         E1=0.D0

        do 802 k=1,MM
 802     E1=E1+RMM(k)*RMM(M11+k-1)
c
c      write(*,*) 'E1+En =',E1+En
c
c BSSE ------------
c      if (BSSE) then
c      do i=1,natom
c       Iz(i)=Iaux(i)
c      enddo
c      endif
c
c-- prueba ----------------
c      En=En+0.0D0*(d(1,2)-2.89D0)**2
c--------------------------
c
c     write(*,*) 'matriz overlap'
c     do i=1,MM
c      write(*,*) i,RMM(M5+i-1)
c     enddo
c     do i=1,natom
c      write(*,*) i,r(i,1),r(i,2),r(i,3)
c     enddo
c     pause
      do i=1,M
        Smat(i,i)=Smat(i,i)/2
      enddo

      deallocate(s0s,s1s,s2s,s3s,s4s)
c      deallocate(Iaux) 


      if (igpu.gt.3) natom = natomold

      return;end subroutine
      end module subm_int1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
