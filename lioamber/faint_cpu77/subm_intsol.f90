c-------------------------------------------------------------------
c INTEGRAL CALCULATIONS FOR THE SOLVENT POINT CHARGES ELECTROSTATIC
c INTERACTIONS WITH THE SOLUTE ELECTRONIC DENSITY
c
c Dario Estrin
c Buenos Aires, August. 1994.
c
c 1 e integrals
c using the Obara-Saika recursive method.
c are added into the Fock matrix
c
c It's the same as int.f, but using the solvent atoms partial charges
c
c 
c-------------------------------------------------------------------
      module subm_intsol; contains
      subroutine intsol(E1s,Ens,elec)

      use liotemp   , only: FUNCT
      use garcha_mod, only: RMM, ll, a, c, d, r, pc, nuc, ncont
     >                    , nsol, natom, pi, pi32, rmax, nshell
     >                    , Iz, M, NORM, Md, ntatom
c
      implicit none
      logical :: elec
      real*8  :: Q(3), xi(3)

! Implicits:
      integer :: ns, np, nd
      integer :: i, j, k, ii, jj, ni, nj, j1, j2
      integer :: l, lk, lij, l1, l2, l3, l4, l12, l34
      integer :: MM, MMd, M1, M2, M3, M7, M9, M11

      real*8  :: E1s, Ens, Ese
      real*8  :: sq3, alf, rexp, ccoef, term, temp, tna, tna1, cc
      real*8  :: z2, zij, u, tx, ty, tz
      real*8  :: ss, dd, t1, t2, f1, f2
      real*8  :: dd2, d1s, d2s, d1p, d0s, d0p 
      real*8  :: pj2s, pj1s, pj1p, pj0s, pj0p, pi1p, pi0p
      real*8  :: p3s, p2s, p1s, p0s
c
      real*8, dimension (:), ALLOCATABLE :: s0s,s1s,s2s,s3s,s4s
c distance between pairs of centers
c
c
c      allocate(d(natom,natom))
      allocate(s0s(ntatom),s1s(ntatom),s2s(ntatom),s3s(ntatom)
     > ,s4s(ntatom))
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2
c
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
      M11=M3+MM
c now G
      M7=M11+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c
      E1s=0.0D0
      Ens=0.0D0
      Ese=0.0D0
c      do 50 i=1,natom
c      do 50 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 50   continue

c SOLVENT-SOLUTE ELECTROSTATIC 
c
c        ncerca2=0
c        nlejos=0
c        if(watermod.eq.1) nvecin=2 !tre sites
c        if(watermod.eq.2) nvecin=3 !four sites
c        if(watermod.eq.1.or.watermod.eq.2) then
c      do i=natom+1,natom+nsol
c          ncerca=0
c        do j=natom+1,natom+nsol

c
c          distx=(r(i,1)-r(j,1))**2
c          disty=(r(i,2)-r(j,2))**2
c          distz=(r(i,3)-r(j,3))**2
c          distint=distx+disty+distz
c          if (distint.lt.8.45D0) then
c            ncerca=ncerca+1
c          endif

c         enddo
c          if(ncerca.le.nvecin) then
c            pc(i)=0.D0
c             nlejos=nlejos+1
c           else
c           ncerca2=ncerca2+1
c           endif
c         enddo
c        write(*,*) 'ncerca2=',ncerca2,nlejos,watermod
c         endif


       do 125 j1=1,natom
       do 125 j2=natom+1,nsol+natom
c
c
       tx=r(j1,1)-r(j2,1)
       ty=r(j1,2)-r(j2,2)
       tz=r(j1,3)-r(j2,3)
       dd2=tx**2+ty**2+tz**2
       dd2=sqrt(dd2)
c
        Ens=Ens+Iz(j1)*pc(j2)/dd2

c
 125   continue
       

        if (elec)  then
c
c---------------------------------------------------------------
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      ccoef=c(i,ni)*c(j,nj)
c        write(88,333) i,j,c(i,ni),c(j,nj),Nuc(i),Nuc(j)
       rexp=alf*dd
      if(rexp.lt.rmax) then
     

c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
c
      k=i+((M2-j)*(j-1))/2
c 
       tna=0.D0
      do 202 j1=natom+1,nsol+natom
c
c
c 
       tx=Q(1)-r(j1,1)
       ty=Q(2)-r(j1,2)
       tz=Q(3)-r(j1,3)
c
       u=tx**2 + ty**2 + tz**2
       u=u*zij

       s0s(j1)=pc(j1)*temp*FUNCT(0,u)
       tna=tna-s0s(j1)
 202   continue

c
      term=ccoef*tna
      RMM(M11+k-1)=RMM(M11+k-1)+ term 

      E1s=E1s+RMM(k)*term
       if(E1s.ne.E1s) then
       write(*,*) 'E1s NaN 1'
        stop
c       write(*,*) term, RMM(M11+k-1)
       endif        
c        write(44,*) tna,term,ccoef,RMM(M11+k-1),RMM(k),E1s
        endif
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
       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 302 j1=natom+1,natom+nsol
c
c
c       
c 
       tx=Q(1)-r(j1,1)
       ty=Q(2)-r(j1,2)
       tz=Q(3)-r(j1,3)
c
c
c
       u= tx**2 +ty**2 +tz**2
       u=u*zij

       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
c
 302  continue
c
      ccoef=c(i,ni)*c(j,nj)
c
c l2: different p in the p shell ( x,y,z respectively)
c
      do 305 l2=1,3
        ii=i+l2-1
        t1=Q(l2)-r(Nuc(i),l2)
c ii index , taking into account different components of the shell
        k=ii+((M2-j)*(j-1))/2
c
c loop over nuclei, specific part
       tna=0.D0
      do 303 j1=natom+1,natom+nsol
c
c
       t2=Q(l2)-r(j1,l2)
       term=t1*s0s(j1)-t2*s1s(j1)
       tna=tna-pc(j1)*term
c
 303  continue

        term=ccoef*tna
        RMM(M11+k-1)=RMM(M11+k-1)+term
        E1s=E1s+RMM(k)*term

c       if(E1s.ne.E1s) then
c       write(*,*) 'E1s NaN 2'
c       write(*,*) term, RMM(M11+k-1)
c       endif        
 305    continue
c
      endif
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

       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 402 j1=natom+1,natom+nsol
c
c 

       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij

       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
c
 402  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c loop over partial charges ( specific part)
      do 403 j1=natom+1,natom+nsol
c
c
c
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
c
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)

c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
       do 406 l2=1,lij
       t1=Q(l2)-r(Nuc(j),l2)
       t2=Q(l2)-r(j1,l2)
       tna=t1*p0s-t2*p1s
c
       if (l1.eq.l2) then
        tna=tna+(s0s(j1)-s1s(j1))/z2
       endif
c
        tna1=tna*pc(j1)
c
       ii=i+l1-1
       jj=j+l2-1
       k=ii+((M2-jj)*(jj-1))/2
       term=-tna1*ccoef
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
c       if(E1s.ne.E1s) then
c       write(*,*) 'E1s NaN 3',j1,k
c       write(*,*) term, RMM(M11+k-1)
c       write(*,*) tna1,ccoef
c       write(*,*) pc(j1)
c       stop

c       endif        
 406  continue
c

 403   continue
c ---------------
      endif
 400  continue
c
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
       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over partial charges, part common for all shell
      do 502 j1=natom+1,natom+nsol
c
c
       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
c
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
 502  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 503 j1=natom+1,natom+nsol
c
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
c
      do 506 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(j1,l2)
       tna=t1*p0s-t2*p1s 
c
       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+(s0s(j1)-s1s(j1))/z2
        f1=sq3
       endif
c
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       ii=i+l12-1
c
       k=ii+((M2-j)*(j-1))/2
       cc=ccoef/f1
       term=-cc*tna*pc(j1)
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
c       if(E1s.ne.E1s) then
c       write(*,*) 'E1s NaN 4'
c       write(*,*) term, RMM(M11+k-1)
c       endif        
 506  continue
c
 503  continue
      endif
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
       rexp=alf*dd
       if(rexp.lt.rmax) then
       ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 602 j1=natom+1,natom+nsol
c
c
       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
c
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
       s3s(j1)=temp*FUNCT(3,u)
 602  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 603 j1=natom+1,natom+nsol
c
c
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
c
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(j1,l2)
       pj0s=t1*s0s(j1)-t2*s1s(j1)
       pj1s=t1*s1s(j1)-t2*s2s(j1)
c
       f1=1.D0
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(j1)-s1s(j1))/z2
        d1s=d1s+(s1s(j1)-s2s(j1))/z2
       endif
c
c
      do 606 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(j1,l3)
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
       ii=i+l12-1
       jj=j+l3-1
c
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
       term=-cc*tna*pc(j1)
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
c       if(E1s.ne.E1s) then
c       write(*,*) 'E1s NaN 5'
c       write(*,*) term, RMM(M11+k-1)
c       endif        
 606  continue
c
 603  continue
c
      endif
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
       rexp=alf*dd
      if(rexp.lt.rmax) then
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      do 702 j1=natom+1,natom+nsol
c
       xi(1)=Q(1)-r(j1,1)
       xi(2)=Q(2)-r(j1,2)
       xi(3)=Q(3)-r(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
c
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
       s3s(j1)=temp*FUNCT(3,u)
       s4s(j1)=temp*FUNCT(4,u)
 702  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c Loop over partial charges
      do 703 j1=natom+1,natom+nsol
c
c
c
      do 706 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-r(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
       p3s=t1*s3s(j1)-t2*s4s(j1)
c
      do 706 l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-r(j1,l2)
       pj0s=t1*s0s(j1)-t2*s1s(j1)
       pj1s=t1*s1s(j1)-t2*s2s(j1)
       pj2s=t1*s2s(j1)-t2*s3s(j1)
c
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
c
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(j1)-s1s(j1))/z2
        d1s=d1s+(s1s(j1)-s2s(j1))/z2
        d2s=d2s+(s2s(j1)-s3s(j1))/z2
       endif
c
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
       do 706 l3=1,lij
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-r(j1,l3)
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
        pi0p=pi0p+(s0s(j1)-s1s(j1))/z2
        pi1p=pi1p+(s1s(j1)-s2s(j1))/z2
       endif
c
       if (l2.eq.l3) then
        d0p=d0p+(p0s-p1s)/z2
        d1p=d1p+(p1s-p2s)/z2
        pj0p=pj0p+(s0s(j1)-s1s(j1))/z2
        pj1p=pj1p+(s1s(j1)-s2s(j1))/z2
       endif
c
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      do 706 l4=1,lk
c
       f2=1.D0
       t1=Q(l4)-r(Nuc(j),l4)
       t2=Q(l4)-r(j1,l4)
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
       term=-cc*pc(j1)*tna
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
c       if(E1s.ne.E1s) then
c       write(*,*) 'E1s NaN 6'
c       write(*,*) term, RMM(M11+k-1)
c       endif        
c
 706  continue
 703  continue
      endif
 700  continue


       endif
      deallocate(s0s,s2s,s3s,s4s)

c
 333  format(2(I4,2x),2(F10.4,2x),2(I4,2x))
c       write(*,*) 'E1s=',E1s,Ens
c       if(E1s.ne.E1s) then
c       write(33,*) RMM(1:MM)
c       stop
c       endif       
      return
      end subroutine
      end module subm_intsol
c-------------------------------------------------------------------
