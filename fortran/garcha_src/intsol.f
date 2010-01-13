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
      subroutine intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,E1s)
c
      implicit real*8 (a-h,o-z)
      logical NORM
      INCLUDE 'param'
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(natom)
      dimension r(nt,3),nshell(0:3),pc(nss)
      dimension xi(3)
      dimension RMM(*)
c
      COMMON /TABLE/ STR(880,0:21)
c auxiliar quantities
c
      dimension Q(3),d(ntq,ntq),s0s(nt),s1s(nt),s2s(nt),s3s(nt),
     >   s4s(nt),Ll(3),rr(3)
c distance between pairs of centers
c
c datos TIP4P ----------------
c corresponde a desplazar carga negativa 0.15A en direccion de los H
c      alpha=0.7439762D0
c      alpha2=0.1280119D0
c caso SPC
       alpha=1.00D0
       alpha2=0.00D0
c -----------------------------
c
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
      do 50 i=1,natom
      do 50 j=1,natom
       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
 50   continue
c
c---------------------------------------------------------------
c
      E1s=0.0D0
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
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
c
      k=i+((M2-j)*(j-1))/2
c 
       tna=0.D0
      j1=natom
      do 202 n1=1,Nsol
      do 202 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
c
c test
       dp=(rr(1)-r(j1,1))**2+(rr(2)-r(j1,2))**2+(rr(3)-r(j1,3))**2
       dp=sqrt(dp)
c      write(*,*) 'dd',dp*0.529177
c
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c       
c 
       tx=Q(1)-rr(1)
       ty=Q(2)-rr(2)
       tz=Q(3)-rr(3)
c
       u=tx**2 + ty**2 + tz**2
       u=u*zij

       s0s(j1)=pc(k1)*temp*FUNCT(0,u)
       tna=tna-s0s(j1)
c
 202   continue
c
      term=ccoef*tna
      RMM(M11+k-1)=RMM(M11+k-1)+ term 
      E1s=E1s+RMM(k)*term
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      j1=natom
      do 302 n1=1,Nsol
      do 302 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c       
c 
       tx=Q(1)-rr(1)
       ty=Q(2)-rr(2)
       tz=Q(3)-rr(3)
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
      j1=natom
      do 303 n1=1,Nsol
      do 303 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c       
c 
c
       t2=Q(l2)-rr(l2)
       term=t1*s0s(j1)-t2*s1s(j1)
       tna=tna-pc(k1)*term
c
 303  continue

        term=ccoef*tna
        RMM(M11+k-1)=RMM(M11+k-1)+term
        E1s=E1s+RMM(k)*term
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      j1=natom
      do 402 n1=1,Nsol
      do 402 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c       
c 

       xi(1)=Q(1)-rr(1)
       xi(2)=Q(2)-rr(2)
       xi(3)=Q(3)-rr(3)
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
      j1=natom
      do 403 n1=1,Nsol
      do 403 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c
c
c
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
c
       t2=Q(l1)-rr(l1)
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
       t2=Q(l2)-rr(l2)
       tna=t1*p0s-t2*p1s
c
       if (l1.eq.l2) then
        tna=tna+(s0s(j1)-s1s(j1))/z2
       endif
c
        tna1=tna*pc(k1)
c
       ii=i+l1-1
       jj=j+l2-1
       k=ii+((M2-jj)*(jj-1))/2
       term=-tna1*ccoef
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
 406  continue
c

 403   continue
c ---------------
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over partial charges, part common for all shell
      j1=natom
      do 502 n1=1,Nsol
      do 502 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c
       xi(1)=Q(1)-rr(1)
       xi(2)=Q(2)-rr(2)
       xi(3)=Q(3)-rr(3)
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
      j1=natom
      do 503 n1=1,Nsol
      do 503 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-rr(l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
c
      do 506 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-rr(l2)
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
       term=-cc*tna*pc(k1)
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
 506  continue
c
 503  continue
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      j1=natom
      do 602 n1=1,Nsol
      do 602 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c
       xi(1)=Q(1)-rr(1)
       xi(2)=Q(2)-rr(2)
       xi(3)=Q(3)-rr(3)
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
      j1=natom
      do 603 n1=1,Nsol
      do 603 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c
c
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-rr(l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
c
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-rr(l2)
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
       t2=Q(l3)-rr(l3)
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
       term=-cc*tna*pc(k1)
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
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
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
c
c loop over nuclei, part common for all shell
      j1=natom
      do 702 n1=1,Nsol
      do 702 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c
       xi(1)=Q(1)-rr(1)
       xi(2)=Q(2)-rr(2)
       xi(3)=Q(3)-rr(3)
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
      j1=natom
      do 703 n1=1,Nsol
      do 703 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c
c
      do 706 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-rr(l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
       p3s=t1*s3s(j1)-t2*s4s(j1)
c
      do 706 l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-rr(l2)
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
       t2=Q(l3)-rr(l3)
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
       t2=Q(l4)-rr(l4)
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
       term=-cc*pc(k1)*tna
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
c
 706  continue
 703  continue
 700  continue
c
c     write(*,*) 'E1s=',E1s
      return
      end
c-------------------------------------------------------------------
