c-------------------------------------------------------------------
c Calculation of electrical potential for an arbitrary point of
c the space
c 1 e integrals
c using the Obara-Saika recursive method.
c Y calcula la densidad en los mismos puntos.
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
c Output: F matrix, and S matrix
c-------------------------------------------------------------------
      subroutine elec(NX,NY,NZ,deltax,xmin,ymin,zmin)
      use garcha_mod
c
      implicit real*8 (a-h,o-z)
c auxiliar quantities
      dimension Q(3),xi(3)
      real*8, dimension (:), ALLOCATABLE :: pote
c      parameter (deltax=.1D0)      
c      parameter (buffer=5.D0)
      parameter (Ar0=0.529177D0)
      real*8 VVV,deltax
      integer n_graf
c
      rmax=3.
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
      M15=M11+3*Npoint+1
c
      do 50 i=1,natom
      do 50 j=1,natom
       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
 50   continue
c
       Ntot=NX*NY*NZ
       allocate(pote(Ntot))

       pote=0

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
       ntotal=0
       do 202 nnx=1,NX
       do 202 nny=1,NY
       do 202 nnz=1,NZ

        ntotal=ntotal+1

        xi(1)=Xmin+(nnx-1)*deltax
        xi(2)=Ymin+(nny-1)*deltax
        xi(3)=Zmin+(nnz-1)*deltax
        u=(Q(1)-xi(1))**2+(Q(2)-xi(2))**2+(Q(3)-xi(3))**2
        u=u*zij
        s0s=temp*FUNCT(0,u)
c
      term=ccoef*s0s*RMM(kk)
      pote(ntotal)=pote(ntotal)+ term 

 202  continue
 201  continue
       write(*,*) 'un ss hecho'
 200  continue
       write(*,*) 'elec ss hecho'
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
       ntotal=0
       do 305 nnx=1,NX
       do 305 nny=1,NY
       do 305 nnz=1,NZ

        ntotal=ntotal+1

        xi(1)=Xmin+(nnx-1)*deltax
        xi(2)=Ymin+(nny-1)*deltax
        xi(3)=Zmin+(nnz-1)*deltax
c
       u=(Q(1)-xi(1))**2+(Q(2)-xi(2))**2+(Q(3)-xi(3))**2
       u=u*zij
       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
c
c
c l2: different p in the p shell ( x,y,z respectively)
c
      do 305 l2=1,3
        t1=Q(l2)-r(Nuc(i),l2)
        ii=i+l2-1
c ii index , taking into account different components of the shell
c
        k=ii+((M2-j)*(j-1))/2
c
        tna=t1*s0s-(Q(l2)-xi(l2))*s1s

        kk=ii+((M2-j)*(j-1))/2
        term=ccoef*tna*RMM(kk)

        pote(ntotal)=pote(ntotal)+term
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
       ntotal=0
       do 406 nnx=1,NX
       do 406 nny=1,NY
       do 406 nnz=1,NZ
        
        ntotal=ntotal+1
        
        xi(1)=Xmin+(nnx-1)*deltax
        xi(2)=Ymin+(nny-1)*deltax
        xi(3)=Zmin+(nnz-1)*deltax
c
       u=(Q(1)-xi(1))**2+(Q(2)-xi(2))**2+(Q(3)-xi(3))**2
       u=u*zij

       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
c
c
c
c
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-xi(l1)
       p0s=t1*s0s-t2*s1s
       p1s=t1*s1s-t2*s2s
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
       if (l1.eq.l2) then
        tna=tna+(s0s-s1s)/z2
       endif
c
c
       ii=i+l1-1
       jj=j+l2-1
       k=ii+((M2-jj)*(jj-1))/2
       term=tna*ccoef*RMM(k)
      pote(ntotal)=pote(ntotal)+term
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

       ntotal=0
       do 506 nnx=1,NX
       do 506 nny=1,NY
       do 506 nnz=1,NZ
        
        ntotal=ntotal+1
        
        xi(1)=Xmin+(nnx-1)*deltax
        xi(2)=Ymin+(nny-1)*deltax
        xi(3)=Zmin+(nnz-1)*deltax
c
       u=(Q(1)-xi(1))**2+(Q(2)-xi(2))**2+(Q(3)-xi(3))**2
       u=u*zij
       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
c
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-xi(l1)
       p0s=t1*s0s-t2*s1s
       p1s=t1*s1s-t2*s2s
c
      do 506 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-xi(l2)
       tna=t1*p0s-t2*p1s
c
       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+(s0s-s1s)/z2
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
       term=cc*tna*RMM(k)
      pote(ntotal)=pote(ntotal)+term
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
       ntotal=0
       do 606 nnx=1,NX
       do 606 nny=1,NY
       do 606 nnz=1,NZ
        
        ntotal=ntotal+1
        
        xi(1)=Xmin+(nnx-1)*deltax
        xi(2)=Ymin+(nny-1)*deltax
        xi(3)=Zmin+(nnz-1)*deltax
c
       u=(Q(1)-xi(1))**2+(Q(2)-xi(2))**2+(Q(3)-xi(3))**2
       u=u*zij
c
       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
       s3s=temp*FUNCT(3,u)
c
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=Q(l1)-xi(l1)
       p0s=t1*s0s-t2*s1s
       p1s=t1*s1s-t2*s2s
       p2s=t1*s2s-t2*s3s
c
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=Q(l2)-xi(l2)
       pj0s=t1*s0s-t2*s1s
       pj1s=t1*s1s-t2*s2s
c
       f1=1.D0
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s-s1s)/z2
        d1s=d1s+(s1s-s2s)/z2
       endif
c
      do 606 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=Q(l3)-xi(l3)
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
       term=cc*tna*RMM(k)
      pote(ntotal)=pote(ntotal)+term
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
       ntotal=0
       do 706 nnx=1,NX
       do 706 nny=1,NY
       do 706 nnz=1,NZ

        ntotal=ntotal+1

        xi(1)=Xmin+(nnx-1)*deltax
        xi(2)=Ymin+(nny-1)*deltax
        xi(3)=Zmin+(nnz-1)*deltax
c
c
       u=(Q(1)-xi(1))**2+(Q(2)-xi(2))**2+(Q(3)-xi(3))**2
       u=u*zij

       s0s=temp*FUNCT(0,u)
       s1s=temp*FUNCT(1,u)
       s2s=temp*FUNCT(2,u)
       s3s=temp*FUNCT(3,u)
       s4s=temp*FUNCT(4,u)
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
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s-s1s)/z2
        d1s=d1s+(s1s-s2s)/z2
        d2s=d2s+(s2s-s3s)/z2
       endif
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
       if (l1.eq.l3) then
        d0p=d0p+(pj0s-pj1s)/z2
        d1p=d1p+(pj1s-pj2s)/z2
        pi0p=pi0p+(s0s-s1s)/z2
        pi1p=pi1p+(s1s-s2s)/z2
       endif
c
       if (l2.eq.l3) then
        d0p=d0p+(p0s-p1s)/z2
        d1p=d1p+(p1s-p2s)/z2
        pj0p=pj0p+(s0s-s1s)/z2
        pj1p=pj1p+(s1s-s2s)/z2
       endif
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
       term=cc*tna
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
      pote(ntotal)=pote(ntotal) + term*RMM(k)
c
c
 706  continue
 701  continue
 700  continue
cccccccccccccc Ahora calculamos la densidad electrónica en la grilla
       ntotal=0
       do 806 nnx=1,NX
       do 806 nny=1,NY
       do 806 nnz=1,NZ
        ntotal=ntotal+1

        xi(1)=Xmin+(nnx-1)*deltax
        xi(2)=Ymin+(nny-1)*deltax
        xi(3)=Zmin+(nnz-1)*deltax

       Vvv=0.0D0
       do 11 n=1,natom
       dd=(xi(1)-r(n,1))**2+(xi(2)-r(n,2))**2+(xi(3)-r(n,3))**2
       dd=sqrt(dd)
       Vvv=Vvv-Iz(n)/dd
  11   continue
      pote(ntotal)=pote(ntotal) + Vvv 
 806   continue
c       else

c
c
C     WRITE A FORMATTED CUBEFILE VERY SIMILAR TO THOSE CREATED BY 
C     THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
C     THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
C
C     LINE   FORMAT      CONTENTS
C     ===============================================================
C      1     A           TITLE
C      2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
C      3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
C      4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
C      #ATOMS LINES OF ATOM COORDINATES:
C      ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
C      REST: 6E13.5      CUBE DATA (WITH Z INCREMENT MOVING FASTEST, THEN
C                        Y AND THEN X)
C
C     FOR ORBITAL CUBE FILES, #ATOMS WILL BE < 0 AND THERE WILL BE ONE
C     ADDITIONAL LINE AFTER THE FINAL ATOM GIVING THE NUMBER OF ORBITALS
C     AND THEIR RESPECTIVE NUMBERS. ALSO THE ORBITAL NUMBER WILL BE
C     THE FASTEST MOVING INCREMENT.
C
C     ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
       open(unit=188,file=cube_elec_file)

       write(188,*) 'Elfield'
       write(188,*)
       write(188,678) natom,Xmin,Ymin,Zmin
       write(188,678) NX,deltax,0.0D0,0.0D0 
       write(188,678) NY,0.0D0, deltax,0.0D0 
       write(188,678) NZ,0.0D0,0.0D0,deltax
       do i=1,natom
       write(188,677) Iz(i),0.0D0,r(i,1),r(i,2),r(i,3)
       enddo 
       write(188,679) pote
       close(188)
677    format(I5,4(F12.6))
678    format(I5,3(F12.6))
679    format(6(E13.5))

         deallocate(pote)
      return
      end
c-------------------------------------------------------------------


 
