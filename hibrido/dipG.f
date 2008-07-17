c-------------------------------------------------------------------
c
c Subroutine for evaluating the contributions to the
c gradients of reaction field terms
c See J. Compu. Chem. 13 675 (1992) and
c JACS 113 4776 (1991)
c
c integrals evaluated using Obara Saika method
c Buenos Aires, 11 Enero 1994
c Debugged hasta funciones p
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
c Input :  density basis 
c Output:  dipole moment components 
c-----------------------------------------------------------------
      subroutine dipg(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >               Nel,g,ux,uy,uz,f)
c
      implicit real*8 (a-h,o-z)
      logical NORM
      INCLUDE 'param'
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0,pi5=34.9868366552497108D0)
      dimension r(nt,3),nshell(0:3),Iz(natom)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M)
c
c aux . things
      dimension aux(3),aux1(3),aux2(3),aux3(3),aux4(3),aux5(3)
      dimension aux6(3),aux7(3),aux8(3),aux9(3),aux10(3)
      dimension auxd(3,3)
      dimension Q(3),d(ntq,ntq)
      dimension RMM(*)
c
      dimension f(nt,3)
c
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
c  
      do 50 i=1,natom
      do 50 j=1,natom
       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
 50   continue
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
       Qc=0.0D0
      do 99 i=1,natom
       Qc=Qc+Iz(i)
 99   continue
c
      Qc=Qc-Nel
c fac : for charged species dipole moment depends on the origin
c definition - Using the factor, it is defined with respect to the
c center of charge (important in Reaction Field calculations)
c For neutral systems it doesn't matter
c
      fac=(Qc+Nel)/Nel
c
      g=g/2.54D0
      g1=g*fac
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
      z2=2.0D0*zij
      alf=a(i,ni)*a(j,nj)/zij
      ccoef=c(i,ni)*c(j,nj)
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ssz2=ss/z2
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
      k=i+((M2-j)*(j-1))/2
      cc=ccoef*RMM(k)
c
c derivatives calculation
       ti=2.0D0*a(i,ni)
       tj=2.0D0*a(j,nj)
       tx=Q(1)-r(Nuc(i),1)
       ty=Q(2)-r(Nuc(i),2)
       tz=Q(3)-r(Nuc(i),3)
c
       txj=Q(1)-r(Nuc(j),1)
       tyj=Q(2)-r(Nuc(j),2)
       tzj=Q(3)-r(Nuc(j),3)
c
       cci=cc*ti
       ccx=cci*tx
       ccy=cci*ty
       ccz=cci*tz
cc
c 1-3 derivatives of ux, 4-6 derivatives of uy, 7-9 of uz
       auxd(1,1)=ccx*sxs
       auxd(1,2)=ccy*sxs
       auxd(1,3)=ccz*sxs
c
       auxd(2,1)=ccx*sys
       auxd(2,2)=ccy*sys
       auxd(2,3)=ccz*sys
c
       auxd(3,1)=ccx*szs
       auxd(3,2)=ccy*szs
       auxd(3,3)=ccz*szs
c
       t1=cci*ssz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       f(Nuc(i),1)=f(Nuc(i),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     > uz*auxd(3,1))
       f(Nuc(i),2)=f(Nuc(i),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     > uz*auxd(3,2))
       f(Nuc(i),3)=f(Nuc(i),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c------------------------------------
       ccj=cc*tj
       ccx=ccj*txj
       ccy=ccj*tyj
       ccz=ccj*tzj
cc
       auxd(1,1)=ccx*sxs
       auxd(1,2)=ccy*sxs
       auxd(1,3)=ccz*sxs
c
       auxd(2,1)=ccx*sys
       auxd(2,2)=ccy*sys
       auxd(2,3)=ccz*sys
c
       auxd(3,1)=ccx*szs
       auxd(3,2)=ccy*szs
       auxd(3,3)=ccz*szs
c
       t1=ccj*ssz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       f(Nuc(j),1)=f(Nuc(j),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     >  uz*auxd(3,1))
       f(Nuc(j),2)=f(Nuc(j),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     >   uz*auxd(3,2))
       f(Nuc(j),3)=f(Nuc(j),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
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
      z2=2.D0*zij
      t0=a(i,ni)*a(j,nj)
      alf=t0/zij
c
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ssz2=ss/z2
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c l2: different p in the p shell ( x,y,z respectively)
c
      do 305 l1=1,3
        t1=Q(l1)-r(Nuc(i),l1)
        ps=t1*ss
        aux(1)=t1*sxs
        aux(2)=t1*sys
        aux(3)=t1*szs
c
        psz2=ps/z2
        ii=i+l1-1
c ii index , taking into account different components of the shell
c
        aux(l1)=aux(l1)+ssz2
c
        k=ii+((M2-j)*(j-1))/2
        cc=ccoef*RMM(k)
c 
c derivatives
c
       ti=2.0D0*a(i,ni)
       tj=2.0D0*a(j,nj)
       tx=Q(1)-r(Nuc(i),1)
       ty=Q(2)-r(Nuc(i),2)
       tz=Q(3)-r(Nuc(i),3)
c
       txj=Q(1)-r(Nuc(j),1)
       tyj=Q(2)-r(Nuc(j),2)
       tzj=Q(3)-r(Nuc(j),3)
c
       cci=cc*ti
       ccx=cci*tx
       ccy=cci*ty
       ccz=cci*tz
cc
c 1-3 derivatives of ux, 4-6 derivatives of uy, 7-9 of uz
       auxd(1,1)=ccx*aux(1)
       auxd(1,2)=ccy*aux(1)
       auxd(1,3)=ccz*aux(1)
c
       auxd(2,1)=ccx*aux(2)
       auxd(2,2)=ccy*aux(2)
       auxd(2,3)=ccz*aux(2)
c
       auxd(3,1)=ccx*aux(3)
       auxd(3,2)=ccy*aux(3)
       auxd(3,3)=ccz*aux(3)
c
       t1=cci*psz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=cci/z2-cc
        auxd(1,l1)=auxd(1,l1)+t1*sxs
        auxd(2,l1)=auxd(2,l1)+t1*sys
        auxd(3,l1)=auxd(3,l1)+t1*szs
c
       f(Nuc(i),1)=f(Nuc(i),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     > uz*auxd(3,1))
       f(Nuc(i),2)=f(Nuc(i),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     > uz*auxd(3,2))
       f(Nuc(i),3)=f(Nuc(i),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c------------------------------------
       ccj=cc*tj
       ccx=ccj*txj
       ccy=ccj*tyj
       ccz=ccj*tzj
cc
       auxd(1,1)=ccx*aux(1)
       auxd(1,2)=ccy*aux(1)
       auxd(1,3)=ccz*aux(1)
c
       auxd(2,1)=ccx*aux(2)
       auxd(2,2)=ccy*aux(2)
       auxd(2,3)=ccz*aux(2)
c
       auxd(3,1)=ccx*aux(3)
       auxd(3,2)=ccy*aux(3)
       auxd(3,3)=ccz*aux(3)
c
       t1=ccj*psz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=ccj/z2
        auxd(1,l1)=auxd(1,l1)+t1*sxs
        auxd(2,l1)=auxd(2,l1)+t1*sys
        auxd(3,l1)=auxd(3,l1)+t1*szs
c
c
       f(Nuc(j),1)=f(Nuc(j),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     >  uz*auxd(3,1))
       f(Nuc(j),2)=f(Nuc(j),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     >   uz*auxd(3,2))
       f(Nuc(j),3)=f(Nuc(j),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c
 305   continue
 300   continue
c
c------------------------------------------------------------
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
      t0=a(i,ni)*a(j,nj)
      alf=t0/zij
c
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ssz2=ss/z2
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c

      ccoef=c(i,ni)*c(j,nj)
c
      do 405 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
c
        aux(1)=t1*sxs
        aux(2)=t1*sys
        aux(3)=t1*szs
c
        aux(l1)=aux(l1)+ssz2
        ps=ss*t1
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
       do 405 l2=1,lij
c
       t1=Q(l2)-r(Nuc(j),l2)
       pp=t1*ps
       aux2(1)=t1*sxs
       aux2(2)=t1*sys
       aux2(3)=t1*szs

c
       aux2(l2)=aux2(l2)+ssz2
c
       aux1(1)=aux(1)*t1
       aux1(2)=aux(2)*t1
       aux1(3)=aux(3)*t1
       if (l1.eq.l2) then
        aux1(1)=aux1(1)+sxs/z2
        aux1(2)=aux1(2)+sys/z2
        aux1(3)=aux1(3)+szs/z2
        pp=pp+ss/z2
       endif
c
       aux1(l2)=aux1(l2)+ps/z2
       ii=i+l1-1
       jj=j+l2-1
c
       k=ii+((M2-jj)*(jj-1))/2
c
        cc=ccoef*RMM(k)
c
        ppz2=pp/z2
c derivatives
c
       ti=2.0D0*a(i,ni)
       tj=2.0D0*a(j,nj)
       tx=Q(1)-r(Nuc(i),1)
       ty=Q(2)-r(Nuc(i),2)
       tz=Q(3)-r(Nuc(i),3)
c
       txj=Q(1)-r(Nuc(j),1)
       tyj=Q(2)-r(Nuc(j),2)
       tzj=Q(3)-r(Nuc(j),3)
c
       cci=cc*ti
       ccx=cci*tx
       ccy=cci*ty
       ccz=cci*tz
cc
c 1-3 derivatives of ux, 4-6 derivatives of uy, 7-9 of uz
       auxd(1,1)=ccx*aux1(1)
       auxd(1,2)=ccy*aux1(1)
       auxd(1,3)=ccz*aux1(1)
c
       auxd(2,1)=ccx*aux1(2)
       auxd(2,2)=ccy*aux1(2)
       auxd(2,3)=ccz*aux1(2)
c
       auxd(3,1)=ccx*aux1(3)
       auxd(3,2)=ccy*aux1(3)
       auxd(3,3)=ccz*aux1(3)
c
       t1=cci*ppz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=cci/z2-cc
        auxd(1,l1)=auxd(1,l1)+t1*aux2(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux2(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux2(3)
c
       t1=cci/z2
        auxd(1,l2)=auxd(1,l2)+t1*aux(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux(3)
c
       f(Nuc(i),1)=f(Nuc(i),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     > uz*auxd(3,1))
       f(Nuc(i),2)=f(Nuc(i),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     > uz*auxd(3,2))
       f(Nuc(i),3)=f(Nuc(i),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c------------------------------------
       ccj=cc*tj
       ccx=ccj*txj
       ccy=ccj*tyj
       ccz=ccj*tzj
cc
       auxd(1,1)=ccx*aux1(1)
       auxd(1,2)=ccy*aux1(1)
       auxd(1,3)=ccz*aux1(1)
c
       auxd(2,1)=ccx*aux1(2)
       auxd(2,2)=ccy*aux1(2)
       auxd(2,3)=ccz*aux1(2)
c
       auxd(3,1)=ccx*aux1(3)
       auxd(3,2)=ccy*aux1(3)
       auxd(3,3)=ccz*aux1(3)
c
       t1=ccj*ppz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=ccj/z2
        auxd(1,l1)=auxd(1,l1)+t1*aux2(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux2(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux2(3)
c
       t1=ccj/z2-cc
        auxd(1,l2)=auxd(1,l2)+t1*aux(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux(3)
c
c
c
       f(Nuc(j),1)=f(Nuc(j),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     >  uz*auxd(3,1))
       f(Nuc(j),2)=f(Nuc(j),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     >   uz*auxd(3,2))
       f(Nuc(j),3)=f(Nuc(j),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c
c 
 405  continue
c
 400  continue
c
c------------------------------------------------------------------
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
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       tx=Q(1)-r(Nuc(i),1)
       ty=Q(2)-r(Nuc(i),2)
       tz=Q(3)-r(Nuc(i),3)
c
       txj=Q(1)-r(Nuc(j),1)
       tyj=Q(2)-r(Nuc(j),2)
       tzj=Q(3)-r(Nuc(j),3)
       ti=2.0D0*a(i,ni)
       tj=2.0D0*a(j,nj)
c
c
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      t3a=ss/z2
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c
      ccoef=c(i,ni)*c(j,nj)
c
c
      do 505 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       ps=ss*t1
        aux(1)=t1*sxs
        aux(2)=t1*sys
        aux(3)=t1*szs
c
        aux(l1)=aux(l1)+t3a
c
      do 505 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       ds=t1*ps
       aux1(1)=aux(1)*t1
       aux1(2)=aux(2)*t1
       aux1(3)=aux(3)*t1
c
        aux2(1)=t1*sxs
        aux2(2)=t1*sys
        aux2(3)=t1*szs
c
        aux2(l2)=aux2(l2)+t3a
       
       f1=1.D0
       if (l1.eq.l2) then
        aux1(1)=aux1(1)+sxs/z2
        aux1(2)=aux1(2)+sys/z2
        aux1(3)=aux1(3)+szs/z2
        ds=ds+t3a
        f1=sq3
       endif
c
       aux1(l2)=aux1(l2)+ps/z2
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       ii=i+l12-1
c
       k=ii+((M2-j)*(j-1))/2
       cc = ccoef/f1*RMM(k)
c
c derivatives
c
       cci=cc*ti
       ccx=cci*tx
       ccy=cci*ty
       ccz=cci*tz
cc
       auxd(1,1)=ccx*aux1(1)
       auxd(1,2)=ccy*aux1(1)
       auxd(1,3)=ccz*aux1(1)
c
       auxd(2,1)=ccx*aux1(2)
       auxd(2,2)=ccy*aux1(2)
       auxd(2,3)=ccz*aux1(2)
c
       auxd(3,1)=ccx*aux1(3)
       auxd(3,2)=ccy*aux1(3)
       auxd(3,3)=ccz*aux1(3)
c
       dsz2=ds/z2
       t1=cci*dsz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=cci/z2-cc
        auxd(1,l1)=auxd(1,l1)+t1*aux2(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux2(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux2(3)
c
        auxd(1,l2)=auxd(1,l2)+t1*aux(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux(3)
c
       f(Nuc(i),1)=f(Nuc(i),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     > uz*auxd(3,1))
       f(Nuc(i),2)=f(Nuc(i),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     > uz*auxd(3,2))
       f(Nuc(i),3)=f(Nuc(i),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c------------------------------------
       ccj=cc*tj
       ccx=ccj*txj
       ccy=ccj*tyj
       ccz=ccj*tzj
cc
       auxd(1,1)=ccx*aux1(1)
       auxd(1,2)=ccy*aux1(1)
       auxd(1,3)=ccz*aux1(1)
c
       auxd(2,1)=ccx*aux1(2)
       auxd(2,2)=ccy*aux1(2)
       auxd(2,3)=ccz*aux1(2)
c
       auxd(3,1)=ccx*aux1(3)
       auxd(3,2)=ccy*aux1(3)
       auxd(3,3)=ccz*aux1(3)
c
       t1=ccj*dsz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=ccj/z2
        auxd(1,l1)=auxd(1,l1)+t1*aux2(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux2(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux2(3)
c
        auxd(1,l2)=auxd(1,l2)+t1*aux(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux(3)
c
c
       f(Nuc(j),1)=f(Nuc(j),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     >  uz*auxd(3,1))
       f(Nuc(j),2)=f(Nuc(j),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     >   uz*auxd(3,2))
       f(Nuc(j),3)=f(Nuc(j),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c
 505  continue
 500  continue
c
c-----------------------------------------------------------------------
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
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       tx=Q(1)-r(Nuc(i),1)
       ty=Q(2)-r(Nuc(i),2)
       tz=Q(3)-r(Nuc(i),3)
c
       txj=Q(1)-r(Nuc(j),1)
       tyj=Q(2)-r(Nuc(j),2)
       tzj=Q(3)-r(Nuc(j),3)
c
       ti=2.0D0*a(i,ni)
       tj=2.0D0*a(j,nj)
c
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ssz2=ss/z2
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
      sxz2=Q(1)*ssz2
      syz2=Q(2)*ssz2
      szz2=Q(3)*ssz2
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 605 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
c aux : (pi|r|s)
        aux(1)=t1*sxs
        aux(2)=t1*sys
        aux(3)=t1*szs
        aux(l1)=aux(l1)+ssz2
c
      do 605 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       dijs=t1*pis
c aux1 : (pj|r|s)
        aux1(1)=t1*sxs
        aux1(2)=t1*sys
        aux1(3)=t1*szs
        aux1(l2)=aux1(l2)+ssz2
c
c aux2 : (dij|r|s)
        aux2(1)=t1*aux(1)
        aux2(2)=t1*aux(2)
        aux2(3)=t1*aux(3)
        aux2(l2)=aux2(l2)+pis/z2
       f1=1.D0
c
       if (l1.eq.l2) then
        f1=sq3
        aux2(1)=aux2(1)+sxz2
        aux2(2)=aux2(2)+syz2
        aux2(3)=aux2(3)+szz2
        dijs=dijs+ssz2
       endif
c
      do 605 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       dp=t1*dijs
c
       aux3(1)=aux2(1)*t1
       aux3(2)=aux2(2)*t1
       aux3(3)=aux2(3)*t1
c
       aux4(1)=aux(1)*t1
       aux4(2)=aux(2)*t1
       aux4(3)=aux(3)*t1
c
       aux5(1)=aux1(1)*t1
       aux5(2)=aux1(2)*t1
       aux5(3)=aux1(3)*t1
c


       if (l1.eq.l3) then
        aux3(1)=aux3(1)+aux1(1)/z2
        aux3(2)=aux3(2)+aux1(2)/z2
        aux3(3)=aux3(3)+aux1(3)/z2
        aux4(1)=aux4(1)+sxz2
        aux4(2)=aux4(2)+syz2
        aux4(3)=aux4(3)+szz2
        dp=dp+pjs/z2
       endif
c
       if (l2.eq.l3) then
        aux3(1)=aux3(1)+aux(1)/z2
        aux3(2)=aux3(2)+aux(2)/z2
        aux3(3)=aux3(3)+aux(3)/z2
        aux5(1)=aux5(1)+sxz2
        aux5(2)=aux5(2)+syz2
        aux5(3)=aux5(3)+szz2
        dp=dp+pis/z2
       endif
c
       aux3(l3)=aux3(l3)+dijs/z2
       aux4(l3)=aux4(l3)+pis/z2
       aux5(l3)=aux5(l3)+pjs/z2
c
       dpz2=dp/z2
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1*RMM(k)
c
c derivatives
c
       cci=cc*ti
       ccx=cci*tx
       ccy=cci*ty
       ccz=cci*tz
cc
       auxd(1,1)=ccx*aux3(1)
       auxd(1,2)=ccy*aux3(1)
       auxd(1,3)=ccz*aux3(1)
c
       auxd(2,1)=ccx*aux3(2)
       auxd(2,2)=ccy*aux3(2)
       auxd(2,3)=ccz*aux3(2)
c
       auxd(3,1)=ccx*aux3(3)
       auxd(3,2)=ccy*aux3(3)
       auxd(3,3)=ccz*aux3(3)
c
       t1=cci*dpz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=cci/z2-cc
        auxd(1,l1)=auxd(1,l1)+t1*aux5(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux5(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux5(3)
c
        auxd(1,l2)=auxd(1,l2)+t1*aux4(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux4(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux4(3)
c
       t1=cci/z2
        auxd(1,l3)=auxd(1,l3)+t1*aux2(1)
        auxd(2,l3)=auxd(2,l3)+t1*aux2(2)
        auxd(3,l3)=auxd(3,l3)+t1*aux2(3)
c
       f(Nuc(i),1)=f(Nuc(i),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     > uz*auxd(3,1))
       f(Nuc(i),2)=f(Nuc(i),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     > uz*auxd(3,2))
       f(Nuc(i),3)=f(Nuc(i),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c------------------------------------
       ccj=cc*tj
       ccx=ccj*txj
       ccy=ccj*tyj
       ccz=ccj*tzj
cc
       auxd(1,1)=ccx*aux3(1)
       auxd(1,2)=ccy*aux3(1)
       auxd(1,3)=ccz*aux3(1)
c
       auxd(2,1)=ccx*aux3(2)
       auxd(2,2)=ccy*aux3(2)
       auxd(2,3)=ccz*aux3(2)
c
       auxd(3,1)=ccx*aux3(3)
       auxd(3,2)=ccy*aux3(3)
       auxd(3,3)=ccz*aux3(3)
c
       t1=ccj*dpz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=ccj/z2
        auxd(1,l1)=auxd(1,l1)+t1*aux5(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux5(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux5(3)
c
        auxd(1,l2)=auxd(1,l2)+t1*aux4(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux4(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux4(3)
c
        t1=ccj/z2-cc
        auxd(1,l3)=auxd(1,l3)+t1*aux2(1)
        auxd(2,l3)=auxd(2,l3)+t1*aux2(2)
        auxd(3,l3)=auxd(3,l3)+t1*aux2(3)
c
       f(Nuc(j),1)=f(Nuc(j),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     >  uz*auxd(3,1))
       f(Nuc(j),2)=f(Nuc(j),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     >   uz*auxd(3,2))
       f(Nuc(j),3)=f(Nuc(j),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
 605  continue
 600  continue
c
c ----------------------------------------------------------------
c (d|d) case
c
      do 700 i=ns+np+1,M,6
      do 700 j=ns+np+1,M,6
c
      dd=d(Nuc(i),Nuc(j))
c
c
      do 700 ni=1,ncont(i)
      do 700 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       tx=Q(1)-r(Nuc(i),1)
       ty=Q(2)-r(Nuc(i),2)
       tz=Q(3)-r(Nuc(i),3)
c
       txj=Q(1)-r(Nuc(j),1)
       tyj=Q(2)-r(Nuc(j),2)
       tzj=Q(3)-r(Nuc(j),3)
c
       ti=2.0D0*a(i,ni)
       tj=2.0D0*a(j,nj)
c
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      ssz2=ss/z2
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c
      sxz2=Q(1)*ssz2
      syz2=Q(2)*ssz2
      szz2=Q(3)*ssz2
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 705 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
c aux : (pi|r|s)
        aux(1)=t1*sxs
        aux(2)=t1*sys
        aux(3)=t1*szs
        aux(l1)=aux(l1)+ssz2
c
      do 705 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       dijs=t1*pis
c aux1 : (pj|r|s)
        aux1(1)=t1*sxs
        aux1(2)=t1*sys
        aux1(3)=t1*szs
        aux1(l2)=aux1(l2)+ssz2
c
c aux2 : (dij|r|s)
        aux2(1)=t1*aux(1)
        aux2(2)=t1*aux(2)
        aux2(3)=t1*aux(3)
        aux2(l2)=aux2(l2)+pis/z2
       f1=1.D0
c
       if (l1.eq.l2) then
        f1=sq3
        aux2(1)=aux2(1)+sxz2
        aux2(2)=aux2(2)+syz2
        aux2(3)=aux2(3)+szz2
        dijs=dijs+ssz2
       endif
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 705 l3=1,lij
c
       t1=Q(l3)-r(Nuc(j),l3)
       dp=t1*dijs
       pipk=t1*pis
       pjpk=t1*pjs
c
c aux3 : (dij|r|pk)
       aux3(1)=aux2(1)*t1
       aux3(2)=aux2(2)*t1
       aux3(3)=aux2(3)*t1
c aux4 : (pi|r|pk)
       aux4(1)=aux(1)*t1
       aux4(2)=aux(2)*t1
       aux4(3)=aux(3)*t1
c aux5 : (pj|r|pk)
       aux5(1)=aux1(1)*t1
       aux5(2)=aux1(2)*t1
       aux5(3)=aux1(3)*t1
c aux10 (s|r|pk)
       aux10(1)=sxs*t1
       aux10(2)=sys*t1
       aux10(3)=szs*t1
c
       if (l1.eq.l3) then
        dp=dp+pjs/z2
        pipk=pipk+ss/z2
        aux3(1)=aux3(1)+aux1(1)/z2
        aux3(2)=aux3(2)+aux1(2)/z2
        aux3(3)=aux3(3)+aux1(3)/z2
        aux4(1)=aux4(1)+sxz2
        aux4(2)=aux4(2)+syz2
        aux4(3)=aux4(3)+szz2
       endif
c
       if (l2.eq.l3) then
        dp=dp+pis/z2
        pjpk=pjpk+ss/z2
        aux3(1)=aux3(1)+aux(1)/z2
        aux3(2)=aux3(2)+aux(2)/z2
        aux3(3)=aux3(3)+aux(3)/z2
        aux5(1)=aux5(1)+sxz2
        aux5(2)=aux5(2)+syz2
        aux5(3)=aux5(3)+szz2
       endif
c
       aux3(l3)=aux3(l3)+dijs/z2
       aux4(l3)=aux4(l3)+pis/z2
       aux5(l3)=aux5(l3)+pjs/z2
c
       aux10(l3)=aux10(l3)+ssz2
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,l1*(l1-1)/2-l3*(l3-1)/2+l2)
      endif
c
       do 705 l4=1,lk
       t1=Q(l4)-r(Nuc(j),l4)
       dsd=t1*dp
c aux3 : used here for (d|r|d)
       aux6(1)=aux3(1)*t1
       aux6(2)=aux3(2)*t1
       aux6(3)=aux3(3)*t1
c aux7 (pi|r|d)
       aux7(1)=aux4(1)*t1
       aux7(2)=aux4(2)*t1
       aux7(3)=aux4(3)*t1
c aux8 (pj|r|d)
       aux8(1)=aux5(1)*t1
       aux8(2)=aux5(2)*t1
       aux8(3)=aux5(3)*t1
c aux9 (d|r|pl)
       aux9(1)=aux2(1)*t1
       aux9(2)=aux2(2)*t1
       aux9(3)=aux2(3)*t1
c
       if (l1.eq.l4) then
        dsd=dsd+pjpk/z2
        aux6(1)=aux6(1)+aux5(1)/z2
        aux6(2)=aux6(2)+aux5(2)/z2
        aux6(3)=aux6(3)+aux5(3)/z2
c
        aux9(1)=aux9(1)+aux1(1)/z2
        aux9(2)=aux9(2)+aux1(2)/z2
        aux9(3)=aux9(3)+aux1(3)/z2
c
        aux7(1)=aux7(1)+aux10(1)/z2
        aux7(2)=aux7(2)+aux10(2)/z2
        aux7(3)=aux7(3)+aux10(3)/z2
       endif
c
       if (l2.eq.l4) then
        dsd=dsd+pipk/z2
        aux6(1)=aux6(1)+aux4(1)/z2
        aux6(2)=aux6(2)+aux4(2)/z2
        aux6(3)=aux6(3)+aux4(3)/z2
c
        aux9(1)=aux9(1)+aux(1)/z2
        aux9(2)=aux9(2)+aux(2)/z2
        aux9(3)=aux9(3)+aux(3)/z2
c
        aux8(1)=aux8(1)+aux10(1)/z2
        aux8(2)=aux8(2)+aux10(2)/z2
        aux8(3)=aux8(3)+aux10(3)/z2
       endif
c
       f2=1.D0
       if (l3.eq.l4) then
       f2=sq3
        dsd=dsd+dijs/z2
        aux6(1)=aux6(1)+aux2(1)/z2
        aux6(2)=aux6(2)+aux2(2)/z2
        aux6(3)=aux6(3)+aux2(3)/z2
c
        aux7(1)=aux7(1)+aux(1)/z2
        aux7(2)=aux7(2)+aux(2)/z2
        aux7(3)=aux7(3)+aux(3)/z2
c
        aux8(1)=aux8(1)+aux1(1)/z2
        aux8(2)=aux8(2)+aux1(2)/z2
        aux8(3)=aux8(3)+aux1(3)/z2
       endif
c
       dsdz2=dsd/z2
       aux6(l4)=aux6(l4)+dp/z2
       aux7(l4)=aux7(l4)+pipk/z2
       aux8(l4)=aux8(l4)+pjpk/z2
       aux9(l4)=aux9(l4)+dijs/z2
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/(f1*f2)*RMM(k)
c
c derivatives -------------------------------------------------
c
       cci=cc*ti
       ccx=cci*tx
       ccy=cci*ty
       ccz=cci*tz
cc
       auxd(1,1)=ccx*aux6(1)
       auxd(1,2)=ccy*aux6(1)
       auxd(1,3)=ccz*aux6(1)
c
       auxd(2,1)=ccx*aux6(2)
       auxd(2,2)=ccy*aux6(2)
       auxd(2,3)=ccz*aux6(2)
c
       auxd(3,1)=ccx*aux6(3)
       auxd(3,2)=ccy*aux6(3)
       auxd(3,3)=ccz*aux6(3)
c
       t1=cci*dsdz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=cci/z2-cc
        auxd(1,l1)=auxd(1,l1)+t1*aux8(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux8(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux8(3)
c
        auxd(1,l2)=auxd(1,l2)+t1*aux7(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux7(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux7(3)
c
       t1=cci/z2
        auxd(1,l3)=auxd(1,l3)+t1*aux9(1)
        auxd(2,l3)=auxd(2,l3)+t1*aux9(2)
        auxd(3,l3)=auxd(3,l3)+t1*aux9(3)
c
        auxd(1,l4)=auxd(1,l4)+t1*aux3(1)
        auxd(2,l4)=auxd(2,l4)+t1*aux3(2)
        auxd(3,l4)=auxd(3,l4)+t1*aux3(3)
c
       f(Nuc(i),1)=f(Nuc(i),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     > uz*auxd(3,1))
       f(Nuc(i),2)=f(Nuc(i),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     > uz*auxd(3,2))
       f(Nuc(i),3)=f(Nuc(i),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
c------------------------------------
       ccj=cc*tj
       ccx=ccj*txj
       ccy=ccj*tyj
       ccz=ccj*tzj
cc
       auxd(1,1)=ccx*aux6(1)
       auxd(1,2)=ccy*aux6(1)
       auxd(1,3)=ccz*aux6(1)
c
       auxd(2,1)=ccx*aux6(2)
       auxd(2,2)=ccy*aux6(2)
       auxd(2,3)=ccz*aux6(2)
c
       auxd(3,1)=ccx*aux6(3)
       auxd(3,2)=ccy*aux6(3)
       auxd(3,3)=ccz*aux6(3)
c
       t1=ccj*dsdz2
       auxd(1,1)=auxd(1,1)+t1
       auxd(2,2)=auxd(2,2)+t1
       auxd(3,3)=auxd(3,3)+t1
c
       t1=ccj/z2
        auxd(1,l1)=auxd(1,l1)+t1*aux8(1)
        auxd(2,l1)=auxd(2,l1)+t1*aux8(2)
        auxd(3,l1)=auxd(3,l1)+t1*aux8(3)
c
        auxd(1,l2)=auxd(1,l2)+t1*aux7(1)
        auxd(2,l2)=auxd(2,l2)+t1*aux7(2)
        auxd(3,l2)=auxd(3,l2)+t1*aux7(3)
c
        t1=ccj/z2-cc
        auxd(1,l3)=auxd(1,l3)+t1*aux9(1)
        auxd(2,l3)=auxd(2,l3)+t1*aux9(2)
        auxd(3,l3)=auxd(3,l3)+t1*aux9(3)
c
        auxd(1,l4)=auxd(1,l4)+t1*aux3(1)
        auxd(2,l4)=auxd(2,l4)+t1*aux3(2)
        auxd(3,l4)=auxd(3,l4)+t1*aux3(3)
c
       f(Nuc(j),1)=f(Nuc(j),1)+g1*(ux*auxd(1,1)+uy*auxd(2,1)+
     >  uz*auxd(3,1))
       f(Nuc(j),2)=f(Nuc(j),2)+g1*(ux*auxd(1,2)+uy*auxd(2,2)+
     >   uz*auxd(3,2))
       f(Nuc(j),3)=f(Nuc(j),3)+g1*(ux*auxd(1,3)+uy*auxd(2,3)+
     > uz*auxd(3,3))
c
 705  continue
 700  continue
c--------------------------------------------------------------

c ----------------------------------------------------------------

c Nuclear contributions
      do n=1,natom
      f(n,1)=f(n,1)-Iz(n)*g*ux
      f(n,2)=f(n,2)-Iz(n)*g*uy
      f(n,3)=f(n,3)-Iz(n)*g*uz
c
      enddo
c
      return
      end
