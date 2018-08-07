c-------------------------------------------------------------------
c
c Subroutine for evaluating the reaction field contribution to
c the one electron matrix elements
c
c Dipole moment should be given in input (3 components)
c
c integrals evaluated using Obara Saika method
c 19-1-93
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
      module subm_intfld; contains
      subroutine intfld(g,ux,uy,uz)
      use garcha_mod, only: RMM, a, c, d, r, nuc, ncont, nshell
     >                    , Iz, OPEN, NORM, M, Md, natom, pi32
c
c
      implicit none
      real*8, intent(in) :: g
      real*8, intent(in) :: ux, uy, uz
c aux . things
      real*8  :: aux(3) , aux1(3), aux2(3), aux3(3), aux4(3)
      real*8  :: aux5(3), aux6(3), Q(3)

! Implicits:
      integer :: MM, MM2, MMd, Md2
      integer :: M1, M2, M3, M5
      integer :: l1, l2, l3, l4, l12, l34
      integer :: ns, np, nd, ni, nj, i, j, k, ii, jj

      real*8  :: sq3, term, alf, ccoef, cc, f1, f2
      real*8  :: ss, ps, pp, pis, pjs, dp, dd, dijs
      real*8  :: sxs, sys, szs, uxa, uya, uza
      real*8  :: t0, t1, ti, tj, z2, zij
c
c
c pointers
c
      MM=M*(M+1)/2
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M
c first P
      M1=1
c now Pnew
c     M3=M1+MM
c now S, F also uses the same position after S was used
c Fock matrix at this point
c     M5=M3+MM
c first P
      M1=1
c now F alpha
      M3=M1+MM
c now S, F beta also uses the same position after S was used
      M5=M3+MM
c-----------------------------------------------------
c Nuclear component of dipole moment
      uxa=0.0D0
      uya=0.0D0
      uza=0.0D0
      do 99 i=1,natom
       uxa=uxa+Iz(i)*r(i,1)
       uya=uya+Iz(i)*r(i,2)
       uza=uza+Iz(i)*r(i,3)
 99   continue
c
c
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
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
c-----------------------------------------
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
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
      k=i+((M2-j)*(j-1))/2
c
c
      term=ccoef*(sxs*ux+sys*uy+szs*uz)
      RMM(M5+k-1)=RMM(M5+k-1)+g*term
      if (OPEN) then
       RMM(M3+k-1)=RMM(M3+k-1)+g*term
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
        ii=i+l1-1
c ii index , taking into account different components of the shell
c
        aux(l1)=aux(l1)+ss/z2
c
        k=ii+((M2-j)*(j-1))/2
c
      term=ccoef*(aux(1)*ux+aux(2)*uy+aux(3)*uz)
      RMM(M5+k-1)=RMM(M5+k-1)+g*term
c
      if (OPEN) then
       RMM(M3+k-1)=RMM(M3+k-1)+g*term
      endif
 305   continue
 300   continue
c
c------------------------------------------------------------
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
      t0=a(i,ni)*a(j,nj)
      alf=t0/zij
c
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
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
        aux(l1)=aux(l1)+ss/z2
        ps=ss*t1
c
       do 405 l2=1,3
c
       t1=Q(l2)-r(Nuc(j),l2)
       pp=t1*ps
c
       aux1(1)=aux(1)*t1
       aux1(2)=aux(2)*t1
       aux1(3)=aux(3)*t1
c
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
c      eliminated
       if(ii.ge.jj) then
       k=ii+((M2-jj)*(jj-1))/2
c
c
      term=ccoef*(aux1(1)*ux+aux1(2)*uy+aux1(3)*uz)
      RMM(M5+k-1)=RMM(M5+k-1)+g*term
c
      if (OPEN) then
       RMM(M3+k-1)=RMM(M3+k-1)+g*term
      endif
c
       endif
 405  continue
c
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
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 505 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       ps=ss*t1
        aux(1)=t1*sxs
        aux(2)=t1*sys
        aux(3)=t1*szs
c
        aux(l1)=aux(l1)+ss/z2
c
      do 505 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       aux1(1)=aux(1)*t1
       aux1(2)=aux(2)*t1
       aux1(3)=aux(3)*t1
c
       f1=1.D0
       if (l1.eq.l2) then
        aux1(1)=aux1(1)+sxs/z2
        aux1(2)=aux1(2)+sys/z2
        aux1(3)=aux1(3)+szs/z2
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
       cc = ccoef/f1
c
      term=cc*(aux1(1)*ux+aux1(2)*uy+aux1(3)*uz)
      RMM(M5+k-1)=RMM(M5+k-1)+g*term
      if (OPEN) then
       RMM(M3+k-1)=RMM(M3+k-1)+g*term
      endif
c
c
 505  continue
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
      ti=a(i,ni)/zij
      tj=a(j,nj)/zij
      Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
      Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
      Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
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
        aux(l1)=aux(l1)+ss/z2
c
      do 605 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       dijs=t1*pis
c aux1 : (pj|r|s)
        aux1(1)=t1*sxs
        aux1(2)=t1*sys
        aux1(3)=t1*szs
        aux1(l2)=aux1(l2)+ss/z2
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
        aux2(1)=aux2(1)+sxs/z2
        aux2(2)=aux2(2)+sys/z2
        aux2(3)=aux2(3)+szs/z2
        dijs=dijs+ss/z2
       endif
c
      do 605 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
c
       aux3(1)=aux2(1)*t1
       aux3(2)=aux2(2)*t1
       aux3(3)=aux2(3)*t1
c

       if (l1.eq.l3) then
        aux3(1)=aux3(1)+aux1(1)/z2
        aux3(2)=aux3(2)+aux1(2)/z2
        aux3(3)=aux3(3)+aux1(3)/z2
       endif
c
       if (l2.eq.l3) then
        aux3(1)=aux3(1)+aux(1)/z2
        aux3(2)=aux3(2)+aux(2)/z2
        aux3(3)=aux3(3)+aux(3)/z2
       endif
c
       aux3(l3)=aux3(l3)+dijs/z2
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
c
c
      term=cc*(aux3(1)*ux+aux3(2)*uy+aux3(3)*uz)
      RMM(M5+k-1)=RMM(M5+k-1)+g*term
      if (OPEN) then
       RMM(M3+k-1)=RMM(M3+k-1)+g*term
      endif
c
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
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
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
        aux(l1)=aux(l1)+ss/z2
c
      do 705 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       pjs=ss*t1
       dijs=t1*pis
c aux1 : (pj|r|s)
        aux1(1)=t1*sxs
        aux1(2)=t1*sys
        aux1(3)=t1*szs
        aux1(l2)=aux1(l2)+ss/z2
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
        aux2(1)=aux2(1)+sxs/z2
        aux2(2)=aux2(2)+sys/z2
        aux2(3)=aux2(3)+szs/z2
        dijs=dijs+ss/z2
       endif
c
      do 705 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       dp=t1*dijs
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
c
       if (l1.eq.l3) then
        dp=dp+pjs/z2
        aux3(1)=aux3(1)+aux1(1)/z2
        aux3(2)=aux3(2)+aux1(2)/z2
        aux3(3)=aux3(3)+aux1(3)/z2
        aux4(1)=aux4(1)+sxs/z2
        aux4(2)=aux4(2)+sys/z2
        aux4(3)=aux4(3)+szs/z2
       endif
c
       if (l2.eq.l3) then
        dp=dp+pis/z2
        aux3(1)=aux3(1)+aux(1)/z2
        aux3(2)=aux3(2)+aux(2)/z2
        aux3(3)=aux3(3)+aux(3)/z2
        aux5(1)=aux5(1)+sxs/z2
        aux5(2)=aux5(2)+sys/z2
        aux5(3)=aux5(3)+szs/z2
       endif
c
       aux3(l3)=aux3(l3)+dijs/z2
       aux4(l3)=aux4(l3)+pis/z2
       aux5(l3)=aux5(l3)+pjs/z2
c
       do 705 l4=1,l3
       t1=Q(l4)-r(Nuc(j),l4)
c aux3 : used here for (d|r|d)
       aux6(1)=aux3(1)*t1
       aux6(2)=aux3(2)*t1
       aux6(3)=aux3(3)*t1
c
       if (l1.eq.l4) then
        aux6(1)=aux6(1)+aux5(1)/z2
        aux6(2)=aux6(2)+aux5(2)/z2
        aux6(3)=aux6(3)+aux5(3)/z2
       endif
c
       if (l2.eq.l4) then
        aux6(1)=aux6(1)+aux4(1)/z2
        aux6(2)=aux6(2)+aux4(2)/z2
        aux6(3)=aux6(3)+aux4(3)/z2
       endif
c
       f2=1.D0
       if (l3.eq.l4) then
       f2=sq3
        aux6(1)=aux6(1)+aux2(1)/z2
        aux6(2)=aux6(2)+aux2(2)/z2
        aux6(3)=aux6(3)+aux2(3)/z2
       endif
c
       aux6(l4)=aux6(l4)+dp/z2
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       if (ii.ge.jj) then
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/(f1*f2)
c
c
      term=cc*(aux6(1)*ux+aux6(2)*uy+aux6(3)*uz)
      RMM(M5+k-1)=RMM(M5+k-1)+g*term
      if (OPEN) then
       RMM(M3+k-1)=RMM(M3+k-1)+g*term
      endif
c
       endif
 705  continue
 700  continue
c--------------------------------------------------------------
c--------------------------------------------------------------
c
      return
      end subroutine
      end module subm_intfld
