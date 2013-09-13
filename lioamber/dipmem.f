c-------------------------------------------------------------------
c
c Subroutine for calculating dipole moment
c ONLY FOR NEUTRAL SYSTEMS
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
      subroutine dipmem(ux,uy,uz,elmu)
      use garcha_mod
c
      dimension aux(3),aux1(3),aux2(3),aux3(3),aux4(3),aux5(3)
      dimension aux6(3)
      dimension Q(3)
c
      dimension elmu(3,M*(M+1)/2)
c
      elmu=0
c
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
c      do 50 i=1,natom
c      do 50 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 50   continue
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
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
      ccoef=c(i,ni)*c(j,nj)
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
c
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
      k=i+((M2-j)*(j-1))/2
      elmu(1,k)=elmu(1,k) + ccoef*sxs
      elmu(2,k)=elmu(2,k) + ccoef*sys
      elmu(3,k)=elmu(3,k) + ccoef*szs
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
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
        
       do jjj=1,3
       elmu(jjj,k)=elmu(jjj,k) + ccoef*aux(jjj)
       enddo

c       cc=ccoef*RMM(k)
c       ux=ux+cc*aux(1)
c       uy=uy+cc*aux(2)
c       uz=uz+cc*aux(3)
c
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
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
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
       aux1(1)=aux(1)*t1
       aux1(2)=aux(2)*t1
       aux1(3)=aux(3)*t1
c
       if (l1.eq.l2) then
        aux1(1)=aux1(1)+sxs/z2
        aux1(2)=aux1(2)+sys/z2
        aux1(3)=aux1(3)+szs/z2
       endif
c
       aux1(l2)=aux1(l2)+ps/z2
       ii=i+l1-1
       jj=j+l2-1
c
c      eliminated
       if(ii.ge.jj) then
       k=ii+((M2-jj)*(jj-1))/2
       
       do jjj=1,3
       elmu(jjj,k)=elmu(jjj,k) + ccoef*aux1(jjj)
       enddo


c       cc=RMM(k)*ccoef
c       ux=ux+cc*aux1(1)
c       uy=uy+cc*aux1(2)
c       uz=uz+cc*aux1(3)
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
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

       do jjj=1,3
       elmu(jjj,k)=elmu(jjj,k) + ccoef*aux1(jjj)/f1
       enddo
c       cc = RMM(k)*ccoef/f1
c
c       ux=ux+cc*aux1(1)
c       uy=uy+cc*aux1(2)
c       uz=uz+cc*aux1(3)
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
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

       do jjj=1,3
       elmu(jjj,k)=elmu(jjj,k) + ccoef*aux3(jjj)/f1
       enddo
c       cc=ccoef/f1*RMM(k)
c
c       ux=ux+cc*aux3(1)
c       uy=uy+cc*aux3(2)
c       uz=uz+cc*aux3(3)
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
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      alf2=2.D0*alf
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

       do jjj=1,3
       elmu(jjj,k)=elmu(jjj,k) + ccoef*aux6(jjj)/(f1*f2)
       enddo
c       cc=ccoef/(f1*f2)*RMM(k)
c
c       ux=ux+cc*aux6(1)
c       uy=uy+cc*aux6(2)
c       uz=uz+cc*aux6(3)
       endif
 705  continue
 700  continue
c--------------------------------------------------------------
c--------------------------------------------------------------
c
 999  continue
      uxa=0.0D0
      uya=0.0D0
      uza=0.0D0
       Qc=0.0D0
      do 99 i=1,natom
       Qc=Qc+Iz(i)
       uxa=uxa+Iz(i)*r(i,1)
       uya=uya+Iz(i)*r(i,2)
       uza=uza+Iz(i)*r(i,3)
 99   continue
c
      Qc=Qc-Nel
c
      if (sol) then
       do n=1,Nsol
       do k=1,natsol
       Qc=Qc+pc(k)
       enddo
       enddo
      endif

        




c fac : for charged species dipole moment depends on the origin
c definition - Using the factor, it is defined with respect to the
c center of charge (important in Reaction Field calculations)
c For neutral systems it doesn't matter
c
      fac=(Qc+Nel)/Nel
c test
c     fac=1
c
      ux=uxa-ux*fac
      uy=uya-uy*fac
      uz=uza-uz*fac
c
c
c Contributions from the solvent (calculated from point charges)
c
c     if (sol) then
c     n1=natom
c     do n=1,Nsol
c      do k=1,natsol
c      n1=n1+1
c
c      ux=ux +pc(k)*r(n1,1)
c      uy=uy+ pc(k)*r(n1,2)
c      uz=uz+ pc(k)*r(n1,3)
c     enddo
c     enddo
c
c     endif
c
      ux=ux*2.54D0
      uy=uy*2.54D0
      uz=uz*2.54D0
c
c u in Debyes
c     u=sqrt(ux**2+uy**2+uz**2)
c
c     write(*,*)
c     write(*,*) 'DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
c     write(*,900) ux,uy,uz,u
c     write(*,*)
c u in Debyes
c
c900  format(3(F10.4,2x),2x,F10.4)
      return
      end
