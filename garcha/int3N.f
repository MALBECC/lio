c-------------------------------------------------------------------
c Integrals subroutine -Third part
c 2 e integrals, 3 index : wavefunction and density fitting functions
c All of them are calculated
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
c Input : G ,F,  standard basis and density basis
c F comes, computed the 1 electron part, and here the
c Coulomb part is added, without storing the integrals
c Output: F updated with Coulomb part, also Coulomb energy
c F also updated with exchange correlation part, also energy
c is updated
c this subroutine calls the fitting for exchange-correlation
c-----------------------------------------------------------------
      subroutine int3N(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,E2,Ex,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
c
      implicit real*8 (a-h,o-z)
      logical NORM,dens,OPEN,SVD,ATRHO,integ
      logical VCINP,DIRECT,EXTR,SHFT,write
      integer nopt,iconst,igrid,igrid2
      INCLUDE 'param'
      parameter(pi52=34.9868366552497108D0,pi=3.14159265358979312D0)
      parameter (rmax=25.D0)
      dimension r(nt,3),nshelld(0:3),nshell(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(natom)
c
      dimension Q(3),W(3),Rc(ngd),af(ngd),FF(ngd),P(ngd)
      dimension d(ntq,ntq),Jx(ng),Ll(3)
      dimension RMM(*),X(Md,Md)
c scratch space
c
c auxiliars
      dimension B(ngd,3),aux(ngd)
      COMMON /TABLE/ STR(880,0:21)
c     common /HF/ nopt,OPEN,NMAX,NCO,ATRHO,VCINP,DIRECT,
c    >             IDAMP,EXTR,SHFT,SHI,GOLD,told,write,Nunp
      common /Sys/ SVD,iconst
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
      common /Nc/ Ndens
      common /coef/ af
      common /coef2/ B
      common /intg1/ e_(50,3),wang(50)
      common /intg2/ e_2(116,3),wang2(116),Nr(0:54),e3(194,3),wang3(194)      
c
c------------------------------------------------------------------
c now 16 loops for all combinations, first 2 correspond to 
c wavefunction basis, the third correspond to the density fit
c Rc(k) is constructed adding t(i,j,k)*P(i,j)
c cf(k) , variationally obtained fitting coefficient, is
c obtained by adding R(i)*G-1(i,k)
c if the t(i,j,k) were not stored, then in order to evaluate
c the corresponding part of the Fock matrix, they should be
c calculated again.
c V(i,j) obtained by adding af(k) * t(i,j,k)
c
c----------------------------------------------------------------
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
      pi32=pi**1.50000000000000000D0
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md
c  pointers
c
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
c
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, also F later
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
      M15=M13+M
c least squares
      M17=M15+MM
c vectors of MO
      M18=M17+MMd
c
      NCOa=NCO
      NCOb=NCO+Nunp
c end ------------------------------------------------
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2
c
      do 2 i=1,M
 2     Jx(i)=(M2-i)*(i-1)/2
c
      do 5 i=1,natom
      do 5 j=1,natom
       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
 5    continue
c
c
       if (OPEN) then
      do  k=1,Md
       aux(k)=af(k)+B(k,3)
      enddo
c 
       else
       do k=1,Md
        aux(k)=0.0D0
       enddo
       endif

      do 6 k=1,Md
       af(k)=af(k)+B(k,2)
 6     Rc(k)=0.D0
c     
c----------------------------------------------------------------
c
c------------------------------------------------------------------
c AF OF PREVIOUS ITERATION ARE USED IN THIS SUBROUTINE
c
c Calculation of all integrals , for constructing the
c Fock matrix
c Previously energy was computed, and coefficients for
c the fit were generated 
c------------------------------------------------------------------
c Second Part: initializations
c
      do 215 k=1,MM
 215   RMM(M5+k-1)=RMM(M11+k-1)
c
      if (OPEN) then
      do 219 k=1,MM
 219   RMM(M3+k-1)=RMM(M11+k-1)
      endif
c
c-------------------------------------------------------------------
c (ss|s)
c
      do 310 i=1,ns
      do 310 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
      kk=i+Jx(j)
c
      do 310 ni=1,ncont(i)
      do 310 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
c
       if (rexp.gt.rmax) goto 312
       sks=pi52*exp(-rexp)/zij
c
      do 311 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 311 nk=1,ncontd(k)
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      u=ad(k,nk)*zij/t0*dpc
      t1=ad(k,nk)*sqrt(t0)
      term=sks/t1*FUNCT(0,u)
      term=term*ccoef
c
      RMM(M5+kk-1)=RMM(M5+kk-1)+af(k)*term
      RMM(M3+kk-1)=RMM(M3+kk-1)+aux(k)*term
      Rc(k)=Rc(k)+RMM(kk)*term
 311   continue
 312   continue
 310   continue
c
c-------------------------------------------------------------
c (ps|s)
c
      do 320 i=ns+1,ns+np,3
      do 320 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 320 ni=1,ncont(i)
      do 320 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 322
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
c
       sks=pi52*exp(-rexp)/zij
c
      do 321 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 321 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      u=ad(k,nk)*ti*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
c
      do 325 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      term=ccoef*(t1*sss+t2*ss1s)
      ii=i+l1-1
c
c
      kk=ii+k1
      RMM(M5+kk-1)=RMM(M5+kk-1)+af(k)*term
      RMM(M3+kk-1)=RMM(M3+kk-1)+aux(k)*term
      Rc(k)=Rc(k)+RMM(kk)*term
 325   continue
 321   continue
 322   continue
 320   continue
c
c-------------------------------------------------------------
c (pp|s)
      do 330 i=ns+1,ns+np,3
      do 330 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 330 ni=1,ncont(i)
      do 330 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 332
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 331 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 331 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      ro=ad(k,nk)*ti
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-tj*ss1s)/z2
c
      do 335 l1=1,3
c
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 335 l2=1,lij
      t1=Q(l2)-r(Nuc(j),l2)
      t2=W(l2)-Q(l2)
      term=t1*ps+t2*p1s
c
      if (l1.eq.l2) then
       term=term+ta
      endif
c
      ii=i+l1-1
      jj=j+l2-1
c
      term=term*ccoef
c
      kk=ii+Jx(jj)
      RMM(M5+kk-1)=RMM(M5+kk-1)+af(k)*term
      RMM(M3+kk-1)=RMM(M3+kk-1)+aux(k)*term
      Rc(k)=Rc(k)+RMM(kk)*term
 335   continue
 331   continue
 332   continue
 330   continue
c
c-------------------------------------------------------------
c
c (ds|s)
      do 340 i=ns+np+1,M,6
      do 340 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 340 ni=1,ncont(i)
      do 340 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 342
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 341 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 341 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-roz*ss1s)/z2
c
      do 345 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
c
      do 345 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      term=t1*ps+t2*p1s
c
      f1=1.
      if (l1.eq.l2) then
       term=term+ta
       f1=sq3
      endif
c
      cc=ccoef/f1
      term=term*cc
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      kk=ii+k1
      RMM(M5+kk-1)=RMM(M5+kk-1)+af(k)*term
      RMM(M3+kk-1)=RMM(M3+kk-1)+aux(k)*term
      Rc(k)=Rc(k)+RMM(kk)*term
 345   continue
 341   continue
 342   continue
 340   continue
c
c-------------------------------------------------------------
c (dp|s)
      do 350 i=ns+np+1,M,6
      do 350 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 350 ni=1,ncont(i)
      do 350 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 352
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 351 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 351 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
c
      do 355 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      t5=(ps-roz*p1s)/z2
      p2s=t1*ss2s+t2*ss3s
c
      do 355 l2=1,l1
c
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      t6=(pjs-roz*pj1s)/z2
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
c
      f1=1.
      if (l1.eq.l2) then
       f1=sq3
       ds=ds+t3
       d1s=d1s+t4
      endif
c
      do 355 l3=1,3
c
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      term=t1*ds+t2*d1s
c
      if (l1.eq.l3) then
       term=term+t6
      endif
c
      if (l2.eq.l3) then
       term=term+t5
      endif
c
      l12=Ll(l1)+l2
      ii=i+l12-1
      jj=j+l3-1
c
      cc=ccoef/f1
      term=term*cc
c
      kk=ii+Jx(jj)
      RMM(M5+kk-1)=RMM(M5+kk-1)+af(k)*term
      RMM(M3+kk-1)=RMM(M3+kk-1)+aux(k)*term
      Rc(k)=Rc(k)+RMM(kk)*term
 355   continue
 351   continue
 352   continue
 350   continue
c
c-------------------------------------------------------------
c
c (dd|s)
      do 360 i=ns+np+1,M,6
      do 360 j=ns+np+1,i,6
c
      dd=d(Nuc(i),Nuc(j))
c
      do 360 ni=1,ncont(i)
      do 360 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 362
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 361 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 361 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
c
      do 365 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t6=(ps-roz*p1s)/z2
      t7=(p1s-roz*p2s)/z2
c
      do 365 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      t8=(pjs-roz*pj1s)/z2
      t9=(pj1s-roz*pj2s)/z2
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
c
      f1=1.
      if (l1.eq.l2) then
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t5
       f1=sq3
      endif
c
      t12=(ds-roz*d1s)/z2
c
c test now
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 365 l3=1,lij
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      pip=t1*ps+t2*p1s
      pi1p=t1*p1s+t2*p2s
      pjp=t1*pjs+t2*pj1s
      pj1p=t1*pj1s+t2*pj2s
      dp=t1*ds+t2*d1s
      d1p=t1*d1s+t2*d2s
c
      if (l1.eq.l3) then
       pip=pip+t3
       pi1p=pi1p+t4
       dp=dp+t8
       d1p=d1p+t9
      endif
c
      if (l2.eq.l3) then
       pjp=pjp+t3
       pj1p=pj1p+t4
       dp=dp+t6
       d1p=d1p+t7
      endif
c
      t10=(pjp-roz*pj1p)/z2
      t11=(pip-roz*pi1p)/z2
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      do 365 l4=1,lk
c
      t1=Q(l4)-r(Nuc(j),l4)
      t2=W(l4)-Q(l4)
      term=t1*dp+t2*d1p
c
      if (l1.eq.l4) then
       term=term+t10
      endif
c
      if (l2.eq.l4) then
       term=term+t11
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       term=term+t12
       f2=sq3
      endif
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
c
      ii=i+l12-1
      jj=j+l34-1
c
      cc=ccoef/(f1*f2)
      term=term*cc
c
      kk=ii+Jx(jj)
      RMM(M5+kk-1)=RMM(M5+kk-1)+af(k)*term
      RMM(M3+kk-1)=RMM(M3+kk-1)+aux(k)*term
      Rc(k)=Rc(k)+RMM(kk)*term
 365   continue
 361   continue
 362   continue
 360   continue
c
c-------------------------------------------------------------
c
c (ss|p)
      do 370 i=1,ns
      do 370 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
      kn=i+Jx(j)
c
      do 370 ni=1,ncont(i)
      do 370 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 372
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
c
       sks=pi52*exp(-rexp)/zij
c
      do 371 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 371 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      u=ad(k,nk)*ti*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
c
      do 375 l1=1,3
      t1=W(l1)-r(Nucd(k),l1)
      term=ccoef*t1*ss1s
      kk=k+l1-1
c
c
      Rc(kk)=Rc(kk)+RMM(kn)*term
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
 375   continue
 371   continue
 372   continue
 370   continue
c
c-------------------------------------------------------------
c
c (ps|p)
      do 380 i=ns+1,ns+np,3
      do 380 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 380 ni=1,ncont(i)
      do 380 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 382
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 381 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 381 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      u=ad(k,nk)*ti*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=ss1s/z2a
c
      do 385 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      p1s=t1*ss1s+t2*ss2s
c
      do 385 l2=1,3
      t2=W(l2)-r(Nucd(k),l2)
      term=t2*p1s
c
      if (l1.eq.l2) then
       term=term+ta
      endif
c
      ii=i+l1-1
      kk=k+l2-1
      term=term*ccoef
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 385   continue
 381   continue
 382   continue
 380   continue
c
c-------------------------------------------------------------
c
c (pp|p)
      do 390 i=ns+1,ns+np,3
      do 390 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 390 ni=1,ncont(i)
      do 390 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 392
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 391 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 391 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      t3=(ss1s-roz*ss2s)/z2
c
      do 395 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      t5=p1s/z2a
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 395 l2=1,lij
      t1=Q(l2)-r(Nuc(j),l2)
      t2=W(l2)-Q(l2)
      spj=t1*ss1s+t2*ss2s
      t4=spj/z2a
      pp=t1*p1s+t2*p2s
c
      if (l1.eq.l2) then
       pp=pp+t3
      endif
c
      do 395 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
       term=t1*pp
c
      if (l1.eq.l3) then
       term=term+t4
      endif
c
      if (l2.eq.l3) then
       term=term+t5
      endif
c
      ii=i+l1-1
      jj=j+l2-1
      kk=k+l3-1
c
      term=term*ccoef
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 395   continue
 391   continue
 392   continue
 390   continue
c
c-------------------------------------------------------------
c
c (ds|p)
      do 400 i=ns+np+1,M,6
      do 400 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 400 ni=1,ncont(i)
      do 400 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 402
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 401 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 401 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      t3=(ss1s-roz*ss2s)/z2
c
      do 405 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      p1s=t1*ss1s+t2*ss2s
      t5=p1s/z2a
      p2s=t1*ss2s+t2*ss3s
c
      do 405 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pj1s=t1*ss1s+t2*ss2s
      t6=pj1s/z2a
      ds=t1*p1s+t2*p2s
c
      f1=1.
      if (l1.eq.l2) then
       ds=ds+t3
       f1=sq3
      endif
c
      do 405 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      term=t1*ds
c
      if (l1.eq.l3) then
       term=term+t6
      endif
c
      if (l2.eq.l3) then
       term=term+t5
      endif
c
      l12=Ll(l1)+l2
      ii=i+l12-1
      kk=k+l3-1
c    
      cc=ccoef/f1
      term=term*cc
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 405   continue
 401   continue
 402   continue
 400   continue
c
c-------------------------------------------------------------
c
c (dp|p)
      do 410 i=ns+np+1,M,6
      do 410 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 410 ni=1,ncont(i)
      do 410 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 412
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 411 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 411 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(ss1s-roz*ss2s)/z2
      t4=(ss2s-roz*ss3s)/z2
c
      do 415 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      t5=(p1s-roz*p2s)/z2
      p3s=t1*ss3s+t2*ss4s
c
      do 415 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      t6=(pj1s-roz*pj2s)/z2
c
      f1=1.
      if (l1.eq.l2) then
       d1s=d1s+t3
       d2s=d2s+t4
       f1=sq3
      endif
c
      t9=d1s/z2a
c
      do 415 l3=1,3
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      d1p=t1*d1s+t2*d2s
      pi1p=t1*p1s+t2*p2s
      pj1p=t1*pj1s+t2*pj2s
c   
      if (l1.eq.l3) then
       d1p=d1p+t6
       pi1p=pi1p+t3
      endif
c
      if (l2.eq.l3) then
       d1p=d1p+t5
       pj1p=pj1p+t3
      endif
c
      t7=pi1p/z2a
      t8=pj1p/z2a
c
      do 415 l4=1,3
      t1=W(l4)-r(Nucd(k),l4)
      term=t1*d1p
c
      if (l1.eq.l4) then
       term=term+t8
      endif
c
      if (l2.eq.l4) then
       term=term+t7
      endif
c
      if (l3.eq.l4) then
       term=term+t9
      endif
c     
      l12=Ll(l1)+l2
      ii=i+l12-1
      jj=j+l3-1
      kk=k+l4-1
c
      cc=ccoef/f1
      term=term*cc
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 415   continue
 411   continue
 412   continue
 410   continue
c
c-------------------------------------------------------------
c
c (dd|p)
      do 420 i=ns+np+1,M,6
      do 420 j=ns+np+1,i,6
c
      dd=d(Nuc(i),Nuc(j))
c
      do 420 ni=1,ncont(i)
      do 420 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 422
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 421 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 421 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      ta=(sss-roz*ss1s)/z2
      t3=(ss1s-roz*ss2s)/z2
      t4=(ss2s-roz*ss3s)/z2
      t5=(ss3s-roz*ss4s)/z2
c
      do 425 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      t6=(p1s-roz*p2s)/z2
      t7=(p2s-roz*p3s)/z2
c
      do 425 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      t8=(pj1s-roz*pj2s)/z2
      t9=(pj2s-roz*pj3s)/z2
c
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
c
      f1=1.D0
      if (l1.eq.l2) then
       d1s=d1s+t3
       d2s=d2s+t4
       d3s=d3s+t5
       f1=sq3
      endif
c
      t18=(d1s-roz*d2s)/z2
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 425 l3=1,lij
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      d1pk=t1*d1s+t2*d2s
      d2pk=t1*d2s+t2*d3s
      pjpk=t1*pjs+t2*pj1s
      pj1pk=t1*pj1s+t2*pj2s
      pj2pk=t1*pj2s+t2*pj3s
      pipk=t1*ps+t2*p1s
      pi1pk=t1*p1s+t2*p2s
      pi2pk=t1*p2s+t2*p3s
      spk=t1*sss+t2*ss1s
      s1pk=t1*ss1s+t2*ss2s
      s2pk=t1*ss2s+t2*ss3s
      t10=(s1pk-roz*s2pk)/z2
c
c
      if (l1.eq.l3) then
       d1pk=d1pk+t8
       d2pk=d2pk+t9
       pipk=pipk+ta
       pi1pk=pi1pk+t3
       pi2pk=pi2pk+t4
      endif
c
      if (l2.eq.l3) then
       d1pk=d1pk+t6
       d2pk=d2pk+t7
       pjpk=pjpk+ta
       pj1pk=pj1pk+t3
       pj2pk=pj2pk+t4
      endif
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      t16=(pj1pk-roz*pj2pk)/z2
      t17=(pi1pk-roz*pi2pk)/z2
c
      do 425 l4=1,lk
      t1=Q(l4)-r(Nuc(j),l4)
      t2=W(l4)-Q(l4)
      d1d=t1*d1pk+t2*d2pk
c     
      pjdkl=t1*pj1pk+t2*pj2pk
      pidkl=t1*pi1pk+t2*pi2pk
      d1pl=t1*d1s+t2*d2s
c
c
      if (l1.eq.l4) then
       d1d=d1d+t16
       pidkl=pidkl+t10
       d1pl=d1pl+t8
      endif
c
      if (l2.eq.l4) then
       d1d=d1d+t17
       pjdkl=pjdkl+t10
       d1pl=d1pl+t6
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       d1d=d1d+t18
       pjdkl=pjdkl+t8
       pidkl=pidkl+t6
       f2=sq3
      endif
c
      t11=pjdkl/z2a
      t12=pidkl/z2a
      t13=d1pl/z2a
      t14=d1pk/z2a
      do 425 l5=1,3
c
      t1=W(l5)-r(Nucd(k),l5)
      term=t1*d1d
c
      if (l1.eq.l5) then
       term=term+t11
      endif
c
      if (l2.eq.l5) then
       term=term+t12
      endif
c
      if (l3.eq.l5) then
       term=term+t13
      endif
c
      if (l4.eq.l5) then
       term=term+t14
      endif
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
      ii=i+l12-1
      jj=j+l34-1
      kk=k+l5-1
c
      cc=ccoef/(f1*f2)
      term=term*cc
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 425   continue
 421   continue
 422   continue
 420   continue
c
c-------------------------------------------------------------
c
c (ss|d)
      do 430 i=1,ns
      do 430 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
      kn=i+Jx(j)
c
      do 430 ni=1,ncont(i)
      do 430 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 432
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 431 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 431 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=ti
c
      zc=2.D0*ad(k,nk)
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-roz*ss1s)/zc
c
      do 435 l1=1,3
      t1=W(l1)-r(Nucd(k),l1)
      ss1p=t1*ss2s
c
      do 435 l2=1,l1
      t1=W(l2)-r(Nucd(k),l2)
      term=t1*ss1p
c
      f1=1.D0
      if (l1.eq.l2) then
       term=term+ta
       f1=sq3
      endif
c
      cc=ccoef/f1
      term=term*cc
c
      l12=Ll(l1)+l2
      kk=k+l12-1
c
c
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 435   continue
 431   continue
 432   continue
 430   continue
c
c-------------------------------------------------------------
c
c (ps|d)
      do 440 i=ns+1,ns+np,3
      do 440 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 440 ni=1,ncont(i)
      do 440 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 442
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 441 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 441 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      zc=2.D0*ad(k,nk)
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      roz=ti
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      t3=ss2s/z2a
      ss3s=t2*FUNCT(3,u)
c
      do 445 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      t5=(ps-roz*p1s)/zc
c
      do 445 l2=1,3
      t1=W(l2)-r(Nucd(k),l2)
      sspj=t1*ss2s
      pispj=t1*p2s
c
      if (l1.eq.l2) then
       pispj=pispj+t3
      endif
c
      t4=sspj/z2a
      do 445 l3=1,l2
      t1=W(l3)-r(Nucd(k),l3)
      term=t1*pispj
c
      if (l1.eq.l3) then
       term=term+t4
      endif
c
      f1=1.D0
      if (l2.eq.l3) then
       term=term+t5
       f1=sq3
      endif
c
c
      l23=l2*(l2-1)/2+l3
      ii=i+l1-1
      kk=k+l23-1
c     
      cc=ccoef/f1
      term=term*cc
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 445  continue
 441  continue
 442  continue
 440  continue
c
c-------------------------------------------------------------
c
c (pp|d)
      do 450 i=ns+1,ns+np,3
      do 450 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 450 ni=1,ncont(i)
      do 450 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 452
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 451 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 451 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      zc=2.D0*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t6=ss2s/z2a
c
      do 455 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t8=p2s/z2a
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 455 l2=1,lij
      t1=Q(l2)-r(Nuc(j),l2)
      t2=W(l2)-Q(l2)
      pijs=t1*ps+t2*p1s
      pij1s=t1*p1s+t2*p2s
      pij2s=t1*p2s+t2*p3s
      spjs=t1*ss1s+t2*ss2s
      sp2js=t1*ss2s+t2*ss3s
      t7=sp2js/z2a
c
      if (l1.eq.l2) then
       pijs=pijs+t3
       pij1s=pij1s+t4
       pij2s=pij2s+t5
      endif
c
      t11=(pijs-ti*pij1s)/zc
c
      do 455 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      pp1p=t1*pij2s
      spjpk=t1*sp2js
      pispk=t1*p2s
c
      if (l1.eq.l3) then
       pp1p=pp1p+t7
       pispk=pispk+t6
      endif
c
      if (l2.eq.l3) then
       pp1p=pp1p+t8
       spjpk=spjpk+t6
      endif
c
      t9=spjpk/z2a
      t10=pispk/z2a
c
      do 455 l4=1,l3
      t1=W(l4)-r(Nucd(k),l4)
      term=t1*pp1p
c
      if (l1.eq.l4) then
       term=term+t9
      endif
c
      if (l2.eq.l4) then
       term=term+t10
      endif
c
      f1=1.D0
      if (l3.eq.l4) then
       term=term+t11
       f1=sq3
      endif
c
      ii=i+l1-1
      jj=j+l2-1
      l34=l3*(l3-1)/2+l4
      kk=k+l34-1
c
      cc=ccoef/f1
      term=term*cc
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 455   continue
 451   continue
 452   continue
 450   continue
c
c-------------------------------------------------------------
c
c (ds|d)
      do 460 i=ns+np+1,M,6
      do 460 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 460 ni=1,ncont(i)
      do 460 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 462
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 461 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 461 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t6=ss2s/z2a
c
      do 465 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t7=p2s/z2a
c
      do 465 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      t8=pj2s/z2a
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t5
       f1=sq3
      endif
c
      t11=(ds-ti*d1s)/zc
      do 465 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      ds1p=t1*d2s
      pjs1pk=t1*pj2s
      pis1pk=t1*p2s
c
      if (l1.eq.l3) then
       ds1p=ds1p+t8
       pis1pk=pis1pk+t6
      endif
c
      if (l2.eq.l3) then
       ds1p=ds1p+t7
       pjs1pk=pjs1pk+t6
      endif
c
      t9=pjs1pk/z2a
      t10=pis1pk/z2a
c
      do 465 l4=1,l3
       t1=W(l4)-r(Nucd(k),l4)
       term=t1*ds1p
c
      if (l1.eq.l4) then
       term=term+t9
      endif
c
      if (l2.eq.l4) then
       term=term+t10
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       term=term+t11
       f2=sq3
      endif
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
c
      ii=i+l12-1
      kk=k+l34-1
c
      cc=ccoef/(f1*f2)
      term=term*cc
c  
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 465   continue
 461   continue
 462   continue
 460   continue
c
c-------------------------------------------------------------
c
c (dp|d)
      do 470 i=ns+np+1,M,6
      do 470 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 470 ni=1,ncont(i)
      do 470 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 472
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 471 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 471 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
c
      do 475 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
c
      do 475 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+(sss-roz*ss1s)/z2
       d1s=d1s+(ss1s-roz*ss2s)/z2
       d2s=d2s+(ss2s-roz*ss3s)/z2
       d3s=d3s+(ss3s-roz*ss4s)/z2
       f1=sq3
      endif
c
      do 475 l3=1,3
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      spks=t1*ss2s+t2*ss3s
      dp=t1*ds+t2*d1s
      d1p=t1*d1s+t2*d2s
      d2p=t1*d2s+t2*d3s
      pi1p=t1*p1s+t2*p2s
      pi2p=t1*p2s+t2*p3s
      pj1p=t1*pj1s+t2*pj2s
      pj2p=t1*pj2s+t2*pj3s
c
      if (l1.eq.l3) then
       dp=dp+(pjs-roz*pj1s)/z2
       d1p=d1p+(pj1s-roz*pj2s)/z2
       d2p=d2p+(pj2s-roz*pj3s)/z2
       pi1p=pi1p+(ss1s-roz*ss2s)/z2
       pi2p=pi2p+(ss2s-roz*ss3s)/z2
      endif
c
      if (l2.eq.l3) then
       dp=dp+(ps-roz*p1s)/z2
       d1p=d1p+(p1s-roz*p2s)/z2
       d2p=d2p+(p2s-roz*p3s)/z2
       pj1p=pj1p+(ss1s-roz*ss2s)/z2
       pj2p=pj2p+(ss2s-roz*ss3s)/z2
      endif
c
      do 475 l4=1,3
      t1=W(l4)-r(Nucd(k),l4)
      dp1p=t1*d2p
      pjpkpl=t1*pj2p
      pipkpl=t1*pi2p
      dspl=t1*d2s
c
      if (l1.eq.l4) then
       dp1p=dp1p+pj2p/z2a
       pipkpl=pipkpl+spks/z2a
       dspl=dspl+pj2s/z2a
      endif
c
      if (l2.eq.l4) then
       dp1p=dp1p+pi2p/z2a
       pjpkpl=pjpkpl+spks/z2a
       dspl=dspl+p2s/z2a
      endif
c
      if (l3.eq.l4) then
       dp1p=dp1p+d2s/z2a
       pipkpl=pipkpl+p2s/z2a
       pjpkpl=pjpkpl+pj2s/z2a
      endif
c
      do 475 l5=1,l4
      t1=W(l5)-r(Nucd(k),l5)
      term=t1*dp1p
c
      if (l1.eq.l5) then
       term=term+pjpkpl/z2a
      endif
c
      if (l2.eq.l5) then
       term=term+pipkpl/z2a
      endif
c
      if (l3.eq.l5) then
       term=term+dspl/z2a
      endif
c
      f2=1.D0
      if (l4.eq.l5) then
       term=term+(dp-ro*d1p/ad(k,nk))/zc
       f2=sq3
      endif
c
      l12=Ll(l1)+l2
      ii=i+l12-1
      jj=j+l3-1
      l45=l4*(l4-1)/2+l5
      kk=k+l45-1
c
      cc=ccoef/(f1*f2)
      term=term*cc
c
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 475   continue
 471   continue
 472   continue
 470   continue
c
c-------------------------------------------------------------
c
c (dd|d)
      do 480 i=ns+np+1,M,6
      do 480 j=ns+np+1,i,6
c
      ddi=d(Nuc(i),Nuc(j))
c
      do 480 ni=1,ncont(i)
      do 480 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*ddi
       if (rexp.gt.rmax) goto 482
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 481 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 481 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      ss6s=t2*FUNCT(6,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t6=(ss3s-roz*ss4s)/z2
      t6b=(ss4s-roz*ss5s)/z2
c
      do 485 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      p5s=t1*ss5s+t2*ss6s
c
      t7=(ps-roz*p1s)/z2
      t8=(p1s-roz*p2s)/z2
      t9=(p2s-roz*p3s)/z2
      t10=(p3s-roz*p4s)/z2
c
      do 485 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      pj4s=t1*ss4s+t2*ss5s
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
      d4s=t1*p4s+t2*p5s
c
      t11=(pjs-roz*pj1s)/z2
      t12=(pj1s-roz*pj2s)/z2
      t13=(pj2s-roz*pj3s)/z2
      t14=(pj3s-roz*pj4s)/z2
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t5
       d3s=d3s+t6
       d4s=d4s+t6b
       f1=sq3
      endif
c
      t16=(ds-roz*d1s)/z2
      t17=(d1s-roz*d2s)/z2
      t18=(d2s-roz*d3s)/z2
      t22a=d2s/z2a
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 485 l3=1,lij
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      dpk=t1*ds+t2*d1s
      d1pk=t1*d1s+t2*d2s
      d2pk=t1*d2s+t2*d3s
      d3pk=t1*d3s+t2*d4s
      pjpk=t1*pjs+t2*pj1s
      pj1pk=t1*pj1s+t2*pj2s
      pj2pk=t1*pj2s+t2*pj3s
      pj3pk=t1*pj3s+t2*pj4s
      pipk=t1*ps+t2*p1s
      pi1pk=t1*p1s+t2*p2s
      pi2pk=t1*p2s+t2*p3s
      pi3pk=t1*p3s+t2*p4s
      spk=t1*sss+t2*ss1s
      s1pk=t1*ss1s+t2*ss2s
      s2pk=t1*ss2s+t2*ss3s
      s3pk=t1*ss3s+t2*ss4s
c
      t15=(s2pk-roz*s3pk)/z2
c
      if (l1.eq.l3) then
       dpk=dpk+t11
       d1pk=d1pk+t12
       d2pk=d2pk+t13
       d3pk=d3pk+t14
       pipk=pipk+t3
       pi1pk=pi1pk+t4
       pi2pk=pi2pk+t5
       pi3pk=pi3pk+t6
      endif
c
      if (l2.eq.l3) then
       dpk=dpk+t7
       d1pk=d1pk+t8
       d2pk=d2pk+t9
       d3pk=d3pk+t10
       pjpk=pjpk+t3
       pj1pk=pj1pk+t4
       pj2pk=pj2pk+t5
       pj3pk=pj3pk+t6
      endif
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      t20=pj2pk/z2a
      t21=pi2pk/z2a
      t22=d2pk/z2a
c
      t24=(pjpk-roz*pj1pk)/z2
      t25=(pj1pk-roz*pj2pk)/z2
      t26=(pj2pk-roz*pj3pk)/z2
      t27=(pipk-roz*pi1pk)/z2
      t28=(pi1pk-roz*pi2pk)/z2
      t29=(pi2pk-roz*pi3pk)/z2
c
      do 485 l4=1,lk
      t1=Q(l4)-r(Nuc(j),l4)
      t2=W(l4)-Q(l4)
      dd=t1*dpk+t2*d1pk
      d1d=t1*d1pk+t2*d2pk
      d2d=t1*d2pk+t2*d3pk
c
      pjdkl=t1*pj2pk+t2*pj3pk
      pidkl=t1*pi2pk+t2*pi3pk
      d2pl=t1*d2s+t2*d3s
c
      sdkl=t1*s2pk+t2*s3pk
      pj2pl=t1*pj2s+t2*pj3s
      pi2pl=t1*p2s+t2*p3s
c
      if (l1.eq.l4) then
       dd=dd+t24
       d1d=d1d+t25
       d2d=d2d+t26
       pidkl=pidkl+t15
       d2pl=d2pl+t13
       pi2pl=pi2pl+t5
      endif
c
      if (l2.eq.l4) then
       dd=dd+t27
       d1d=d1d+t28
       d2d=d2d+t29
       pjdkl=pjdkl+t15
       d2pl=d2pl+t9
       pj2pl=pj2pl+t5
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       sdkl=sdkl+t5
       dd=dd+t16
       d1d=d1d+t17
       d2d=d2d+t18
       pjdkl=pjdkl+t13
       pidkl=pidkl+t9
       f2=sq3
      endif
c
      t30=pjdkl/z2a
      t40=pidkl/z2a
      t50=sdkl/z2a
      t60=pj2pl/z2a
      t70=pi2pl/z2a
      t80=d2pl/z2a
      t23=(dd-ti*d1d)/zc
c
      do 485 l5=1,3
c
      t1=W(l5)-r(Nucd(k),l5)
      ddp=t1*d2d
      pjdklp=t1*pjdkl
      pidklp=t1*pidkl
      dijplp=t1*d2pl
      dijpkp=t1*d2pk
c
      if (l1.eq.l5) then
       ddp=ddp+t30
       pidklp=pidklp+t50
       dijplp=dijplp+t60
       dijpkp=dijpkp+t20
      endif
c
      if (l2.eq.l5) then
       ddp=ddp+t40
       pjdklp=pjdklp+t50
       dijplp=dijplp+t70
       dijpkp=dijpkp+t21
      endif
c
      if (l3.eq.l5) then
       ddp=ddp+t80
       pjdklp=pjdklp+t60
       pidklp=pidklp+t70
       dijpkp=dijpkp+t22a
      endif
c
      if (l4.eq.l5) then
       ddp=ddp+t22
       pjdklp=pjdklp+t20
       pidklp=pidklp+t21
       dijplp=dijplp+t22a
      endif
c
      t31=pjdklp/z2a
      t41=pidklp/z2a
      t51=dijplp/z2a
      t61=dijpkp/z2a
c
      do 485 l6=1,l5
c
      t1=W(l6)-r(Nucd(k),l6)
      term=t1*ddp
c
      if (l1.eq.l6) then
       term=term+t31
      endif
c
      if (l2.eq.l6) then
       term=term+t41
      endif
c
      if (l3.eq.l6) then
       term=term+t51
      endif
c
      if (l4.eq.l6) then
       term=term+t61
      endif
c
      f3=1.D0
      if (l5.eq.l6) then
       term=term+t23
       f3=sq3
      endif
c
      cc=ccoef/(f1*f2*f3)
      term=term*cc
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
      l56=Ll(l5)+l6
      ii=i+l12-1
      jj=j+l34-1
      kk=k+l56-1
c
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af(kk)*term
      RMM(M3+kn-1)=RMM(M3+kn-1)+aux(kk)*term
      Rc(kk)=Rc(kk)+RMM(kn)*term
 485  continue
 481  continue
 482  continue
 480  continue
c
c----------------------------------------------------------------
c
      Ndens=2
c
c call fit for exchange correlation routine
      if (integ) then
       do i=1,Md
        B(i,2)=0.0D0
        B(i,1)=0.0D0
        B(i,3)=0.0D0
       enddo
c
       else
c
      if (dens) then
      call exch(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,
     > Md,ncontd,nshelld,Nucd,cd,ad,RMM,NCOa,NCOb,M17,B)
      else
c
      call exch2(OPEN,Iz,natom,RMM,nshelld,M,Md,M17,NCOa,NCOb,B)
c
      endif
      endif
c
      Ex=0.D0
c end of Least squares option ----------
c
c
c-------------------------------------------------------------
c CALCULATION OF FITTING COEFFICIENTS ( TO BE USED IN NEXT ITERATION)
c
c Constraint that integrated fitted density = N electrons
c
cSVD PART  ------------------------------------------------
c
      if (SVD) then
c
      MMp=Md*(Md+1)/2
      do 199 k=1,MMp
 199   RMM(M9+k-1)=0.0D0
      do 208 k=1,Md
       af(k)=0.0D0
 208   RMM(M9+k-1)=Rc(k)
c
       k1=0
       do 116 j=1,Md
       do 116 i=j,Md
c
       k1=k1+1
c
        X(i,j)=RMM(M7+k1-1)
        X(j,i)=X(i,j)
 116   continue
c
      M10=M9+Md
      M12=M10+Md
      Md3=3*Md
c ESSL OPTION ------------------------
#ifdef essl
      CALL DGESVF(2,X,Md,RMM(M9),Md,1,RMM(M10),
     >             Md,Md,RMM(M12),Md3)
      imax=idamax(Md,RMM(M10),1)
      ss=RMM(M10+imax-1)
      tau=0.22204D-15*Md*ss*10
c     write(*,*) ss,tau
c
      CALL DGESVS(X,Md,RMM(M9),Md,1,RMM(M10),af,Md,Md,Md,tau)
#endif
c---------------------------------------
c LAPACK OPTION
#ifdef pack
      do i=1,Md
        af(i)=Rc(i)
      enddo
c
      Md5=5*Md
      rcond=1.0D-06
      call dgelss(Md,Md,1,X,Md,af,Md,RMM(M9),rcond,irank,RMM(M10),
     >            Md5,info)
c
c     
#endif
c------------------------------------------
c
c END SVD PART --
c
c if SVD.eq.false, then goes to Normal equation method, with or without
c constraint
      else
c---------------------------------------
c NORMAL EQUATION PART
      if (iconst.eq.1) then
c P : integrals of fitting functions --------------
c
      do 294 k=1,Md
 294   P(k)=0.0D0
c
      do 295 k=1,nsd
       do 295 nk=1,ncontd(k)
c
 295   P(k)=P(k)+cd(k,nk)/sqrt(ad(k,nk)**3)
c
c p functions
c all 0
c
c d functions ----
      do 297 k=nsd+npd+1,Md,6
       do 297 nk=1,ncontd(k)
c
        t0=cd(k,nk)/(sqrt(ad(k,nk)**3)*2.D0*ad(k,nk))
       do 297 l1=1,3
       do 297 l2=1,l1
        kk=k+Ll(l1)+l2-1
c
        if (l1.eq.l2) then
         f1=sq3
         P(kk)=P(kk)+t0/f1
        endif
c
 297    continue
c
        do 298 k=1,Md
 298     P(k)=P(k)*pi32
c--------------------------------------
      do 300 m1=1,Md
       FF(m1)=0.0D0
       do 301 k=1,m1-1
 301   FF(m1)=FF(m1)+P(k)*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
       do 302 k=m1,Md
 302   FF(m1)=FF(m1)+P(k)*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
 300  continue
c
      r0=0.0D0
      r1=0.0D0
      do 900 m1=1,Md
       r0=r0+FF(m1)*Rc(m1)
       r1=r1+FF(m1)*P(m1)
 900  continue
c
      Nel=2*NCO+Nunp
      bda=(Nel-r0)/r1
c
c
      do 200 m1=1,Md
      af(m1)=0.D0
      do 201 k=1,m1-1
 201  af(m1)=af(m1)+(Rc(k)+bda*P(k))*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
      do 202 k=m1,Md
 202  af(m1)=af(m1)+(Rc(k)+bda*P(k))*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
       ss9=ss9+af(m1)*P(m1)
 200  continue

c no constraint applied
      else
      do 1200 m1=1,Md
      af(m1)=0.D0
      do 1201 k=1,m1-1
 1201  af(m1)=af(m1)+Rc(k)*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
      do 1202 k=m1,Md
 1202  af(m1)=af(m1)+Rc(k)*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
 1200  continue

      endif
      endif
c----------------------------------------------------------------
c
      Ea=0.D0
      Eb=0.D0
c
      do 210 m1=1,Md
       Ex=Ex+B(m1,1)*Rc(m1)
       Ea=Ea+af(m1)*Rc(m1)
      do 211 k=1,m1
 211   Eb=Eb+af(k)*af(m1)*RMM(M7+m1+(2*Md-k)*(k-1)/2-1)
      do 212 k=m1+1,Md
 212   Eb=Eb+af(k)*af(m1)*RMM(M7+k+(2*Md-m1)*(m1-1)/2-1)
 210  continue
c
c
c Preparation of the coefficients for next iteration
c
       if (integ) then
        NCOa=NCO
        NCOb=NCO+Nunp
        write(*,*) 'int3N'
#ifdef GPU
        call exchnum_gpu(NORM, natom, r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >    M18,M5,NCOa,Ex,nopt,Iexch, igrid2, e_, e_2, e3, wang, wang2,
     >    wang3,Ndens, 1)        
#else      
       call EXCHFOCK(OPEN,NORM,natom,Iz,Nuc,ncont,nshell,a,c,r,
     >               M,M18,NCOa,NCOb,RMM,Ex)
       write(*,*) 'energia fortran (double)',Ex
#endif
       endif
c
c
c
      E2=Ea-Eb/2.D0+Ex
c
c------------------------------------------------------------------
c
      return
      end
