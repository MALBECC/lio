c Integrals subroutine -Third part c 2 e integrals, 3 index : wavefunction and density fitting functions
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
       module subm_int3mem; contains
       subroutine int3mem()

       use liotemp   , only: FUNCT
       use garcha_mod, only: RMM, ll, cool, cools, kkind, kkinds
     >                     , nuc, nucd, a, c, d, r, ad, cd, natomc
     >                     , nns, nnp, nnd, nnps, nnpp, nnpd, jatc
     >                     , ncont, ncontd, nshell, nshelld, M, Md
     >                     , rmax, rmaxs, pi52, Nunp, NORM, NCO
     >                     , kknums, kknumd
c
c
       implicit none
c
       real*8  :: Q(3), W(3)
       integer, dimension(:), allocatable :: Jx
c scratch space

! Eliminating implicits:
       real*8  :: alf, cc, ccoef, f1, f2, f3
       real*8  :: rexp, ro, roz, sq3, term, u, ddij
       real*8  :: tii, tjj, NCOb, NCOa
       real*8  :: z2, z2a, zc, zij

       real*8  :: d1s, d1p, d1d, d1pk, d1pl
       real*8  :: d2s, d2p, d2d, d2pl, d2pk
       real*8  :: d3s, d3pk, d4s
       real*8  :: ds, ds1p, dspl, dp, dpc, dpk, dp1p
       real*8  :: dd, ddp, dijplp, dijpkp

       real*8  :: ps, pp, pp1p, p1s, p2s, p3s, p4s, p5s
       real*8  :: pi1p, pi1pk, pi2p, pi2pk, pi2pl, pi3pk
       real*8  :: pijs, pij1s, pij2s, pispj, pispk, pis1pk
       real*8  :: pip, pipk, pipkpl, pidkl, pidklp
       real*8  :: pjs, pjs1pk
       real*8  :: pj1s, pj1p, pj1pk, pj2s, pj2p, pj2pk, pj2pl
       real*8  :: pj3s, pj3pk, pj4s
       real*8  :: pjp, pjpk, pjpkpl, pjdklp,pjdkl

       real*8  :: sss, sks, spk, spjs, sspj, spjpk, sp2js
       real*8  :: ss1s, ss2s, ss3s, ss4s, ss5s, ss6s
       real*8  :: ss1p, s1pk, s2pk, s3pk, spks, spj, sdkl

       real*8  :: ta, ti, tj
       real*8  :: t0, t1, t2, t3, t4, t5, t6, t6b, t7, t8, t9
       real*8  :: t10, t11, t12, t13, t14, t15, t16, t17, t18
       real*8  :: t20, t21, t22, t22a, t23, t24, t25, t26
       real*8  :: t27, t28, t29, t30, t31
       real*8  :: t40, t41, t50, t51, t60, t61, t70, t80

       integer :: MM, MMd, Md2
       integer :: M1, M2, M3, M5, M7, M9, M11, M13, M15, M17, M18
       integer :: i, ii, j, jj, k, kk, ni, nj, nk, kn, k1
       integer :: id, iki, jki, ix, kknan, knan, kknumsmax
       integer :: ns, nsd, nd, ndd, np, npd
       integer :: l, lk, l1, l2, l3, l4, l5, l6
       integer :: lij, l12, l23, l34, l45, l56


c auxiliars
       logical :: fato, fato2
       allocate (Jx(M))
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
c
c------------------------------------------------------------------
c # of calls
c      Ncall=Ncall+1
c      Ndens=1
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
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
*  Mmem is the pointer to the RAM memory fo two-e integrals
*  Mmem = M20 in CLOSED SHELL case,  Mmem = M23 in OPEN SHELL case
        NCOa=NCO
        NCOb=NCO+Nunp
c
c end ------------------------------------------------
      if (NORM) then
      sq3=dsqrt(3.D0)
      else
      sq3=1.D0
      sq3=1.D0
      endif
c
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2


c
      do 2 i=1,M
 2     Jx(i)=(M2-i)*(i-1)/2
c
c      do 5 i=1,natom
c      do 5 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 5    continue
c
c      do 6 k=1,Md
c 6     Rc(k)=0.D0

      ix=0
      kknumd=0
      kknums=0
c para dimecionar cool
      do 3114 i=1,ns


          do 3114 knan=1,natomc(nuc(i))
           j=nnps(jatc(knan,nuc(i)))-1
          do 3114 kknan=1,nns(jatc(knan,nuc(i)))
            j=j+1
            if(j.le.i) then

          kk=i+Jx(j)
          dd=d(Nuc(i),Nuc(j))
c
      fato=.true.
      fato2=.true.

      do 310 ni=1,ncont(i)
      do 310 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then
           if (rexp.lt.rmaxs) then
               if (fato) then
                 kknumd = kknumd+1
                 fato=.false.
               endif
           else
               if (fato2) then
                kknums = kknums+1
                fato2=.false.
               endif
           endif
       endif
310    continue
       endif
3114   continue
      do 320 i=ns+1,ns+np,3
       do 320 knan=1,natomc(nuc(i))

         j=nnps(jatc(knan,nuc(i)))-1

      do 320 kknan=1,nns(jatc(knan,nuc(i)))
        j=j+1


c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)

      fato=.true.
      fato2=.true.
c1
      do 320 ni=1,ncont(i)
      do 320 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then
       if (rexp.lt.rmaxs) then
       if (fato) then
        do iki=1,3
        kknumd=kknumd+1
         enddo
        fato=.false.
        endif
        else
       if (fato2) then
        do iki=1,3
        kknums=kknums+1
         enddo
        fato2=.false.
        endif

         endif
        endif
320    continue

      do 3333 i=ns+1,ns+np,3
         do 3333 knan=1,natomc(nuc(i))
         j=nnpp(jatc(knan,nuc(i)))-3
      do 3333 kknan=1,nnp(jatc(knan,nuc(i))),3
        j=j+3
       if(j.le.i) then
        fato=.true.
        fato2=.true.
      dd=d(Nuc(i),Nuc(j))
      do 330 ni=1,ncont(i)
      do 330 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then
         if (rexp.lt.rmaxs) then
           if (fato) then
             if(i.eq.j) then
                do iki=1,3
                do jki=1,iki
                  kknumd=kknumd+1
                enddo
                enddo
             else
                do iki=1,3
                do  jki=1,3
                 kknumd=kknumd+1
                enddo
                enddo
             endif
            fato=.false.
           endif
         else
         if (fato2) then
             if(i.eq.j) then
                do iki=1,3
                do jki=1,iki
                  kknums=kknums+1
                enddo
                enddo
             else
                do iki=1,3
                do  jki=1,3
                 kknums=kknums+1
                enddo
                enddo
             endif
            fato2=.false.
           endif


       endif
       endif     
330    continue
        endif
3333   continue   
      do 340 i=ns+np+1,M,6
         do 340 knan=1,natomc(nuc(i))
         j=nnps(jatc(knan,nuc(i)))-1
      do 340 kknan=1,nns(jatc(knan,nuc(i)))
        j=j+1
       fato=.true.
       fato2=.true.
      k1=Jx(j)
      dd=d(Nuc(i),Nuc(j))
      do 340 ni=1,ncont(i)
      do 340 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then
         if (rexp.lt.rmaxs) then
            if (fato) then
               do iki=1,6
               kknumd=kknumd+1
               enddo
              fato=.false.
            endif
         else
            if (fato2) then
              do iki=1,6
                kknums=kknums+1
              enddo
              fato2=.false.
            endif
         endif
        endif
340    continue

      do 350 i=ns+np+1,M,6
         do 350 knan=1,natomc(nuc(i))
         j=nnpp(jatc(knan,nuc(i)))-3
      do 350 kknan=1,nnp(jatc(knan,nuc(i))),3
        j=j+3
       fato=.true.
       fato2=.true.
      dd=d(Nuc(i),Nuc(j))
      do 350 ni=1,ncont(i)
      do 350 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then
       if (rexp.lt.rmaxs) then
       if (fato) then
        do iki=1,6
          do  jki=1,3
             kknumd=kknumd+1
            enddo
         enddo
        fato=.false.
       endif
       else
       if (fato2) then
        do iki=1,6
          do  jki=1,3
             kknums=kknums+1
            enddo
         enddo
        fato2=.false.
       endif

       endif
       endif
350    continue
      do 3666 i=ns+np+1,M,6
         do 3666 knan=1,natomc(nuc(i))
         j=nnpd(jatc(knan,nuc(i)))-6
      do 3666 kknan=1,nnd(jatc(knan,nuc(i))),6
        j=j+6
       if(j.le.i) then
      fato=.true.
      fato2=.true.
      dd=d(Nuc(i),Nuc(j))
      do 360 ni=1,ncont(i)
      do 360 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then
       if (rexp.lt.rmaxs) then
       if (fato) then
         fato=.false.
          if(i.eq.j) then
             do iki=1,6
               do jki=1,iki
                 kknumd=kknumd+1
               enddo
             enddo
          else
            do iki=1,6
              do jki=1,6
                kknumd=kknumd+1
              enddo
             enddo
         endif
       endif
       else
       if (fato2) then
         fato2=.false.
          if(i.eq.j) then
             do iki=1,6
               do jki=1,iki
                 kknums=kknums+1
               enddo
             enddo
          else
            do iki=1,6
              do jki=1,6
                kknums=kknums+1
              enddo
             enddo
         endif
       endif

       endif
       endif
360     continue
         endif
3666    continue

      
        allocate(cool(kknumd*Md),kkind(kknumd))
        allocate(cools(kknums*Md),kkinds(kknums))
        write(*,*) 'kknum y kknums',kknumd,kknums

         kknumsmax=kknums
        cool=0
        cools=0
        kknumd=0
        kknums=0

c------------------------------------------------------------------
c (ss|s)
c
      do 114 i=1,ns


          do 114 knan=1,natomc(nuc(i))
           j=nnps(jatc(knan,nuc(i)))-1
          do 114 kknan=1,nns(jatc(knan,nuc(i)))
            j=j+1
            if(j.le.i) then
           
          kk=i+Jx(j)
          dd=d(Nuc(i),Nuc(j))
c
      fato=.true.
      fato2=.true.

      do 10 ni=1,ncont(i)
      do 10 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then
          if (rexp.lt.rmaxs) then
            if (fato) then
              kknumd = kknumd+1
              kkind(kknumd)=kk
              fato=.false. 
            endif
          else
            if (fato2) then
               kknums = kknums+1 
              if(kknumsmax.lt.kknums) stop '1'
               kkinds(kknums)=kk
              fato2=.false.
         endif


       endif
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 11 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 11 nk=1,ncontd(k)
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      u=ad(k,nk)*zij/t0*dpc
      t1=ad(k,nk)*dsqrt(t0)
      term=sks/t1*FUNCT(0,u)
      term=term*ccoef
c
           if (rexp.lt.rmaxs) then
             id = (kknumd-1)*md+k
             cool(id) = cool(id) + term
           else
             id = (kknums-1)*md+k
            cools(id) = cools(id) + real(term)
          endif

      ix=ix+1
 11   continue

c----(ss|p)
      do 71 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 71 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      tii=zij/t0
      tjj=ad(k,nk)/t0
      W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
      W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
      W(3)=tii*Q(3)+tjj*r(Nucd(k),3)
c
c
      u=ad(k,nk)*tii*dpc
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
c
      do 75 l1=1,3
      t1=W(l1)-r(Nucd(k),l1)
      term=ccoef*t1*ss1s
      kn=k+l1-1
c
c
      if (rexp.lt.rmaxs) then
       id = (kknumd-1)*Md+kn
      cool(id) = cool(id) + term
      else
       id = (kknums-1)*Md+kn
      cools(id) = cools(id) + real(term)
      endif 
  
       ix=ix+1
 75   continue
 71   continue


c---- (ss|d)

      do 131 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 131 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      tii=zij/t0
      tjj=ad(k,nk)/t0
      W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
      W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
      W(3)=tii*Q(3)+tjj*r(Nucd(k),3)
c
      zc=2.D0*ad(k,nk)
      roz=tii
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-roz*ss1s)/zc
c
      do 135 l1=1,3
      t1=W(l1)-r(Nucd(k),l1)
      ss1p=t1*ss2s
c
      do 135 l2=1,l1
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
      kn=k+l12-1
c
c
c
          if (rexp.lt.rmaxs) then
             id = (kknumd-1)*Md+kn
            cool(id) = cool(id) + term
         else
            id = (kknums-1)*Md+kn
           cools(id) = cools(id) + real(term)
         endif

      ix=ix+1

 135   continue
 131   continue



       endif
 10   continue
      endif
114   continue
c
c-------------------------------------------------------------
c (ps|s)
c
      do 20 i=ns+1,ns+np,3
       do 20 knan=1,natomc(nuc(i))

         j=nnps(jatc(knan,nuc(i)))-1

      do 20 kknan=1,nns(jatc(knan,nuc(i)))
        j=j+1


c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)

      fato=.true.
      fato2=.true.
c1
      do 20 ni=1,ncont(i)
      do 20 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then   
       if (rexp.lt.rmaxs) then   
       if (fato) then
        do iki=1,3
        kknumd=kknumd+1
         kkind(kknumd)=i+k1+iki-1
         enddo
        fato=.false.
       endif
       else
       if (fato2) then
        do iki=1,3
        kknums=kknums+1
              if(kknumsmax.lt.kknums) stop '2'
         kkinds(kknums)=i+k1+iki-1
         enddo
        fato2=.false.
       endif

       endif
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij

c
      do 21 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 21 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      tii=zij/t0
      tjj=ad(k,nk)/t0
      W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
      W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
      W(3)=tii*Q(3)+tjj*r(Nucd(k),3)
c
c
      u=ad(k,nk)*tii*dpc
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1

      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
c
      do 21 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      term=ccoef*(t1*sss+t2*ss1s)
c      ii=i+l1-1
c
c
c       kk=ii+k1
      if(rexp.lt.rmaxs) then
       kk=l1-1+kknumd-2
      id = (kk-1)*md+k
      cool(id) = cool(id) + term
      else
       kk=l1-1+kknums-2
      id = (kk-1)*md+k
      cools(id) = cools(id) + real(term)

      endif


      ix=ix+1
 21   continue
c------(ps|p)

      do 81 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2



      do 81 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.*t0
c
      tii=zij/t0
      tjj=ad(k,nk)/t0
      W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
      W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
      W(3)=tii*Q(3)+tjj*r(Nucd(k),3)
c
c
      u=ad(k,nk)*tii*dpc
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=ss1s/z2a
c
      do 85 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      p1s=t1*ss1s+t2*ss2s
c
      do 85 l2=1,3
      t2=W(l2)-r(Nucd(k),l2)
      term=t2*p1s
c
      if (l1.eq.l2) then
       term=term+ta
      endif
c
c      ii=i+l1-1
      kk=k+l2-1
      term=term*ccoef
c
      if(rexp.lt.rmaxs) then
      kn=l1-1+kknumd-2
      id = (kn-1)*Md+kk
      cool(id) = cool(id) + term
       else
      kn=l1-1+kknums-2
      id = (kn-1)*Md+kk
      cools(id) = cools(id) + real(term)

      endif
      ix=ix+1
c      write(*,*) 'ninteg',ninteg
 85   continue
 81   continue

c-------(ps|d)

      do 141 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c

      do 141 nk=1,ncontd(k)
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
      roz=ti
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      t3=ss2s/z2a
      ss3s=t2*FUNCT(3,u)
c
      do 145 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      t5=(ps-roz*p1s)/zc
c
      do 145 l2=1,3
      t1=W(l2)-r(Nucd(k),l2)
      sspj=t1*ss2s
      pispj=t1*p2s
c
      t4=sspj/z2a
      if (l1.eq.l2) then
       pispj=pispj+t3
      endif
c
      do 145 l3=1,l2
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
      l23=l2*(l2-1)/2+l3
c      ii=i+l1-1
      kk=k+l23-1
c     
      cc=ccoef/f1
      term=term*cc
c
c      kn=ii+k1
      if(rexp.lt.rmaxs) then
      kn=l1-1+kknumd-2
      id = (kn-1)*Md+kk
      cool(id) = cool(id) + term
      else
      kn=l1-1+kknums-2
      id = (kn-1)*Md+kk
      cools(id) = cools(id) + real(term)
      endif

      ix=ix+1
 145  continue
 141  continue



      endif
 20   continue


c
c-------------------------------------------------------------
c (pp|s)
      do 333 i=ns+1,ns+np,3
         do 333 knan=1,natomc(nuc(i))
         j=nnpp(jatc(knan,nuc(i)))-3
      do 333 kknan=1,nnp(jatc(knan,nuc(i))),3
        j=j+3
       if(j.le.i) then
        fato=.true.
        fato2=.true.
      dd=d(Nuc(i),Nuc(j))
      do 30 ni=1,ncont(i)
      do 30 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then   
       if (rexp.lt.rmaxs) then   
       if (fato) then
          if(i.eq.j) then
            do iki=1,3
             do jki=1,iki
              kknumd=kknumd+1
               kkind(kknumd)=i+iki-1+Jx(j+jki-1)
            enddo
           enddo       
          else 
            do iki=1,3
            do  jki=1,3
             kknumd=kknumd+1
             kkind(kknumd)=i+iki-1+Jx(j+jki-1)
            enddo
            enddo
          endif
        fato=.false.
       endif

       else
        if (fato2) then
          if(i.eq.j) then
            do iki=1,3
             do jki=1,iki
              kknums=kknums+1
              if(kknumsmax.lt.kknums) stop '3'
               kkinds(kknums)=i+iki-1+Jx(j+jki-1)
            enddo
           enddo
          else
            do iki=1,3
            do  jki=1,3
             kknums=kknums+1
              if(kknumsmax.lt.kknums) stop '4'
             kkinds(kknums)=i+iki-1+Jx(j+jki-1)
            enddo
            enddo
          endif
        fato2=.false.
       endif

       endif
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 31 k=1,nsd
c
         
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 31 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-tj*ss1s)/z2
c
      ii=0
      do 35 l1=1,3
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
      do 35 l2=1,lij
      t1=Q(l2)-r(Nuc(j),l2)
      t2=W(l2)-Q(l2)
      term=t1*ps+t2*p1s
c
      if (l1.eq.l2) then
       term=term+ta
      endif
c
c      ii=i+l1-1
c      jj=j+l2-1
c
      term=term*ccoef
      ii=ii+1
c
c      kk=ii+Jx(jj)
      if (rexp.lt.rmaxs) then   
      if(i.eq.j) then
      kk=ii+kknumd-6
      else
      kk=ii+kknumd-9      
      endif
      id = (kk-1)*Md+k
      cool(id) = cool(id) + term
      else
      if(i.eq.j) then
      kk=ii+kknums-6
      else
      kk=ii+kknums-9
      endif
      id = (kk-1)*Md+k
      cools(id) = cools(id) + real(term)
      endif
      ix=ix+1
 35   continue
 31   continue

c (pp|p)

      do 91 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 91 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      t3=(ss1s-roz*ss2s)/z2
c
      ii=0
      do 95 l1=1,3
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
      do 95 l2=1,lij
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
       ii=ii+1
      do 95 l3=1,3
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
c     ii=i+l1-1
c     jj=j+l2-1
      kk=k+l3-1
c
      term=term*ccoef
c
c     kn=ii+Jx(jj)
      if (rexp.lt.rmaxs) then   
      if(i.eq.j) then
      kn=ii+kknumd-6
      else
      kn=ii+kknumd-9      
      endif
      id = (kn-1)*Md+kk
      cool(id) = cool(id) + term

      else
      if(i.eq.j) then
      kn=ii+kknums-6
      else
      kn=ii+kknums-9
      endif
      id = (kn-1)*Md+kk
      cools(id) = cools(id) + real(term)


      endif
      ix=ix+1
c      write(*,*) 'ninteg',ninteg
 95   continue
 91   continue


c(pp|d)

      do 151 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 151 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
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
      ii=0
c
      do 155 l1=1,3
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
      do 155 l2=1,lij
      t1=Q(l2)-r(Nuc(j),l2)
      t2=W(l2)-Q(l2)
      pijs=t1*ps+t2*p1s
      pij1s=t1*p1s+t2*p2s
      pij2s=t1*p2s+t2*p3s
      spjs=t1*ss1s+t2*ss2s
      sp2js=t1*ss2s+t2*ss3s
      t7=sp2js/z2a
c
       ii=ii+1
c
      if (l1.eq.l2) then
       pijs=pijs+t3
       pij1s=pij1s+t4
       pij2s=pij2s+t5
      endif
c
      t11=(pijs-ti*pij1s)/zc
c
      do 155 l3=1,3
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
      do 155 l4=1,l3
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
c      ii=i+l1-1
c      jj=j+l2-1
      l34=l3*(l3-1)/2+l4
      kk=k+l34-1
c
      cc=ccoef/f1
      term=term*cc
c
c      kn=ii+Jx(jji)
      if (rexp.lt.rmaxs) then

      if(i.eq.j) then
      kn=ii+kknumd-6
      else
      kn=ii+kknumd-9      
      endif
      id = (kn-1)*Md+kk

      cool(id) = cool(id) + term
      else

      if(i.eq.j) then
      kn=ii+kknums-6
      else
      kn=ii+kknums-9
      endif
      id = (kn-1)*Md+kk

      cools(id) = cools(id) + real(term)
      endif
      ix=ix+1

 155   continue
 151   continue


      endif
 30   continue
      endif
333   continue
c
c-------------------------------------------------------------
c
c (ds|s)
      do 40 i=ns+np+1,M,6
         do 40 knan=1,natomc(nuc(i))
         j=nnps(jatc(knan,nuc(i)))-1
      do 40 kknan=1,nns(jatc(knan,nuc(i)))
        j=j+1
       fato=.true.
       fato2=.true.
      k1=Jx(j)
      dd=d(Nuc(i),Nuc(j))
      do 40 ni=1,ncont(i)
      do 40 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then   
       if (rexp.lt.rmaxs) then   
       if (fato) then
        do iki=1,6
             kknumd=kknumd+1
             kkind(kknumd)=i+iki-1+k1
         enddo
        fato=.false.
       endif
       else

       if (fato2) then
        do iki=1,6
             kknums=kknums+1
              if(kknumsmax.lt.kknums) stop '5'
             kkinds(kknums)=i+iki-1+k1
         enddo
        fato2=.false.
       endif

       endif

c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 41 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 41 nk=1,ncontd(k)
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
      roz=tj
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-roz*ss1s)/z2
c
      do 45 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
c
      do 45 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      term=t1*ps+t2*p1s
c
      f1=1.D0
      if (l1.eq.l2) then
       term=term+ta
       f1=sq3
      endif
c
      cc=ccoef/f1
      term=term*cc
      l12=Ll(l1)+l2
c      ii=i+l12-1
c
c      kk=ii+k1
       if (rexp.lt.rmaxs) then   
      kk=l12-1+kknumd-5
      id = (kk-1)*Md+k
      cool(id) = cool(id) + term
      else
      kk=l12-1+kknums-5
      id = (kk-1)*Md+k
      cools(id) = cools(id) + real(term)

      endif
      ix=ix+1

 45   continue
 41   continue

c---(ds|p)
       do 101 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 101 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      t3=(ss1s-roz*ss2s)/z2
c
      do 105 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      p1s=t1*ss1s+t2*ss2s
      t5=p1s/z2a
      p2s=t1*ss2s+t2*ss3s
c
      do 105 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pj1s=t1*ss1s+t2*ss2s
      t4=pj1s/z2a
      ds=t1*p1s+t2*p2s
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+t3
       f1=sq3
      endif
c
      do 105 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      term=t1*ds

      if (l1.eq.l3) then
       term=term+t4
      endif
c
      if (l2.eq.l3) then
       term=term+t5
      endif
c
      l12=Ll(l1)+l2
c      ii=i+l12-1
      kk=k+l3-1
c 
      cc=ccoef/f1
      term=term*cc
c
c      kn=ii+k1
       if (rexp.lt.rmaxs) then
      kn=l12-1+kknumd-5
      id = (kn-1)*Md+kk
      cool(id) = cool(id) + term
      else
      kn=l12-1+kknums-5
      id = (kn-1)*Md+kk
      cools(id) = cools(id) + real(term)
      endif


      ix=ix+1
 105   continue
 101   continue
c------(ds|d)
      do 161 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 161 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
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
      do 165 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t7=p2s/z2a
c
      do 165 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      t8=pj2s/z2a
      f1=1.
      if (l1.eq.l2) then
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t5
       f1=sq3
      endif
c
      t11=(ds-ti*d1s)/zc
      do 165 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      ds1p=t1*d2s
      pis1pk=t1*p2s
      pjs1pk=t1*pj2s
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
      do 165 l4=1,l3
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
c      ii=i+l12-1
      kk=k+l34-1
c
      cc=ccoef/(f1*f2)
      term=term*cc
c      kn=ii+k1
      if (rexp.lt.rmaxs) then
      kn=l12-1+kknumd-5
      id = (kn-1)*Md+kk
      cool(id) = cool(id) + term
      else
      kn=l12-1+kknums-5
      id = (kn-1)*Md+kk
      cools(id) = cools(id) + real(term)

      endif


      ix=ix+1
 165   continue
 161   continue


 42   continue
      endif
 40   continue
c
c-------------------------------------------------------------
c (dp|s)
      do 50 i=ns+np+1,M,6
         do 50 knan=1,natomc(nuc(i))
         j=nnpp(jatc(knan,nuc(i)))-3
      do 50 kknan=1,nnp(jatc(knan,nuc(i))),3
        j=j+3
       fato=.true.
       fato2=.true.
      dd=d(Nuc(i),Nuc(j))
      do 50 ni=1,ncont(i)
      do 50 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.lt.rmax) then   
       if (rexp.lt.rmaxs) then   
       if (fato) then
        do iki=1,6
          do  jki=1,3
             kknumd=kknumd+1
             kkind(kknumd)=i+iki-1+Jx(j+jki-1)

            enddo
         enddo
        fato=.false.
       endif

       else
       if (fato2) then
        do iki=1,6
          do  jki=1,3
             kknums=kknums+1
              if(kknumsmax.lt.kknums) stop '6'
             kkinds(kknums)=i+iki-1+Jx(j+jki-1)
            enddo
         enddo
        fato2=.false.
       endif



       endif
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 51 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 51 nk=1,ncontd(k)
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
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      ii=0
c
      do 55 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      t5=(ps-roz*p1s)/z2
      p2s=t1*ss2s+t2*ss3s
c
      do 55 l2=1,l1
c
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      t6=(pjs-roz*pj1s)/z2
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
c
      f1=1.D0
      if (l1.eq.l2) then
       f1=sq3
       ds=ds+t3
       d1s=d1s+t4
      endif
c
      do 55 l3=1,3
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
c      ii=i+l12-1
c      jj=j+l3-1
c
       ii=ii+1

      cc=ccoef/f1
      term=term*cc
c
c      kk=ii+Jx(jj)
       if (rexp.lt.rmaxs) then   
       kk=ii+kknumd-18
      id =(kk-1)*Md+k
      cool(id) = cool(id) + term
      else
       kk=ii+kknums-18
      id =(kk-1)*Md+k
      cools(id) = cools(id) + real(term)

      endif
      ix=ix+1
 55   continue
 51   continue
c-----(dp|p)

      do 111 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 111 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(ss1s-roz*ss2s)/z2
      t4=(ss2s-roz*ss3s)/z2
      ii=0
c
      do 115 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      t5=(p1s-roz*p2s)/z2
      p3s=t1*ss3s+t2*ss4s
c
      do 115 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      t6=(pj1s-roz*pj2s)/z2
c
      f1=1.D0
      if (l1.eq.l2) then
       d1s=d1s+t3
       d2s=d2s+t4
       f1=sq3
      endif
      t9=d1s/z2a
c
      do 115 l3=1,3
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
      ii=ii+1
      do 115 l4=1,3
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
c      ii=i+l12-1
c      jj=j+l3-1
      kk=k+l4-1
c
      cc=ccoef/f1
      term=term*cc
c
c      kn=ii+Jx(jj)
       if (rexp.lt.rmaxs) then   
       kn=ii+kknumd-18
      id =(kn-1)*Md+kk
      cool(id) = cool(id) + term
      else
       kn=ii+kknums-18
      id =(kn-1)*Md+kk
      cools(id) = cools(id) + real(term)
      endif
      ix=ix+1
 115   continue
 111   continue

c-----(dp|d)

      do 171 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 171 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      ii=0
c
      do 175 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
c
      do 175 l2=1,l1
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
      f1=1.
      if (l1.eq.l2) then
       ds=ds+(sss-roz*ss1s)/z2
       d1s=d1s+(ss1s-roz*ss2s)/z2
       d2s=d2s+(ss2s-roz*ss3s)/z2
       d3s=d3s+(ss3s-roz*ss4s)/z2
       f1=sq3
      endif
c
      do 175 l3=1,3
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
ci
      ii=ii+1
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
      do 175 l4=1,3
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
      do 175 l5=1,l4
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
c      ii=i+l12-1
c      jj=j+l3-1
      l45=l4*(l4-1)/2+l5
      kk=k+l45-1
c
      cc=ccoef/(f1*f2)
      term=term*cc
c
c
c      kn=ii+Jx(jj)
       if (rexp.lt.rmaxs) then   
       kn=ii+kknumd-18
      id =(kn-1)*Md+kk
      cool(id) = cool(id) + term
       else
       kn=ii+kknums-18
      id =(kn-1)*Md+kk
      cools(id) = cools(id) + real(term)
       endif
      ix=ix+1
 175   continue
 171   continue

      endif
 50   continue
c
c-------------------------------------------------------------
c
c (dd|s)
      do 666 i=ns+np+1,M,6
         do 666 knan=1,natomc(nuc(i))
         j=nnpd(jatc(knan,nuc(i)))-6
      do 666 kknan=1,nnd(jatc(knan,nuc(i))),6
        j=j+6
       if(j.le.i) then
      fato=.true.
      fato2=.true.
      ddij=d(Nuc(i),Nuc(j))
      do 60 ni=1,ncont(i)
      do 60 nj=1,ncont(j)
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*ddij
       if (rexp.lt.rmax) then   
       if (rexp.lt.rmaxs) then   
       if (fato) then
         fato=.false.
          if(i.eq.j) then
             do iki=1,6
               do jki=1,iki
                 kknumd=kknumd+1
                 ii=i+iki-1
                 jj=j+jki-1
                 kkind(kknumd)=ii+Jx(jj)
               enddo
             enddo
          else
            do iki=1,6
              do jki=1,6
                kknumd=kknumd+1
                kkind(kknumd)=i+iki-1+Jx(j+jki-1)

              enddo
             enddo
         endif
       endif
       else
       if (fato2) then
         fato2=.false.
          if(i.eq.j) then
             do iki=1,6
               do jki=1,iki
                 kknums=kknums+1
              if(kknumsmax.lt.kknums) stop 'a'
                 ii=i+iki-1
                 jj=j+jki-1
                 kkinds(kknums)=ii+Jx(jj)
               enddo
             enddo
          else
            do iki=1,6
              do jki=1,6
                kknums=kknums+1
              if(kknumsmax.lt.kknums) stop 'b'
                kkinds(kknums)=i+iki-1+Jx(j+jki-1)

              enddo
             enddo
         endif
       endif


       endif
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 61 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 61 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      ii=0
c
      do 65 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t6=(ps-roz*p1s)/z2
      t7=(p1s-roz*p2s)/z2
c
      do 65 l2=1,l1
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
      f1=1.D0
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
      do 65 l3=1,lij
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
      do 65 l4=1,lk
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
c      ii=i+l12-1
c      jj=j+l34-1
c
       ii=ii+1
      cc=ccoef/(f1*f2)
      term=term*cc
c
c      kk=ii+Jx(jj)
       if (rexp.lt.rmaxs) then   
       if(j.eq.i) then
       kk=ii+kknumd-21
       else
       kk=ii+kknumd-36
        endif
      id = (kk-1)*Md+k

      cool(id) = cool(id) + term
      else
       if(j.eq.i) then
       kk=ii+kknums-21
       else
       kk=ii+kknums-36
        endif
      id = (kk-1)*Md+k

      cools(id) = cools(id) + real(term)

      endif
      ix=ix+1
 65   continue
 61   continue

c----(dd|p)

      do 121 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 121 nk=1,ncontd(k)
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
      t1=ad(k,nk)*dsqrt(t0)
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
       ii=0
c
      do 125 l1=1,3
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
      do 125 l2=1,l1
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
      do 125 l3=1,lij
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
c      write(*,*) lk,l3,Ll(l1),ll(l3),l1,l2 
c
      t16=(pj1pk-roz*pj2pk)/z2
      t17=(pi1pk-roz*pi2pk)/z2
c
      do 125 l4=1,lk
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
       ii=ii+1
      do 125 l5=1,3
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
      if (l3.eq.l5) then
       term=term+t13
      endif
c
      if (l4.eq.l5) then
       term=term+t14
      endif
c
cx      l12=Ll(l1)+l2
c      l34=Ll(l3)+l4
c      ii=i+l12-1
c      jj=j+l34-1
      kk=k+l5-1
     
c
      cc=ccoef/(f1*f2)
      term=term*cc
c
c      kn=ii+Jx(jj)
       if (rexp.lt.rmaxs) then   
       if(j.eq.i) then
       kn=ii+kknumd-21
       else
       kn=ii+kknumd-36
        endif
      id = (kn-1)*Md+kk

      cool(id) = cool(id) + term
      else
       if(j.eq.i) then
       kn=ii+kknums-21
       else
       kn=ii+kknums-36
        endif
      id = (kn-1)*Md+kk

      cools(id) = cools(id) + real(term)
      endif
 
      ix=ix+1
 125   continue
 121   continue

c-----(dd|d)

      do 181 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 181 nk=1,ncontd(k)
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
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*dsqrt(t0)
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
       ii=0
c
      do 185 l1=1,3
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
      do 185 l2=1,l1
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
      do 185 l3=1,lij
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
      do 185 l4=1,lk
       ii=ii+1
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
      t30=pjdkl/z2a
      t40=pidkl/z2a
      t50=sdkl/z2a
      t60=pj2pl/z2a
      t70=pi2pl/z2a
      t80=d2pl/z2a
      t23=(dd-ti*d1d)/zc
c
      do 185 l5=1,3
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
      do 185 l6=1,l5
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
c---------
      cc=ccoef/(f1*f2*f3)
      term=term*cc
c
c      l12=Ll(l1)+l2
c      l34=Ll(l3)+l4
      l56=Ll(l5)+l6
c      ii=i+l12-1
c      jj=j+l34-1
      kk=k+l56-1
c
c      kn=ii+Jx(jj)
       if (rexp.lt.rmaxs) then
       if(j.eq.i) then
       kn=ii+kknumd-21
       else
       kn=ii+kknumd-36
        endif
      id = (kn-1)*Md+kk
      cool(id) = cool(id) + term

      else
       if(j.eq.i) then
       kn=ii+kknums-21
       else
       kn=ii+kknums-36
        endif
      id = (kn-1)*Md+kk
      cools(id) = cools(id) + real(term)


      endif
c      write(98,*) cool(id),id
      ix=ix+1
 185  continue
 181  continue

      endif
 60   continue
      endif
666   continue
c
c-------------------------------------------------------------
c       write(*,*) 'kknum',kknumd,MM
*
c------------------------------------------------------------------
c-------------------------------------------------------------
*
      deallocate (Jx)
      return
      end subroutine
      end module subm_int3mem
