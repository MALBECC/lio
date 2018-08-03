module subm_int3mem
contains
subroutine int3mem()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Integrals subroutine - 2e integrals with 3 indexes : wavefunction and density!
! fitting functions, calculated using the Obara-Saika recursive method.        !
! Input : G ,F,  standard basis and density basis                              !
! F has already the 1e terms, and here the Coulomb part is added without       !
! storing the integrals separately.                                            !
! Output: F updated with Coulomb part and Coulomb energy is calculated.        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use liotemp   , only: FUNCT
   use garcha_mod, only: cool, cools, kkind, kkinds, Nuc, Nucd, a, c, &
                         d, r, ad, cd, natomc, nns, nnp, nnd, nnps, nnpp, nnpd,&
                         jatc, ncont, ncontd, nshell, nshelld, M, Md, rmax,    &
                         rmaxs, pi52, NORM, kknums, kknumd
   implicit none
   double precision  :: Q(3), W(3)
   integer, dimension(:), allocatable :: Jx

   ! Eliminating implicits:
   double precision  :: alf, cc, ccoef, f1, f2, f3
   double precision  :: rexp, ro, roz, sq3, term, u, ddij
   double precision  :: tii, tjj
   double precision  :: z2, z2a, zc, Zij

   double precision  :: d1s, d1p, d1d, d1pk, d1pl
   double precision  :: d2s, d2p, d2d, d2pl, d2pk
   double precision  :: d3s, d3pk, d4s
   double precision  :: ds, ds1p, dspl, dp, dpc, dpk, dp1p
   double precision  :: dd, ddp, dijplp, dijpkp

   double precision  :: ps, pp, pp1p, p1s, p2s, p3s, p4s, p5s
   double precision  :: pi1p, pi1pk, pi2p, pi2pk, pi2pl, pi3pk
   double precision  :: pijs, pij1s, pij2s, pispj, pispk, pis1pk
   double precision  :: pip, pipk, pipkpl, pidkl, pidklp
   double precision  :: pjs, pjs1pk
   double precision  :: pj1s, pj1p, pj1pk, pj2s, pj2p, pj2pk, pj2pl
   double precision  :: pj3s, pj3pk, pj4s
   double precision  :: pjp, pjpk, pjpkpl, pjdklp,pjdkl

   double precision  :: sss, sks, spk, spjs, sspj, spjpk, sp2js
   double precision  :: ss1s, ss2s, ss3s, ss4s, ss5s, ss6s
   double precision  :: ss1p, s1pk, s2pk, s3pk, spks, spj, sdkl

   double precision  :: ta, ti, tj
   double precision  :: t0, t1, t2, t3, t4, t5, t6, t6b, t7, t8, t9
   double precision  :: t10, t11, t12, t13, t14, t15, t16, t17, t18
   double precision  :: t20, t21, t22, t22a, t23, t24, t25, t26
   double precision  :: t27, t28, t29, t30, t31
   double precision  :: t40, t41, t50, t51, t60, t61, t70, t80

   integer :: M2, Ll(3)
   integer :: i, ii, j, jj, k, kk, ni, nj, nk, kn, k1
   integer :: id, iki, jki, ix, kknan, knan, kknumsmax
   integer :: ns, nsd, nd, ndd, np, npd
   integer :: l, lk, l1, l2, l3, l4, l5, l6
   integer :: lij, l12, l23, l34, l45, l56

   logical :: fato, fato2

   allocate (Jx(M))
   ns  = nshell(0) ; np  = nshell(1) ; nd  = nshell(2)
   nsd = nshelld(0); npd = nshelld(1); ndd = nshelld(2)

   M2=2*M
   Md2=2*Md

   sq3=1.D0
   if (NORM) sq3 = dsqrt(3.D0)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1-1) / 2
   enddo
   do ifunct = 1, M
      Jx(ifunct) = (M2-ifunct) * (ifunct-1) / 2
   enddo

   ix     = 0
   kknumd = 0
   kknums = 0

   ! Search for the dimensions of cool/cools.
   do ifunct = 1, ns
   do knan   = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         if (jfunct .le. ifunct) then

            kk = ifunct + Jx(jfunct)
            dd = d(Nuc(ifunct),Nuc(jfunct))
            fato  = .true.
            fato2 = .true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               rexp = a(ifunct,nci) * a(jfunct,ncj) * dd / &
                      (a(ifunct,nci) + a(jfunct,ncj))

               if (rexp.lt.rmax) then
               if (rexp.lt.rmaxs) then
                  if (fato) then
                     kknumd = kknumd +1
                     fato   = .false.
                  endif
               else
                  if (fato2) then
                     kknums = kknums +1
                     fato2  = .false.
                  endif
               endif
               endif
            enddo
            enddo
         endif
      enddo
   enddo
   enddo

   do ifunct = ns+1, ns+np, 3
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         fato  = .true.
         fato2 = .true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij  = a(ifunct,nci) + a(jfunct,ncj)
            tj   = a(jfunct,ncj) / Zij
            rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

            if (rexp.lt.rmax) then
            if (rexp.lt.rmaxs) then
               if (fato) then
                  do iki=1,3
                     kknumd=kknumd+1
                  enddo
                  fato   = .false.
               endif
            else
               if (fato2) then
                  do iki=1,3
                     kknums=kknums+1
                  enddo
                  fato2  = .false.
               endif
            endif
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   do ifunct = ns+1, ns+np, 3
   do knan = 1, natomc(Nuc(ifunct))
      j=nnpp(jatc(knan,Nuc(i)))-3
      do kknan=1,nnp(jatc(knan,Nuc(i))),3
         j=j+3
         if (jfunct .le. ifunct) then
            fato=.true.
            fato2=.true.
            dd = d(Nuc(ifunct),Nuc(jfunct))
            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij = a(ifunct,nci) + a(jfunct,ncj)
               z2=2.D0*Zij
               ti = a(ifunct,nci) / Zij
               tj = a(jfunct,ncj) / Zij
               alf = a(ifunct,nci) * tj
               rexp = alf * dd
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
                     fato   = .false.
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
                     fato2  = .false.
                  endif
               endif
               endif
            enddo
            enddo
         endif
      enddo
   enddo
   enddo

   do i=ns+np+1,M,6
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1
      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         fato=.true.
         fato2=.true.
         k1=Jx(jfunct)
         dd = d(Nuc(ifunct),Nuc(jfunct))
         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij = a(ifunct,nci) + a(jfunct,ncj)
            z2=2.D0*Zij
            ti = a(ifunct,nci) / Zij
            tj = a(jfunct,ncj) / Zij
            alf = a(ifunct,nci) * tj
            rexp = alf * dd
            if (rexp.lt.rmax) then
            if (rexp.lt.rmaxs) then
               if (fato) then
                  do iki=1,6
                     kknumd=kknumd+1
                  enddo
                  fato   = .false.
               endif
            else
               if (fato2) then
                  do iki=1,6
                     kknums=kknums+1
                  enddo
                  fato2  = .false.
               endif
            endif
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   do i=ns+np+1,M,6
   do knan = 1, natomc(Nuc(ifunct))
      j=nnpp(jatc(knan,Nuc(i)))-3
      do kknan=1,nnp(jatc(knan,Nuc(i))),3
         j=j+3
         fato=.true.
         fato2=.true.
         dd = d(Nuc(ifunct),Nuc(jfunct))
         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij = a(ifunct,nci) + a(jfunct,ncj)
            z2=2.D0*Zij
            ti = a(ifunct,nci) / Zij
            tj = a(jfunct,ncj) / Zij
            alf = a(ifunct,nci) * tj
            rexp = alf * dd
            if (rexp.lt.rmax) then
            if (rexp.lt.rmaxs) then
               if (fato) then
                  do iki=1,6
                  do  jki=1,3
                     kknumd=kknumd+1
                  enddo
                  enddo
                  fato   = .false.
               endif
            else
               if (fato2) then
                  do iki=1,6
                  do  jki=1,3
                     kknums=kknums+1
                  enddo
                  enddo
                  fato2  = .false.
               endif
            endif
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   do i=ns+np+1,M,6
   do knan = 1, natomc(Nuc(ifunct))
      j=nnpd(jatc(knan,Nuc(i)))-6
      do kknan=1,nnd(jatc(knan,Nuc(i))),6
         j=j+6
         if (jfunct .le. ifunct) then
            fato=.true.
            fato2=.true.
            dd = d(Nuc(ifunct),Nuc(jfunct))
            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij = a(ifunct,nci) + a(jfunct,ncj)
               z2=2.D0*Zij
               ti = a(ifunct,nci) / Zij
               tj = a(jfunct,ncj) / Zij
               alf = a(ifunct,nci) * tj
               rexp = alf * dd
               if (rexp.lt.rmax) then
               if (rexp.lt.rmaxs) then
                  if (fato) then
                     fato   = .false.
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
                     fato2  = .false.
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
            enddo
            enddo
         endif
      enddo
   enddo
   enddo

   if (allocated(cool))   deallocate(cool)
   if (allocated(kkind))  deallocate(kkind)
   if (allocated(cools))  deallocate(cools)
   if (allocated(kkinds)) deallocate(kkinds)
   allocate(cool(kknumd*Md), cools(kknums*Md))
   allocate(kkind(kknumd)  , kkinds(kknums))
   ! End of cool dimensions

   kknumsmax=kknums
   cool=0
   cools=0
   kknumd=0
   kknums=0

   ! Start of integrals.
   ! (ss|s)
   do ifunct = 1, ns
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         if (jfunct .le. ifunct) then

            kk=i+Jx(jfunct)
            dd = d(Nuc(ifunct),Nuc(jfunct))

            fato=.true.
            fato2=.true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij = a(ifunct,nci) + a(jfunct,ncj)
               ti = a(ifunct,nci) / Zij
               tj = a(jfunct,ncj) / Zij
               alf = a(ifunct,nci) * tj
               rexp = alf * dd
               if (rexp.lt.rmax) then
                  if (rexp.lt.rmaxs) then
                     if (fato) then
                        kknumd = kknumd+1
                        kkind(kknumd)=kk
                        fato   = .false.
                     endif
                  else
                     if (fato2) then
                        kknums = kknums+1
                        if(kknumsmax.lt.kknums) stop '1'
                        kkinds(kknums)=kk
                        fato2  = .false.
                     endif
                  endif

                  Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
                  Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
                  Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)

                  sks=pi52*exp(-rexp)/Zij

                  do k=1,nsd
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij
                        u=ad(k,nk)*Zij/t0*dpc
                        t1=ad(k,nk)*dsqrt(t0)
                        term=sks/t1*FUNCT(0,u)
                        term=term*ccoef

                        if (rexp.lt.rmaxs) then
                           id = (kknumd-1)*md+k
                           cool(id) = cool(id) + term
                        else
                           id = (kknums-1)*md+k
                           cools(id) = cools(id) + real(term)
                        endif

                        ix=ix+1
                     enddo
                  enddo

                  !----(ss|p)
                  do k=nsd+1,nsd+npd,3
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+ &
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij

                        tii=Zij/t0
                        tjj=ad(k,nk)/t0
                        W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
                        W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
                        W(3)=tii*Q(3)+tjj*r(Nucd(k),3)

                        u=ad(k,nk)*tii*dpc
                        t1=ad(k,nk)*dsqrt(t0)
                        t2=sks/t1
                        sss=t2*FUNCT(0,u)
                        ss1s=t2*FUNCT(1,u)

                        do l1=1,3
                           t1=W(l1)-r(Nucd(k),l1)
                           term=ccoef*t1*ss1s
                           kn=k+l1-1
                           if (rexp.lt.rmaxs) then
                              id = (kknumd-1)*Md+kn
                              cool(id) = cool(id) + term
                           else
                              id = (kknums-1)*Md+kn
                              cools(id) = cools(id) + real(term)
                           endif

                           ix=ix+1
                        enddo
                     enddo
                  enddo


                  !---- (ss|d)
                  do k=nsd+npd+1,Md,6
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij
                        tii=Zij/t0
                        tjj=ad(k,nk)/t0
                        W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
                        W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
                        W(3)=tii*Q(3)+tjj*r(Nucd(k),3)

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

                        do l1=1,3
                           t1=W(l1)-r(Nucd(k),l1)
                           ss1p=t1*ss2s

                           do l2=1,l1
                              t1=W(l2)-r(Nucd(k),l2)
                              term=t1*ss1p

                              f1=1.D0
                              if (l1.eq.l2) then
                                 term=term+ta
                                 f1=sq3
                              endif

                              cc=ccoef/f1
                              term=term*cc

                              l12=Ll(l1)+l2
                              kn=k+l12-1

                              if (rexp.lt.rmaxs) then
                                 id = (kknumd-1)*Md+kn
                                 cool(id) = cool(id) + term
                              else
                                 id = (kknums-1)*Md+kn
                                 cools(id) = cools(id) + real(term)
                              endif

                              ix=ix+1
                           enddo
                        enddo
                     enddo
                  enddo

               endif
            enddo
            enddo
         endif
      enddo
   enddo
   enddo

   ! (ps|s)
   do ifunct = ns+1, ns+np, 3
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         dd = d(Nuc(ifunct),Nuc(jfunct))
         k1=Jx(jfunct)
         fato=.true.
         fato2=.true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij = a(ifunct,nci) + a(jfunct,ncj)
            ti = a(ifunct,nci) / Zij
            tj = a(jfunct,ncj) / Zij
            alf = a(ifunct,nci) * tj
            rexp = alf * dd
            if (rexp.lt.rmax) then
               if (rexp.lt.rmaxs) then
                  if (fato) then
                     do iki=1,3
                        kknumd=kknumd+1
                        kkind(kknumd)=i+k1+iki-1
                     enddo
                     fato   = .false.
                  endif
               else
                  if (fato2) then
                     do iki=1,3
                        kknums=kknums+1
                        if(kknumsmax.lt.kknums) stop '2'
                        kkinds(kknums)=i+k1+iki-1
                     enddo
                     fato2  = .false.
                  endif
               endif

               Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
               Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
               Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)

               sks=pi52*exp(-rexp)/Zij

               do k=1,nsd
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij

                     tii=Zij/t0
                     tjj=ad(k,nk)/t0
                     W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
                     W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
                     W(3)=tii*Q(3)+tjj*r(Nucd(k),3)


                     u=ad(k,nk)*tii*dpc
                     t1=ad(k,nk)*dsqrt(t0)
                     t2=sks/t1

                     sss=t2*FUNCT(0,u)
                     ss1s=t2*FUNCT(1,u)

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        term=ccoef*(t1*sss+t2*ss1s)
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
                     enddo
                  enddo
               enddo
               !------(ps|p)
               do k=nsd+1,nsd+npd,3
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij
                     z2a=2.*t0

                     tii=Zij/t0
                     tjj=ad(k,nk)/t0
                     W(1)=tii*Q(1)+tjj*r(Nucd(k),1)
                     W(2)=tii*Q(2)+tjj*r(Nucd(k),2)
                     W(3)=tii*Q(3)+tjj*r(Nucd(k),3)

                     u=ad(k,nk)*tii*dpc
                     t1=ad(k,nk)*dsqrt(t0)
                     t2=sks/t1
                     sss=t2*FUNCT(0,u)
                     ss1s=t2*FUNCT(1,u)
                     ss2s=t2*FUNCT(2,u)
                     ta=ss1s/z2a

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        p1s=t1*ss1s+t2*ss2s

                        do l2=1,3
                           t2=W(l2)-r(Nucd(k),l2)
                           term=t2*p1s

                           if (l1.eq.l2) then
                              term=term+ta
                           endif
                           kk=k+l2-1
                           term=term*ccoef

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
                        enddo
                     enddo
                  enddo
               enddo

               !-------(ps|d)
               do k=nsd+npd+1,Md,6
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij
                     zc=2.D0*ad(k,nk)
                     z2a=2.D0*t0

                     ti=Zij/t0
                     tj=ad(k,nk)/t0
                     W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                     W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                     W(3)=ti*Q(3)+tj*r(Nucd(k),3)

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

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        ps=t1*sss+t2*ss1s
                        p1s=t1*ss1s+t2*ss2s
                        p2s=t1*ss2s+t2*ss3s
                        t5=(ps-roz*p1s)/zc

                        do l2=1,3
                           t1=W(l2)-r(Nucd(k),l2)
                           sspj=t1*ss2s
                           pispj=t1*p2s

                           t4=sspj/z2a
                           if (l1.eq.l2) then
                              pispj=pispj+t3
                           endif

                           do l3=1,l2
                              t1=W(l3)-r(Nucd(k),l3)
                              term=t1*pispj

                              if (l1.eq.l3) then
                                 term=term+t4
                              endif

                              f1=1.D0
                              if (l2.eq.l3) then
                                 term=term+t5
                                 f1=sq3
                              endif

                              l23=l2*(l2-1)/2+l3

                              kk=k+l23-1

                              cc=ccoef/f1
                              term=term*cc
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
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   ! (pp|s)
   do ifunct = ns+1, ns+np, 3
   do knan = 1, natomc(Nuc(ifunct))
      j=nnpp(jatc(knan,Nuc(i)))-3

      do kknan=1,nnp(jatc(knan,Nuc(i))),3
         j=j+3

         if (jfunct .le. ifunct) then
            fato=.true.
            fato2=.true.
            dd = d(Nuc(ifunct),Nuc(jfunct))
            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij = a(ifunct,nci) + a(jfunct,ncj)
               z2=2.D0*Zij
               ti = a(ifunct,nci) / Zij
               tj = a(jfunct,ncj) / Zij
               alf = a(ifunct,nci) * tj
               rexp = alf * dd
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
                        fato   = .false.
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
                        fato2  = .false.
                     endif
                  endif

                  Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
                  Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
                  Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
                  sks=pi52*exp(-rexp)/Zij

                  do k=1,nsd
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij

                        ti=Zij/t0
                        tj=ad(k,nk)/t0
                        W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                        W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                        W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                        ro=ad(k,nk)*ti
                        u=ro*dpc
                        t1=ad(k,nk)*dsqrt(t0)
                        t2=sks/t1
                        sss=t2*FUNCT(0,u)
                        ss1s=t2*FUNCT(1,u)
                        ss2s=t2*FUNCT(2,u)
                        ta=(sss-tj*ss1s)/z2

                        ii=0
                        do l1=1,3
                           t1=Q(l1)-r(Nuc(i),l1)
                           t2=W(l1)-Q(l1)
                           ps=t1*sss+t2*ss1s
                           p1s=t1*ss1s+t2*ss2s

                           lij=3
                           if (i.eq.j) then
                              lij=l1
                           endif

                           do l2=1,lij
                              t1=Q(l2)-r(Nuc(j),l2)
                              t2=W(l2)-Q(l2)
                              term=t1*ps+t2*p1s

                              if (l1.eq.l2) then
                                 term=term+ta
                              endif

                              term=term*ccoef
                              ii=ii+1
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
                           enddo
                        enddo
                     enddo
                  enddo

                  ! (pp|p)
                  do k=nsd+1,nsd+npd,3
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij
                        z2a=2.*t0

                        ti=Zij/t0
                        tj=ad(k,nk)/t0
                        W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                        W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                        W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                        roz=tj
                        ro=roz*Zij
                        u=ro*dpc
                        t1=ad(k,nk)*dsqrt(t0)
                        t2=sks/t1
                        sss=t2*FUNCT(0,u)
                        ss1s=t2*FUNCT(1,u)
                        ss2s=t2*FUNCT(2,u)
                        ss3s=t2*FUNCT(3,u)
                        t3=(ss1s-roz*ss2s)/z2

                        ii=0
                        do l1=1,3
                           t1=Q(l1)-r(Nuc(i),l1)
                           t2=W(l1)-Q(l1)
                           p1s=t1*ss1s+t2*ss2s
                           p2s=t1*ss2s+t2*ss3s
                           t5=p1s/z2a

                           lij=3
                           if (i.eq.j) then
                              lij=l1
                           endif

                           do l2=1,lij
                              t1=Q(l2)-r(Nuc(j),l2)
                              t2=W(l2)-Q(l2)
                              spj=t1*ss1s+t2*ss2s
                              t4=spj/z2a
                              pp=t1*p1s+t2*p2s

                              if (l1.eq.l2) then
                                 pp=pp+t3
                              endif

                              ii=ii+1
                              do l3=1,3
                                 t1=W(l3)-r(Nucd(k),l3)
                                 term=t1*pp

                                 if (l1.eq.l3) then
                                    term=term+t4
                                 endif

                                 if (l2.eq.l3) then
                                    term=term+t5
                                 endif
                                 kk=k+l3-1

                                 term=term*ccoef
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
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo

                  !(pp|d)
                  do k=nsd+npd+1,Md,6
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij
                        z2a=2.*t0

                        ti=Zij/t0
                        tj=ad(k,nk)/t0
                        W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                        W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                        W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                        roz=tj

                        ro=roz*Zij
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

                        do l1=1,3
                           t1=Q(l1)-r(Nuc(i),l1)
                           t2=W(l1)-Q(l1)
                           ps=t1*sss+t2*ss1s
                           p1s=t1*ss1s+t2*ss2s
                           p2s=t1*ss2s+t2*ss3s
                           p3s=t1*ss3s+t2*ss4s
                           t8=p2s/z2a

                           lij=3
                           if (i.eq.j) then
                              lij=l1
                           endif

                           do l2=1,lij
                              t1=Q(l2)-r(Nuc(j),l2)
                              t2=W(l2)-Q(l2)
                              pijs=t1*ps+t2*p1s
                              pij1s=t1*p1s+t2*p2s
                              pij2s=t1*p2s+t2*p3s
                              spjs=t1*ss1s+t2*ss2s
                              sp2js=t1*ss2s+t2*ss3s
                              t7=sp2js/z2a

                              ii=ii+1

                              if (l1.eq.l2) then
                                 pijs=pijs+t3
                                 pij1s=pij1s+t4
                                 pij2s=pij2s+t5
                              endif

                              t11=(pijs-ti*pij1s)/zc

                              do l3=1,3
                                 t1=W(l3)-r(Nucd(k),l3)
                                 pp1p=t1*pij2s
                                 spjpk=t1*sp2js
                                 pispk=t1*p2s

                                 if (l1.eq.l3) then
                                    pp1p=pp1p+t7
                                    pispk=pispk+t6
                                 endif

                                 if (l2.eq.l3) then
                                    pp1p=pp1p+t8
                                    spjpk=spjpk+t6
                                 endif

                                 t9=spjpk/z2a
                                 t10=pispk/z2a

                                 do l4=1,l3
                                    t1=W(l4)-r(Nucd(k),l4)
                                    term=t1*pp1p

                                    if (l1.eq.l4) then
                                       term=term+t9
                                    endif

                                    if (l2.eq.l4) then
                                       term=term+t10
                                    endif

                                    f1=1.D0
                                    if (l3.eq.l4) then
                                       term=term+t11
                                       f1=sq3
                                    endif
                                    l34=l3*(l3-1)/2+l4
                                    kk=k+l34-1

                                    cc=ccoef/f1
                                    term=term*cc

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

                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               endif
            enddo
            enddo
         endif
      enddo
   enddo
   enddo

   ! (ds|s)
   do i=ns+np+1,M,6
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         fato=.true.
         fato2=.true.
         k1=Jx(jfunct)
         dd = d(Nuc(ifunct),Nuc(jfunct))
         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij = a(ifunct,nci) + a(jfunct,ncj)
            z2=2.D0*Zij
            ti = a(ifunct,nci) / Zij
            tj = a(jfunct,ncj) / Zij
            alf = a(ifunct,nci) * tj
            rexp = alf * dd
            if (rexp.lt.rmax) then
               if (rexp.lt.rmaxs) then
                  if (fato) then
                     do iki=1,6
                        kknumd=kknumd+1
                        kkind(kknumd)=i+iki-1+k1
                     enddo
                     fato   = .false.
                  endif
               else

                  if (fato2) then
                     do iki=1,6
                        kknums=kknums+1
                        if(kknumsmax.lt.kknums) stop '5'
                        kkinds(kknums)=i+iki-1+k1
                     enddo
                     fato2  = .false.
                  endif
               endif

               Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
               Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
               Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)

               sks=pi52*exp(-rexp)/Zij

               do k=1,nsd
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij

                     ti=Zij/t0
                     tj=ad(k,nk)/t0
                     W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                     W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                     W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                     roz=tj
                     ro=roz*Zij
                     u=ro*dpc
                     t1=ad(k,nk)*dsqrt(t0)
                     t2=sks/t1
                     sss=t2*FUNCT(0,u)
                     ss1s=t2*FUNCT(1,u)
                     ss2s=t2*FUNCT(2,u)
                     ta=(sss-roz*ss1s)/z2

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        ps=t1*sss+t2*ss1s
                        p1s=t1*ss1s+t2*ss2s

                        do l2=1,l1
                           t1=Q(l2)-r(Nuc(i),l2)
                           t2=W(l2)-Q(l2)
                           term=t1*ps+t2*p1s

                           f1=1.D0
                           if (l1.eq.l2) then
                              term=term+ta
                              f1=sq3
                           endif

                           cc=ccoef/f1
                           term=term*cc
                           l12=Ll(l1)+l2
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
                        enddo
                     enddo
                  enddo
               enddo

               !---(ds|p)
               do k=nsd+1,nsd+npd,3
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij
                     z2a=2.*t0

                     ti=Zij/t0
                     tj=ad(k,nk)/t0
                     W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                     W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                     W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                     roz=tj
                     ro=roz*Zij
                     u=ro*dpc
                     t1=ad(k,nk)*dsqrt(t0)
                     t2=sks/t1
                     sss=t2*FUNCT(0,u)
                     ss1s=t2*FUNCT(1,u)
                     ss2s=t2*FUNCT(2,u)
                     ss3s=t2*FUNCT(3,u)
                     t3=(ss1s-roz*ss2s)/z2

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        p1s=t1*ss1s+t2*ss2s
                        t5=p1s/z2a
                        p2s=t1*ss2s+t2*ss3s

                        do l2=1,l1
                           t1=Q(l2)-r(Nuc(i),l2)
                           t2=W(l2)-Q(l2)
                           pj1s=t1*ss1s+t2*ss2s
                           t4=pj1s/z2a
                           ds=t1*p1s+t2*p2s

                           f1=1.D0
                           if (l1.eq.l2) then
                              ds=ds+t3
                              f1=sq3
                           endif

                           do l3=1,3
                              t1=W(l3)-r(Nucd(k),l3)
                              term=t1*ds

                              if (l1.eq.l3) then
                                 term=term+t4
                              endif

                              if (l2.eq.l3) then
                                 term=term+t5
                              endif

                              l12=Ll(l1)+l2
                              kk=k+l3-1
                              cc=ccoef/f1
                              term=term*cc
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
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

               !------(ds|d)
               do k=nsd+npd+1,Md,6
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij
                     z2a=2.D0*t0
                     zc=2.D0*ad(k,nk)

                     ti=Zij/t0
                     tj=ad(k,nk)/t0
                     W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                     W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                     W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                     roz=tj

                     ro=roz*Zij
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

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        ps=t1*sss+t2*ss1s
                        p1s=t1*ss1s+t2*ss2s
                        p2s=t1*ss2s+t2*ss3s
                        p3s=t1*ss3s+t2*ss4s
                        t7=p2s/z2a

                        do l2=1,l1
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

                           t11=(ds-ti*d1s)/zc
                           do l3=1,3
                              t1=W(l3)-r(Nucd(k),l3)
                              ds1p=t1*d2s
                              pis1pk=t1*p2s
                              pjs1pk=t1*pj2s

                              if (l1.eq.l3) then
                                 ds1p=ds1p+t8
                                 pis1pk=pis1pk+t6
                              endif

                              if (l2.eq.l3) then
                                 ds1p=ds1p+t7
                                 pjs1pk=pjs1pk+t6
                              endif

                              t9=pjs1pk/z2a
                              t10=pis1pk/z2a

                              do l4=1,l3
                                 t1=W(l4)-r(Nucd(k),l4)
                                 term=t1*ds1p

                                 if (l1.eq.l4) then
                                    term=term+t9
                                 endif

                                 if (l2.eq.l4) then
                                    term=term+t10
                                 endif

                                 f2=1.D0
                                 if (l3.eq.l4) then
                                    term=term+t11
                                    f2=sq3
                                 endif

                                 l12=Ll(l1)+l2
                                 l34=Ll(l3)+l4
                                 kk=k+l34-1

                                 cc=ccoef/(f1*f2)
                                 term=term*cc
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
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   ! (dp|s)
   do i=ns+np+1,M,6
   do knan = 1, natomc(Nuc(ifunct))
      j=nnpp(jatc(knan,Nuc(i)))-3
      do kknan=1,nnp(jatc(knan,Nuc(i))),3
         j=j+3
         fato=.true.
         fato2=.true.
         dd = d(Nuc(ifunct),Nuc(jfunct))
         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij = a(ifunct,nci) + a(jfunct,ncj)
            z2=2.D0*Zij
            ti = a(ifunct,nci) / Zij
            tj = a(jfunct,ncj) / Zij
            alf = a(ifunct,nci) * tj
            rexp = alf * dd
            if (rexp.lt.rmax) then
               if (rexp.lt.rmaxs) then
                  if (fato) then
                     do iki=1,6
                        do  jki=1,3
                           kknumd=kknumd+1
                           kkind(kknumd)=i+iki-1+Jx(j+jki-1)

                        enddo
                     enddo
                     fato   = .false.
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
                     fato2  = .false.
                  endif
               endif

               Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
               Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
               Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)

               sks=pi52*exp(-rexp)/Zij

               do k=1,nsd
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij

                     ti=Zij/t0
                     tj=ad(k,nk)/t0
                     W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                     W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                     W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                     roz=tj
                     ro=roz*Zij
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

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        ps=t1*sss+t2*ss1s
                        p1s=t1*ss1s+t2*ss2s
                        t5=(ps-roz*p1s)/z2
                        p2s=t1*ss2s+t2*ss3s

                        do l2=1,l1
                           t1=Q(l2)-r(Nuc(i),l2)
                           t2=W(l2)-Q(l2)
                           pjs=t1*sss+t2*ss1s
                           pj1s=t1*ss1s+t2*ss2s
                           t6=(pjs-roz*pj1s)/z2
                           ds=t1*ps+t2*p1s
                           d1s=t1*p1s+t2*p2s

                           f1=1.D0
                           if (l1.eq.l2) then
                              f1=sq3
                              ds=ds+t3
                              d1s=d1s+t4
                           endif

                           do l3=1,3
                              t1=Q(l3)-r(Nuc(j),l3)
                              t2=W(l3)-Q(l3)
                              term=t1*ds+t2*d1s

                              if (l1.eq.l3) then
                                 term=term+t6
                              endif

                              if (l2.eq.l3) then
                                 term=term+t5
                              endif

                              l12=Ll(l1)+l2
                              ii=ii+1

                              cc=ccoef/f1
                              term=term*cc
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
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

               !-----(dp|p)
               do k=nsd+1,nsd+npd,3
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij
                     z2a=2.*t0

                     ti=Zij/t0
                     tj=ad(k,nk)/t0
                     W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                     W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                     W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                     roz=tj

                     ro=roz*Zij
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

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        ps=t1*sss+t2*ss1s
                        p1s=t1*ss1s+t2*ss2s
                        p2s=t1*ss2s+t2*ss3s
                        t5=(p1s-roz*p2s)/z2
                        p3s=t1*ss3s+t2*ss4s

                        do l2=1,l1
                           t1=Q(l2)-r(Nuc(i),l2)
                           t2=W(l2)-Q(l2)
                           d1s=t1*p1s+t2*p2s
                           d2s=t1*p2s+t2*p3s
                           pj1s=t1*ss1s+t2*ss2s
                           pj2s=t1*ss2s+t2*ss3s
                           t6=(pj1s-roz*pj2s)/z2

                           f1=1.D0
                           if (l1.eq.l2) then
                              d1s=d1s+t3
                              d2s=d2s+t4
                              f1=sq3
                           endif
                           t9=d1s/z2a

                           do l3=1,3
                              t1=Q(l3)-r(Nuc(j),l3)
                              t2=W(l3)-Q(l3)
                              d1p=t1*d1s+t2*d2s
                              pi1p=t1*p1s+t2*p2s
                              pj1p=t1*pj1s+t2*pj2s

                              if (l1.eq.l3) then
                                 d1p=d1p+t6
                                 pi1p=pi1p+t3
                              endif

                              if (l2.eq.l3) then
                                 d1p=d1p+t5
                                 pj1p=pj1p+t3
                              endif

                              t7=pi1p/z2a
                              t8=pj1p/z2a

                              ii=ii+1
                              do l4=1,3
                                 t1=W(l4)-r(Nucd(k),l4)
                                 term=t1*d1p

                                 if (l1.eq.l4) then
                                    term=term+t8
                                 endif

                                 if (l2.eq.l4) then
                                    term=term+t7
                                 endif

                                 if (l3.eq.l4) then
                                    term=term+t9
                                 endif

                                 l12=Ll(l1)+l2
                                 kk=k+l4-1

                                 cc=ccoef/f1
                                 term=term*cc
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
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

               !-----(dp|d)
               do k=nsd+npd+1,Md,6
                  dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+ &
                  (Q(3)-r(Nucd(k),3))**2

                  do nk=1,ncontd(k)
                     ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                     t0=ad(k,nk)+Zij
                     z2a=2.D0*t0
                     zc=2.D0*ad(k,nk)

                     ti=Zij/t0
                     tj=ad(k,nk)/t0
                     W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                     W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                     W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                     roz=tj

                     ro=roz*Zij
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

                     do l1=1,3
                        t1=Q(l1)-r(Nuc(i),l1)
                        t2=W(l1)-Q(l1)
                        ps=t1*sss+t2*ss1s
                        p1s=t1*ss1s+t2*ss2s
                        p2s=t1*ss2s+t2*ss3s
                        p3s=t1*ss3s+t2*ss4s
                        p4s=t1*ss4s+t2*ss5s

                        do l2=1,l1
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

                           do l3=1,3
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

                              ii=ii+1
                              if (l1.eq.l3) then
                                 dp=dp+(pjs-roz*pj1s)/z2
                                 d1p=d1p+(pj1s-roz*pj2s)/z2
                                 d2p=d2p+(pj2s-roz*pj3s)/z2
                                 pi1p=pi1p+(ss1s-roz*ss2s)/z2
                                 pi2p=pi2p+(ss2s-roz*ss3s)/z2
                              endif

                              if (l2.eq.l3) then
                                 dp=dp+(ps-roz*p1s)/z2
                                 d1p=d1p+(p1s-roz*p2s)/z2
                                 d2p=d2p+(p2s-roz*p3s)/z2
                                 pj1p=pj1p+(ss1s-roz*ss2s)/z2
                                 pj2p=pj2p+(ss2s-roz*ss3s)/z2
                              endif

                              do l4=1,3
                                 t1=W(l4)-r(Nucd(k),l4)
                                 dp1p=t1*d2p
                                 pjpkpl=t1*pj2p
                                 pipkpl=t1*pi2p
                                 dspl=t1*d2s

                                 if (l1.eq.l4) then
                                    dp1p=dp1p+pj2p/z2a
                                    pipkpl=pipkpl+spks/z2a
                                    dspl=dspl+pj2s/z2a
                                 endif

                                 if (l2.eq.l4) then
                                    dp1p=dp1p+pi2p/z2a
                                    pjpkpl=pjpkpl+spks/z2a
                                    dspl=dspl+p2s/z2a
                                 endif

                                 if (l3.eq.l4) then
                                    dp1p=dp1p+d2s/z2a
                                    pipkpl=pipkpl+p2s/z2a
                                    pjpkpl=pjpkpl+pj2s/z2a
                                 endif

                                 do l5=1,l4
                                    t1=W(l5)-r(Nucd(k),l5)
                                    term=t1*dp1p

                                    if (l1.eq.l5) then
                                       term=term+pjpkpl/z2a
                                    endif

                                    if (l2.eq.l5) then
                                       term=term+pipkpl/z2a
                                    endif

                                    if (l3.eq.l5) then
                                       term=term+dspl/z2a
                                    endif

                                    f2=1.D0
                                    if (l4.eq.l5) then
                                       term=term+(dp-ro*d1p/ad(k,nk))/zc
                                       f2=sq3
                                    endif

                                    l12=Ll(l1)+l2
                                    l45=l4*(l4-1)/2+l5
                                    kk=k+l45-1

                                    cc=ccoef/(f1*f2)
                                    term=term*cc
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

                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   ! (dd|s)
   do i=ns+np+1,M,6
   do knan = 1, natomc(Nuc(ifunct))
      j=nnpd(jatc(knan,Nuc(i)))-6

      do kknan=1,nnd(jatc(knan,Nuc(i))),6
         j=j+6

         if (jfunct .le. ifunct) then
            fato=.true.
            fato2=.true.
            ddij=d(Nuc(i),Nuc(j))
            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij = a(ifunct,nci) + a(jfunct,ncj)
               z2=2.D0*Zij
               ti = a(ifunct,nci) / Zij
               tj = a(jfunct,ncj) / Zij
               alf = a(ifunct,nci) * tj
               rexp = alf * ddij
               if (rexp.lt.rmax) then
                  if (rexp.lt.rmaxs) then
                     if (fato) then
                        fato   = .false.
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
                        fato2  = .false.
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
                  Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
                  Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
                  Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)

                  sks=pi52*exp(-rexp)/Zij

                  do k=1,nsd
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+ &
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij

                        ti=Zij/t0
                        tj=ad(k,nk)/t0
                        W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                        W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                        W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                        roz=tj

                        ro=roz*Zij
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

                        do  l1=1,3
                           t1=Q(l1)-r(Nuc(i),l1)
                           t2=W(l1)-Q(l1)
                           ps=t1*sss+t2*ss1s
                           p1s=t1*ss1s+t2*ss2s
                           p2s=t1*ss2s+t2*ss3s
                           p3s=t1*ss3s+t2*ss4s
                           t6=(ps-roz*p1s)/z2
                           t7=(p1s-roz*p2s)/z2

                           do  l2=1,l1
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

                              f1=1.D0
                              if (l1.eq.l2) then
                                 ds=ds+t3
                                 d1s=d1s+t4
                                 d2s=d2s+t5
                                 f1=sq3
                              endif

                              t12=(ds-roz*d1s)/z2

                              lij=3
                              if (i.eq.j) then
                                 lij=l1
                              endif

                              do  l3=1,lij
                                 t1=Q(l3)-r(Nuc(j),l3)
                                 t2=W(l3)-Q(l3)
                                 pip=t1*ps+t2*p1s
                                 pi1p=t1*p1s+t2*p2s
                                 pjp=t1*pjs+t2*pj1s
                                 pj1p=t1*pj1s+t2*pj2s
                                 dp=t1*ds+t2*d1s
                                 d1p=t1*d1s+t2*d2s

                                 if (l1.eq.l3) then
                                    pip=pip+t3
                                    pi1p=pi1p+t4
                                    dp=dp+t8
                                    d1p=d1p+t9
                                 endif

                                 if (l2.eq.l3) then
                                    pjp=pjp+t3
                                    pj1p=pj1p+t4
                                    dp=dp+t6
                                    d1p=d1p+t7
                                 endif

                                 t10=(pjp-roz*pj1p)/z2
                                 t11=(pip-roz*pi1p)/z2

                                 lk=l3
                                 if (i.eq.j) then
                                    lk=min(l3,Ll(l1)-Ll(l3)+l2)
                                 endif

                                 do  l4=1,lk
                                    t1=Q(l4)-r(Nuc(j),l4)
                                    t2=W(l4)-Q(l4)
                                    term=t1*dp+t2*d1p

                                    if (l1.eq.l4) then
                                       term=term+t10
                                    endif

                                    if (l2.eq.l4) then
                                       term=term+t11
                                    endif

                                    f2=1.D0
                                    if (l3.eq.l4) then
                                       term=term+t12
                                       f2=sq3
                                    endif

                                    l12=Ll(l1)+l2
                                    l34=Ll(l3)+l4
                                    ii=ii+1
                                    cc=ccoef/(f1*f2)
                                    term=term*cc

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
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo

                  !----(dd|p)
                  do k=nsd+1,nsd+npd,3
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij
                        z2a=2.*t0

                        ti=Zij/t0
                        tj=ad(k,nk)/t0
                        W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                        W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                        W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                        roz=tj

                        ro=roz*Zij
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

                        do l1=1,3
                           t1=Q(l1)-r(Nuc(i),l1)
                           t2=W(l1)-Q(l1)
                           ps=t1*sss+t2*ss1s
                           p1s=t1*ss1s+t2*ss2s
                           p2s=t1*ss2s+t2*ss3s
                           p3s=t1*ss3s+t2*ss4s
                           p4s=t1*ss4s+t2*ss5s
                           t6=(p1s-roz*p2s)/z2
                           t7=(p2s-roz*p3s)/z2

                           do l2=1,l1
                              t1=Q(l2)-r(Nuc(i),l2)
                              t2=W(l2)-Q(l2)
                              pjs=t1*sss+t2*ss1s
                              pj1s=t1*ss1s+t2*ss2s
                              pj2s=t1*ss2s+t2*ss3s
                              pj3s=t1*ss3s+t2*ss4s
                              t8=(pj1s-roz*pj2s)/z2
                              t9=(pj2s-roz*pj3s)/z2

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

                              t18=(d1s-roz*d2s)/z2

                              lij=3
                              if (i.eq.j) then
                                 lij=l1
                              endif

                              do l3=1,lij
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

                                 if (l1.eq.l3) then
                                    d1pk=d1pk+t8
                                    d2pk=d2pk+t9
                                    pipk=pipk+ta
                                    pi1pk=pi1pk+t3
                                    pi2pk=pi2pk+t4
                                 endif

                                 if (l2.eq.l3) then
                                    d1pk=d1pk+t6
                                    d2pk=d2pk+t7
                                    pjpk=pjpk+ta
                                    pj1pk=pj1pk+t3
                                    pj2pk=pj2pk+t4
                                 endif

                                 lk=l3
                                 if (i.eq.j) then
                                    lk=min(l3,Ll(l1)-Ll(l3)+l2)
                                 endif
                                 t16=(pj1pk-roz*pj2pk)/z2
                                 t17=(pi1pk-roz*pi2pk)/z2

                                 do l4=1,lk
                                    t1=Q(l4)-r(Nuc(j),l4)
                                    t2=W(l4)-Q(l4)
                                    d1d=t1*d1pk+t2*d2pk

                                    pjdkl=t1*pj1pk+t2*pj2pk
                                    pidkl=t1*pi1pk+t2*pi2pk
                                    d1pl=t1*d1s+t2*d2s

                                    if (l1.eq.l4) then
                                       d1d=d1d+t16
                                       pidkl=pidkl+t10
                                       d1pl=d1pl+t8
                                    endif

                                    if (l2.eq.l4) then
                                       d1d=d1d+t17
                                       pjdkl=pjdkl+t10
                                       d1pl=d1pl+t6
                                    endif

                                    f2=1.D0
                                    if (l3.eq.l4) then
                                       d1d=d1d+t18
                                       pjdkl=pjdkl+t8
                                       pidkl=pidkl+t6
                                       f2=sq3
                                    endif

                                    t11=pjdkl/z2a
                                    t12=pidkl/z2a
                                    t13=d1pl/z2a
                                    t14=d1pk/z2a
                                    ii=ii+1
                                    do l5=1,3

                                       t1=W(l5)-r(Nucd(k),l5)
                                       term=t1*d1d

                                       if (l1.eq.l5) then
                                          term=term+t11
                                       endif

                                       if (l2.eq.l5) then
                                          term=term+t12
                                       endif
                                       if (l3.eq.l5) then
                                          term=term+t13
                                       endif

                                       if (l4.eq.l5) then
                                          term=term+t14
                                       endif
                                       kk=k+l5-1

                                       cc=ccoef/(f1*f2)
                                       term=term*cc
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
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo

                  !-----(dd|d)

                  do k=nsd+npd+1,Md,6
                     dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+&
                     (Q(3)-r(Nucd(k),3))**2

                     do nk=1,ncontd(k)
                        ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
                        t0=ad(k,nk)+Zij
                        z2a=2.D0*t0
                        zc=2.D0*ad(k,nk)

                        ti=Zij/t0
                        tj=ad(k,nk)/t0
                        W(1)=ti*Q(1)+tj*r(Nucd(k),1)
                        W(2)=ti*Q(2)+tj*r(Nucd(k),2)
                        W(3)=ti*Q(3)+tj*r(Nucd(k),3)

                        roz=tj
                        ro=roz*Zij
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

                        do l1=1,3
                           t1=Q(l1)-r(Nuc(i),l1)
                           t2=W(l1)-Q(l1)
                           ps=t1*sss+t2*ss1s
                           p1s=t1*ss1s+t2*ss2s
                           p2s=t1*ss2s+t2*ss3s
                           p3s=t1*ss3s+t2*ss4s
                           p4s=t1*ss4s+t2*ss5s
                           p5s=t1*ss5s+t2*ss6s

                           t7=(ps-roz*p1s)/z2
                           t8=(p1s-roz*p2s)/z2
                           t9=(p2s-roz*p3s)/z2
                           t10=(p3s-roz*p4s)/z2
                           do l2=1,l1
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

                              t11=(pjs-roz*pj1s)/z2
                              t12=(pj1s-roz*pj2s)/z2
                              t13=(pj2s-roz*pj3s)/z2
                              t14=(pj3s-roz*pj4s)/z2

                              f1=1.D0
                              if (l1.eq.l2) then
                                 ds=ds+t3
                                 d1s=d1s+t4
                                 d2s=d2s+t5
                                 d3s=d3s+t6
                                 d4s=d4s+t6b
                                 f1=sq3
                              endif

                              t16=(ds-roz*d1s)/z2
                              t17=(d1s-roz*d2s)/z2
                              t18=(d2s-roz*d3s)/z2
                              t22a=d2s/z2a

                              lij=3
                              if (i.eq.j) then
                                 lij=l1
                              endif

                              do l3=1,lij
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

                                 lk=l3
                                 if (i.eq.j) then
                                    lk=min(l3,Ll(l1)-Ll(l3)+l2)
                                 endif

                                 t20=pj2pk/z2a
                                 t21=pi2pk/z2a
                                 t22=d2pk/z2a

                                 t24=(pjpk-roz*pj1pk)/z2
                                 t25=(pj1pk-roz*pj2pk)/z2
                                 t26=(pj2pk-roz*pj3pk)/z2
                                 t27=(pipk-roz*pi1pk)/z2
                                 t28=(pi1pk-roz*pi2pk)/z2
                                 t29=(pi2pk-roz*pi3pk)/z2

                                 do l4=1,lk
                                    ii=ii+1
                                    t1=Q(l4)-r(Nuc(j),l4)
                                    t2=W(l4)-Q(l4)
                                    dd=t1*dpk+t2*d1pk
                                    d1d=t1*d1pk+t2*d2pk
                                    d2d=t1*d2pk+t2*d3pk

                                    pjdkl=t1*pj2pk+t2*pj3pk
                                    pidkl=t1*pi2pk+t2*pi3pk
                                    d2pl=t1*d2s+t2*d3s

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

                                    if (l2.eq.l4) then
                                       dd=dd+t27
                                       d1d=d1d+t28
                                       d2d=d2d+t29
                                       pjdkl=pjdkl+t15
                                       d2pl=d2pl+t9
                                       pj2pl=pj2pl+t5
                                    endif

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

                                    do l5=1,3

                                       t1=W(l5)-r(Nucd(k),l5)
                                       ddp=t1*d2d
                                       pjdklp=t1*pjdkl
                                       pidklp=t1*pidkl
                                       dijplp=t1*d2pl
                                       dijpkp=t1*d2pk

                                       if (l1.eq.l5) then
                                          ddp=ddp+t30
                                          pidklp=pidklp+t50
                                          dijplp=dijplp+t60
                                          dijpkp=dijpkp+t20
                                       endif

                                       if (l2.eq.l5) then
                                          ddp=ddp+t40
                                          pjdklp=pjdklp+t50
                                          dijplp=dijplp+t70
                                          dijpkp=dijpkp+t21
                                       endif

                                       if (l3.eq.l5) then
                                          ddp=ddp+t80
                                          pjdklp=pjdklp+t60
                                          pidklp=pidklp+t70
                                          dijpkp=dijpkp+t22a
                                       endif

                                       if (l4.eq.l5) then
                                          ddp=ddp+t22
                                          pjdklp=pjdklp+t20
                                          pidklp=pidklp+t21
                                          dijplp=dijplp+t22a
                                       endif

                                       t31=pjdklp/z2a
                                       t41=pidklp/z2a
                                       t51=dijplp/z2a
                                       t61=dijpkp/z2a

                                       do l6=1,l5

                                          t1=W(l6)-r(Nucd(k),l6)
                                          term=t1*ddp

                                          if (l1.eq.l6) then
                                             term=term+t31
                                          endif

                                          if (l2.eq.l6) then
                                             term=term+t41
                                          endif

                                          if (l3.eq.l6) then
                                             term=term+t51
                                          endif

                                          if (l4.eq.l6) then
                                             term=term+t61
                                          endif

                                          f3=1.D0
                                          if (l5.eq.l6) then
                                             term=term+t23
                                             f3=sq3
                                          endif
                                          cc=ccoef/(f1*f2*f3)
                                          term=term*cc
                                          l56=Ll(l5)+l6
                                          kk=k+l56-1
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
                                       enddo
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo
                  enddo
               endif
            enddo
            enddo
         endif
      enddo
   enddo
   enddo
   deallocate (Jx)
   return
end subroutine int3mem
end module subm_int3mem
