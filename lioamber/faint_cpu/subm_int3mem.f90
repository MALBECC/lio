module subm_int3mem
contains
subroutine int3mem()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Integrals subroutine - 2e integrals with 3 indexes : wavefunction and density!
! fitting functions, calculated using the Obara-Saika recursive method.        !
! Input :  Standard basis and density basis                                    !
! Output: Cool and Cools, terms used for Coulomb terms.                        !
!                                                                              !
! Precalculates integral terms and indexes for double (rmax) and single        !
! precision (rmaxs). If a term has a high value (r < rmaxs) it is calculated   !
! in single precision, if it has a lower valule (rmaxs < r < rmax) it is       !
! calculated in double precision. If its value is neglegible (rmax < r), it is !
! not calculated at all.                                                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use liotemp   , only: FUNCT
   use garcha_mod, only: cool, cools, kkind, kkinds, Nuc, Nucd, a, c,          &
                         d, r, ad, cd, natomc, nns, nnp, nnd, nnps, nnpp, nnpd,&
                         jatc, ncont, ncontd, nshell, nshelld, M, Md, rmax,    &
                         rmaxs, NORM, kknums, kknumd
   use constants_mod, only: pi52
   implicit none
   integer, dimension(:), allocatable :: Jx
   double precision  :: Q(3), W(3)
   double precision  :: ccoef, f1, f2, f3, rexp, sq3, term, uf, Z2, Z2a, Zc, Zij

   double precision  :: ta, ti, tj, t0, t1, t2, t3, t4, t5, t6, t6b, t7, t8,   &
                        t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t20,  &
                        t21, t22, t22a, t23, t24, t25, t26, t27, t28, t29, t30,&
                        t31, t40, t41, t50, t51, t60, t61, t70, t80
   double precision  :: sss, sks, spk, spjs, sspj, spjpk, sp2js, ss1s, ss2s,   &
                        ss3s, ss4s, ss5s, ss6s, ss1p, s1pk, s2pk, s3pk, spks,  &
                        spj, sdkl
   double precision  :: ps, pp, pp1p, p1s, p2s, p3s, p4s, p5s, pi1p, pi1pk,    &
                        pi2p, pi2pk, pi2pl, pi3pk, pijs, pij1s, pij2s, pispj,  &
                        pispk, pis1pk, pip, pipk, pipkpl, pidkl, pidklp, pjs,  &
                        pjs1pk, pj1s, pj1p, pj1pk, pj2s, pj2p, pj2pk, pj2pl,   &
                        pj3s, pj3pk, pj4s, pjp, pjpk, pjpkpl, pjdklp,pjdkl
   double precision  :: d0d, d1s, d1p, d1d, d1pk, d1pl, d2s, d2p, d2d, d2pl,   &
                        d2pk, d3s, d3pk, d4s, ds, ds1p, dspl, dp, dpc, dpk,    &
                        dp1p, ddp, dijplp, dijpkp

   integer           :: ns, nsd, nd, ndd, np, npd, kknan, knan, kknumsmax, lk, &
                        lij, l1, l2, l3, l4, l5, l6, l12, l23, l34, l45, l56,  &
                        ifunct, jfunct, kfunct, nci, ncj, nck, lcount,         &
                        cool_ind, Ll(3)
   logical           :: done_sp, done_dp

   allocate (Jx(M))
   ns  = nshell(0) ; np  = nshell(1) ; nd  = nshell(2)
   nsd = nshelld(0); npd = nshelld(1); ndd = nshelld(2)

   sq3=1.D0
   if (NORM) sq3 = dsqrt(3.D0)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1-1) / 2
   enddo
   do ifunct = 1, M
      Jx(ifunct) = (2*M-ifunct) * (ifunct-1) / 2
   enddo

   kknumd = 0
   kknums = 0

   ! Search for the dimensions of cool/cools.
   do ifunct = 1, ns
   do knan   = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         if (jfunct .le. ifunct) then
            done_sp = .true.
            done_dp = .true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               rexp = a(ifunct,nci) * a(jfunct,ncj) * &
                      d(Nuc(ifunct),Nuc(jfunct))    / &
                      (a(ifunct,nci) + a(jfunct,ncj))

               if (rexp .lt. rmax) then
               if (rexp .lt. rmaxs) then
                  if (done_sp) then
                     kknumd  = kknumd +1
                     done_sp = .false.
                  endif
               else
                  if (done_dp) then
                     kknums  = kknums +1
                     done_dp = .false.
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
         done_sp = .true.
         done_dp = .true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) /&
                   (a(ifunct,nci) + a(jfunct,ncj))

            if (rexp .lt. rmax) then
            if (rexp .lt. rmaxs) then
               if (done_sp) then
                  do l1 = 1, 3
                     kknumd = kknumd +1
                  enddo
                  done_sp = .false.
               endif
            else
               if (done_dp) then
                  do l1 = 1, 3
                     kknums = kknums +1
                  enddo
                  done_dp = .false.
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
      jfunct = nnpp(jatc(knan,Nuc(ifunct))) -3

      do kknan = 1, nnp(jatc(knan,Nuc(ifunct))), 3
         jfunct = jfunct +3
         if (jfunct .le. ifunct) then
            done_sp = .true.
            done_dp = .true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               rexp = a(ifunct,nci) * a(jfunct,ncj) * &
                      d(Nuc(ifunct),Nuc(jfunct)) /    &
                      (a(ifunct,nci) + a(jfunct,ncj))

               if (rexp .lt. rmax) then
               if (rexp .lt. rmaxs) then
                  if (done_sp) then
                  if (ifunct .eq. jfunct) then
                     do l1 = 1, 3
                     do l2 = 1, l1
                        kknumd = kknumd +1
                     enddo
                     enddo
                  else
                     do l1 = 1, 3
                     do l2 = 1, 3
                        kknumd = kknumd +1
                     enddo
                     enddo
                  endif
                     done_sp = .false.
                  endif
               else
                  if (done_dp) then
                  if (ifunct .eq. jfunct) then
                     do l1 = 1, 3
                     do l2 = 1, l1
                        kknums = kknums +1
                     enddo
                     enddo
                  else
                     do l1 = 1, 3
                     do l2 = 1, 3
                        kknums = kknums +1
                     enddo
                     enddo
                  endif
                     done_dp = .false.
                  endif
               endif
               endif
            enddo
            enddo
         endif
      enddo
   enddo
   enddo

   do ifunct = ns+np+1, M, 6
   do knan   = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1
      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         done_sp = .true.
         done_dp = .true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) /&
                   (a(ifunct,nci) + a(jfunct,ncj))

            if (rexp .lt. rmax) then
            if (rexp .lt. rmaxs) then
               if (done_sp) then
                  do l1 = 1, 6
                     kknumd = kknumd +1
                  enddo
                  done_sp = .false.
               endif
            else
               if (done_dp) then
                  do l1 = 1, 6
                     kknums = kknums +1
                  enddo
                  done_dp = .false.
               endif
            endif
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   do ifunct = ns+np+1, M, 6
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnpp(jatc(knan,Nuc(ifunct))) -3

      do kknan = 1, nnp(jatc(knan,Nuc(ifunct))), 3
         jfunct = jfunct +3
         done_sp   = .true.
         done_dp  = .true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) /&
                   (a(ifunct,nci) + a(jfunct,ncj))

            if (rexp .lt. rmax) then
            if (rexp .lt. rmaxs) then
               if (done_sp) then
                  do l1 = 1, 6
                  do l2 = 1, 3
                     kknumd = kknumd +1
                  enddo
                  enddo
                  done_sp = .false.
               endif
            else
               if (done_dp) then
                  do l1 = 1, 6
                  do l2 = 1, 3
                     kknums = kknums +1
                  enddo
                  enddo
                  done_dp = .false.
               endif
            endif
            endif
         enddo
         enddo
      enddo
   enddo
   enddo

   do ifunct = ns+np+1, M, 6
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnpd(jatc(knan,Nuc(ifunct))) -6

      do kknan=1, nnd(jatc(knan,Nuc(ifunct))), 6
         jfunct = jfunct +6
         if (jfunct .le. ifunct) then
            done_sp = .true.
            done_dp = .true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               rexp = a(ifunct,nci) * a(jfunct,ncj) * &
                      d(Nuc(ifunct),Nuc(jfunct))    / &
                      (a(ifunct,nci) + a(jfunct,ncj))

               if (rexp .lt. rmax) then
               if (rexp .lt. rmaxs) then
                  if (done_sp) then
                     done_sp = .false.
                     if (ifunct .eq. jfunct) then
                        do l1 = 1, 6
                        do l2 = 1, l1
                           kknumd = kknumd +1
                        enddo
                        enddo
                     else
                        do l1 = 1, 6
                        do l2 = 1, 6
                           kknumd = kknumd +1
                        enddo
                        enddo
                     endif
                  endif
               else
                  if (done_dp) then
                     done_dp = .false.
                     if (ifunct .eq. jfunct) then
                        do l1 = 1, 6
                        do l2 = 1, l1
                           kknums = kknums +1
                        enddo
                        enddo
                     else
                        do l1 = 1, 6
                        do l2 = 1, 6
                           kknums = kknums +1
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

   kknumsmax = kknums
   cool   = 0.0D0; cools  = 0.0
   kknumd = 0    ; kknums = 0

   ! Start of integrals.
   ! (ss|X) terms (X = s,p,d)
   do ifunct = 1, ns
   do knan   = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1

         if (jfunct .le. ifunct) then
            done_sp = .true.
            done_dp = .true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij  = a(ifunct,nci) + a(jfunct,ncj)
               ti   = a(ifunct,nci) / Zij
               tj   = a(jfunct,ncj) / Zij
               rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

               if (rexp .lt. rmax) then
                  if (rexp .lt. rmaxs) then
                     if (done_sp) then
                        kknumd        = kknumd +1
                        kkind(kknumd) = ifunct + Jx(jfunct)
                        done_sp       = .false.
                     endif
                  else
                     if (done_dp) then
                        kknums         = kknums +1
                        kkinds(kknums) = ifunct + Jx(jfunct)
                        done_dp        = .false.
                     endif
                  endif

                  Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
                  Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
                  Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
                  sks  = pi52 * exp(-rexp) / Zij

                  ! (ss|s)
                  do kfunct = 1, nsd
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        t0 = ad(kfunct,nck) + Zij
                        uf = ad(kfunct,nck) * Zij / t0 * dpc
                        term = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck) *&
                               sks * FUNCT(0,uf) / (ad(kfunct,nck) * dsqrt(t0))

                        if (rexp .lt. rmaxs) then
                           cool_ind = (kknumd -1) * Md + kfunct
                           cool(cool_ind) = cool(cool_ind) + term
                        else
                           cool_ind = (kknums -1) * Md + kfunct
                           cools(cool_ind) = cools(cool_ind) + real(term)
                        endif
                     enddo
                  enddo

                  ! (ss|p)
                  do kfunct = nsd+1, nsd+npd, 3
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        t0 = ad(kfunct,nck) + Zij
                        ti = Zij / t0
                        tj = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                        uf   = ad(kfunct,nck) * ti * dpc
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)   * c(ifunct,nci) &
                                  * c(jfunct,ncj) * cd(kfunct,nck)

                        do l1 = 1, 3
                           term = (W(l1) - r(Nucd(kfunct),l1)) * ss1s

                           if (rexp .lt. rmaxs) then
                              cool_ind = (kknumd -1) * Md + kfunct + l1 -1
                              cool(cool_ind) = cool(cool_ind) + term
                           else
                              cool_ind = (kknums -1) * Md + kfunct + l1 -1
                              cools(cool_ind) = cools(cool_ind) + real(term)
                           endif
                        enddo
                     enddo
                  enddo

                  ! (ss|d)
                  do kfunct = nsd+npd+1, Md, 6
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                        t0    = ad(kfunct,nck) + Zij
                        ti    = Zij / t0
                        tj    = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                        uf   = ti * ad(kfunct,nck) * dpc
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)
                        ss2s = t2 * FUNCT(2,uf)
                        ta   = (sss - ti * ss1s) / (2.D0 * ad(kfunct,nck))

                        do l1 = 1, 3
                           ss1p = (W(l1) - r(Nucd(kfunct),l1)) * ss2s

                           do l2 = 1, l1
                              term = (W(l2) - r(Nucd(kfunct),l2)) * ss1p
                              f1   = 1.D0
                              if (l1 .eq. l2) then
                                 term = term + ta
                                 f1   = sq3
                              endif
                              term = term * ccoef / f1

                              l12  = Ll(l1) + l2
                              if (rexp .lt. rmaxs) then
                                 cool_ind = (kknumd -1) * Md + kfunct + l12 -1
                                 cool(cool_ind) = cool(cool_ind) + term
                              else
                                 cool_ind = (kknums -1) * Md + kfunct + l12 -1
                                 cools(cool_ind) = cools(cool_ind) + real(term)
                              endif
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

   ! (ps|X)
   do ifunct = ns+1, ns+np, 3
   do knan   = 1   , natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         done_sp = .true.
         done_dp = .true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij  = a(ifunct,nci) + a(jfunct,ncj)
            ti   = a(ifunct,nci) / Zij
            tj   = a(jfunct,ncj) / Zij
            rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

            if (rexp .lt. rmax) then
               if (rexp .lt. rmaxs) then
                  if (done_sp) then
                     do l1 = 1, 3
                        kknumd = kknumd +1
                        kkind(kknumd) = ifunct + Jx(jfunct) + l1 -1
                     enddo
                     done_sp = .false.
                  endif
               else
                  if (done_dp) then
                     do l1 = 1, 3
                        kknums = kknums +1
                        kkinds(kknums) = ifunct + Jx(jfunct) + l1 -1
                     enddo
                     done_dp = .false.
                  endif
               endif

               Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
               Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
               Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
               sks  = pi52 * exp(-rexp) / Zij

               ! (ps|s)
               do kfunct = 1, nsd
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     uf   = ad(kfunct,nck) * ti * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)

                     do l1 = 1, 3
                        term = ccoef * ((Q(l1) - r(Nuc(ifunct),l1)) * sss + &
                                        (W(l1) - Q(l1)            ) * ss1s)

                        if (rexp .lt. rmaxs) then
                           cool_ind = (l1 + kknumd -4) * Md + kfunct
                           cool(cool_ind) = cool(cool_ind) + term
                        else
                           cool_ind = (l1 + kknums -4) * Md + kfunct
                           cools(cool_ind) = cools(cool_ind) + real(term)
                        endif
                     enddo
                  enddo
               enddo

               ! (ps|p)
               do kfunct = nsd+1, nsd+npd, 3
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     uf   = ad(kfunct,nck) * ti * dpc
                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ta   = ss1s / (2.0D0 * t0)

                     do l1 = 1, 3
                        p1s = (Q(l1) - r(Nuc(ifunct),l1)) * ss1s + &
                              (W(l1) - Q(l1)            ) * ss2s

                        do l2 = 1, 3
                           term = (W(l2) - r(Nucd(kfunct),l2)) * p1s
                           if (l1 .eq. l2) then
                              term = term + ta
                           endif
                           term = term * ccoef

                           if(rexp .lt. rmaxs) then
                              cool_ind = (l1 + kknumd-4) * Md + kfunct + l2 -1
                              cool(cool_ind) = cool(cool_ind) + term
                           else
                              cool_ind = (l1 + kknums-4) * Md + kfunct + l2 -1
                              cools(cool_ind) = cools(cool_ind) + real(term)
                           endif

                        enddo
                     enddo
                  enddo
               enddo

               ! (ps|d)
               do kfunct = nsd+npd+1, Md, 6
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     Zc    = 2.0D0 * ad(kfunct,nck)
                     Z2a   = 2.0D0 * t0
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)


                     t2   = sks/ (ad(kfunct,nck) * dsqrt(t0))
                     uf   = ti * ad(kfunct,nck) * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ss3s = t2 * FUNCT(3,uf)
                     t3   = ss2s / Z2a


                     do l1 = 1, 3
                        t1  = Q(l1) - r(Nuc(ifunct),l1)
                        t2  = W(l1) - Q(l1)
                        ps  = t1 * sss + t2 * ss1s
                        p1s = t1 * ss1s + t2 * ss2s
                        p2s = t1 * ss2s + t2 * ss3s
                        t5  = (ps - ti * p1s) / Zc

                        do l2 = 1, 3
                           t1=W(l2)-r(Nucd(kfunct),l2)
                           sspj = t1 * ss2s
                           pispj = t1 * p2s

                           t4=sspj / Z2a
                           if (l1 .eq. l2) then
                              pispj=pispj+t3
                           endif

                           do l3 = 1, l2
                              f1   = 1.D0
                              term = (W(l3) - r(Nucd(kfunct),l3)) * pispj
                              if (l1 .eq. l3) then
                                 term = term + t4
                              endif
                              if (l2 .eq. l3) then
                                 term = term + t5
                                 f1   = sq3
                              endif
                              term = term * ccoef / f1

                              l23 = Ll(l2) + l3
                              if(rexp .lt. rmaxs) then
                                 cool_ind = (l1 + kknumd -4)*Md + kfunct + l23-1
                                 cool(cool_ind) = cool(cool_ind) + term
                              else
                                 cool_ind = (l1 + kknums -4)*Md + kfunct + l23-1
                                 cools(cool_ind) = cools(cool_ind) + real(term)
                              endif
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

   ! (pp|X)
   do ifunct = ns+1, ns+np, 3
   do knan = 1, natomc(Nuc(ifunct))
      jfunct = nnpp(jatc(knan,Nuc(ifunct))) -3

      do kknan = 1, nnp(jatc(knan,Nuc(ifunct))), 3
         jfunct = jfunct +3

         if (jfunct .le. ifunct) then
            done_sp = .true.
            done_dp = .true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij  = a(ifunct,nci) + a(jfunct,ncj)
               Z2   = 2.D0 * Zij
               ti   = a(ifunct,nci) / Zij
               tj   = a(jfunct,ncj) / Zij
               rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

               if (rexp .lt. rmax) then
                  if (rexp .lt. rmaxs) then
                     if (done_sp) then
                        if (ifunct .eq. jfunct) then
                           do l1 = 1, 3
                           do l2 = 1, l1
                              kknumd = kknumd +1
                              kkind(kknumd) = ifunct + l1-1 + Jx(jfunct+l2-1)
                           enddo
                           enddo
                        else
                           do l1 =1, 3
                           do l2 =1, 3
                              kknumd = kknumd +1
                              kkind(kknumd) = ifunct + l1-1 + Jx(jfunct+l2-1)
                           enddo
                           enddo
                        endif
                        done_sp = .false.
                     endif
                  else
                     if (done_dp) then
                        if (ifunct .eq. jfunct) then
                           do l1 = 1, 3
                           do l2 = 1, l2
                              kknums = kknums +1
                              kkinds(kknums) = ifunct + l1-1 + Jx(jfunct+l2-1)
                           enddo
                           enddo
                        else
                           do l1 = 1, 3
                           do l2 = 1, 3
                              kknums = kknums +1
                              kkinds(kknums) = ifunct + l1-1 + Jx(jfunct+l2-1)
                           enddo
                           enddo
                        endif
                        done_dp = .false.
                     endif
                  endif

                  Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
                  Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
                  Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
                  sks  = pi52 * exp(-rexp) / Zij

                  ! (pp|s)
                  do kfunct = 1, nsd
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                        t0    = ad(kfunct,nck) + Zij
                        ti    = Zij / t0
                        tj    = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        uf   = ad(kfunct,nck) * ti * dpc
                        t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)
                        ss2s = t2 * FUNCT(2,uf)
                        ta   = (sss - tj * ss1s) / Z2

                        lcount = 0
                        do l1 = 1, 3
                           t1  = Q(l1) - r(Nuc(ifunct),l1)
                           t2  = W(l1) - Q(l1)
                           ps  = t1 * sss  + t2 * ss1s
                           p1s = t1 * ss1s + t2 * ss2s

                           lij = 3
                           if (ifunct .eq. jfunct) then
                              lij = l1
                           endif
                           do l2 = 1, lij
                              term = (Q(l2) - r(Nuc(jfunct),l2)) * ps + &
                                     (W(l2) - Q(l2)            ) * p1s
                              if (l1 .eq. l2) then
                                 term = term + ta
                              endif
                              term = term * ccoef

                              lcount = lcount +1
                              if (rexp .lt. rmaxs) then
                                 if (ifunct .eq. jfunct) then
                                    cool_ind = (lcount + kknumd -7) *Md + kfunct
                                 else
                                    cool_ind = (lcount + kknumd -10)*Md + kfunct
                                 endif
                                 cool(cool_ind) = cool(cool_ind) + term
                              else
                                 if (ifunct .eq. jfunct) then
                                    cool_ind = (lcount + kknums -7) *Md + kfunct
                                 else
                                    cool_ind = (lcount + kknums -10)*Md + kfunct
                                 endif
                                 cools(cool_ind) = cools(cool_ind) + real(term)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo

                  ! (pp|p)
                  do kfunct = nsd+1, nsd+npd, 3
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                        t0    = ad(kfunct,nck) + Zij
                        Z2a   = 2.0D0*t0
                        ti    = Zij / t0
                        tj    = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                        uf   = tj * Zij * dpc
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)
                        ss2s = t2 * FUNCT(2,uf)
                        ss3s = t2 * FUNCT(3,uf)
                        t3   = (ss1s - tj * ss2s) / Z2

                        lcount = 0
                        do l1 = 1, 3
                           t1  = Q(l1) - r(Nuc(ifunct),l1)
                           t2  = W(l1) - Q(l1)
                           p1s = t1 * ss1s + t2 * ss2s
                           p2s = t1 * ss2s + t2 * ss3s
                           t5  = p1s / Z2a

                           lij = 3
                           if (ifunct .eq. jfunct) then
                              lij = l1
                           endif

                           do l2 = 1, lij
                              t1  = Q(l2) - r(Nuc(jfunct),l2)
                              t2  = W(l2) - Q(l2)
                              spj = t1 * ss1s + t2 * ss2s
                              pp  = t1 * p1s  + t2 * p2s
                              t4  = spj / Z2a

                              if (l1 .eq. l2) then
                                 pp = pp + t3
                              endif

                              lcount = lcount +1
                              do l3 = 1, 3
                                 t1   = W(l3) - r(Nucd(kfunct),l3)
                                 term = t1 * pp

                                 if (l1 .eq. l3) then
                                    term = term + t4
                                 endif
                                 if (l2 .eq. l3) then
                                    term = term + t5
                                 endif
                                 term = term * ccoef

                                 if (rexp .lt. rmaxs) then
                                    if (ifunct .eq. jfunct) then
                                       cool_ind = (lcount + kknumd -7)  * Md + &
                                                  kfunct + l3 -1
                                    else
                                       cool_ind = (lcount + kknumd -10) * Md + &
                                                  kfunct + l3 -1
                                    endif
                                    cool(cool_ind) = cool(cool_ind) + term

                                 else
                                    if (ifunct .eq. jfunct) then
                                       cool_ind = (lcount + kknums -7)  * Md + &
                                                  kfunct + l3 -1
                                    else
                                       cool_ind = (lcount + kknums -10) * Md + &
                                                  kfunct + l3 -1
                                    endif
                                    cools(cool_ind)=cools(cool_ind) + real(term)
                                 endif
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo

                  !(pp|d)
                  do kfunct = nsd+npd+1, Md, 6
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                        t0    = ad(kfunct,nck) + Zij
                        Z2a   = 2.0D0 * t0
                        ti    = Zij / t0
                        tj    = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        Zc = 2.D0 * ad(kfunct,nck)
                        t1 = ad(kfunct,nck) * dsqrt(t0)
                        t2 = sks / t1
                        uf = tj * Zij * dpc
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)
                        ss2s = t2 * FUNCT(2,uf)
                        ss3s = t2 * FUNCT(3,uf)
                        ss4s = t2 * FUNCT(4,uf)
                        t3 = (sss  - tj * ss1s) / Z2
                        t4 = (ss1s - tj * ss2s) / Z2
                        t5 = (ss2s - tj * ss3s) / Z2
                        t6 = ss2s / Z2a

                        lcount = 0
                        do l1 = 1, 3
                           t1  = Q(l1) - r(Nuc(ifunct),l1)
                           t2  = W(l1) - Q(l1)
                           ps  = t1 * sss  + t2 * ss1s
                           p1s = t1 * ss1s + t2 * ss2s
                           p2s = t1 * ss2s + t2 * ss3s
                           p3s = t1 * ss3s + t2 * ss4s
                           t8  = p2s / Z2a

                           lij = 3
                           if (ifunct .eq. jfunct) then
                              lij = l1
                           endif
                           do l2 = 1, lij
                              t1 = Q(l2) - r(Nuc(jfunct),l2)
                              t2 = W(l2) - Q(l2)
                              pijs  = t1 * ps   + t2 * p1s
                              pij1s = t1 * p1s  + t2 * p2s
                              pij2s = t1 * p2s  + t2 * p3s
                              spjs  = t1 * ss1s + t2 * ss2s
                              sp2js = t1 * ss2s + t2 * ss3s
                              t7    = sp2js / Z2a

                              if (l1 .eq. l2) then
                                 pijs  = pijs  + t3
                                 pij1s = pij1s + t4
                                 pij2s = pij2s + t5
                              endif
                              t11 = (pijs - ti * pij1s) / Zc

                              lcount = lcount +1
                              do l3 = 1, 3
                                 t1    = W(l3) - r(Nucd(kfunct),l3)
                                 pp1p  = t1 * pij2s
                                 spjpk = t1 * sp2js
                                 pispk = t1 * p2s

                                 if (l1 .eq. l3) then
                                    pp1p  = pp1p  + t7
                                    pispk = pispk + t6
                                 endif
                                 if (l2 .eq. l3) then
                                    pp1p  = pp1p  + t8
                                    spjpk = spjpk + t6
                                 endif
                                 t9  = spjpk / Z2a
                                 t10 = pispk / Z2a

                                 do l4 = 1, l3
                                    term = (W(l4) - r(Nucd(kfunct),l4)) * pp1p
                                    f1   = 1.D0

                                    if (l1 .eq. l4) then
                                       term = term + t9
                                    endif
                                    if (l2 .eq. l4) then
                                       term = term + t10
                                    endif
                                    if (l3 .eq. l4) then
                                       term = term + t11
                                       f1   = sq3
                                    endif
                                    term = term * ccoef / f1

                                    l34 = Ll(l3) + l4
                                    if (rexp .lt. rmaxs) then
                                       if (ifunct .eq. jfunct) then
                                          cool_ind = (lcount + kknumd -7 )*Md +&
                                                     kfunct + l34 -1
                                       else
                                          cool_ind = (lcount + kknumd -10)*Md +&
                                                     kfunct + l34 -1
                                       endif
                                       cool(cool_ind) = cool(cool_ind) + term
                                    else
                                       if (ifunct .eq. jfunct) then
                                          cool_ind = (lcount + kknums -7 )*Md +&
                                                     kfunct + l34 -1
                                       else
                                          cool_ind = (lcount + kknums -10)*Md +&
                                                     kfunct + l34 -1
                                       endif
                                       cools(cool_ind) = cools(cool_ind) + &
                                                         real(term)
                                    endif
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

   ! (ds|X)
   do ifunct = ns+np+1, M, 6
   do knan   = 1, natomc(Nuc(ifunct))
      jfunct = nnps(jatc(knan, Nuc(ifunct))) -1

      do kknan = 1, nns(jatc(knan, Nuc(ifunct)))
         jfunct = jfunct +1
         done_sp = .true.
         done_dp = .true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij  = a(ifunct,nci) + a(jfunct,ncj)
            Z2   = 2.D0*Zij
            ti   = a(ifunct,nci) / Zij
            tj   = a(jfunct,ncj) / Zij
            rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

            if (rexp .lt. rmax) then
               if (rexp .lt. rmaxs) then
                  if (done_sp) then
                     do l1 = 1, 6
                        kknumd = kknumd +1
                        kkind(kknumd) = ifunct + l1 -1 + Jx(jfunct)
                     enddo
                     done_sp = .false.
                  endif
               else
                  if (done_dp) then
                     do l1 = 1, 6
                        kknums = kknums +1
                        kkinds(kknums) = ifunct + l1 -1 +Jx(jfunct)
                     enddo
                     done_dp = .false.
                  endif
               endif

               Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
               Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
               Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
               sks  = pi52 * exp(-rexp) / Zij

               ! (ds|s)
               do kfunct = 1, nsd
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     uf   = tj * Zij * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ta   = (sss - tj * ss1s) / Z2

                     do l1 = 1, 3
                        t1  = Q(l1) - r(Nuc(ifunct),l1)
                        t2  = W(l1) - Q(l1)
                        ps  = t1 * sss + t2 * ss1s
                        p1s = t1 * ss1s + t2 * ss2s

                        do l2 = 1, l1
                           term = (Q(l2) - r(Nuc(ifunct),l2)) * ps + &
                                  (W(l2) - Q(l2)            ) * p1s
                           f1   = 1.D0
                           if (l1 .eq. l2) then
                              term = term + ta
                              f1   = sq3
                           endif
                           term = term * ccoef / f1

                           l12 = Ll(l1) + l2
                           if (rexp .lt. rmaxs) then
                              cool_ind = (l12 + kknumd -7) * Md + kfunct
                              cool(cool_ind) = cool(cool_ind) + term
                           else
                              cool_ind = (l12 + kknums -7) * Md + kfunct
                              cools(cool_ind) = cools(cool_ind) + real(term)
                           endif
                        enddo
                     enddo
                  enddo
               enddo

               !---(ds|p)
               do kfunct = nsd+1, nsd+npd, 3
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     Z2a   = 2.0D0 * t0
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     uf   = tj * Zij * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ss3s = t2 * FUNCT(3,uf)
                     t3   = (ss1s - tj * ss2s) / Z2

                     do l1 = 1, 3
                        t1  = Q(l1) - r(Nuc(ifunct),l1)
                        t2  = W(l1) - Q(l1)
                        p1s = t1 * ss1s + t2 * ss2s
                        p2s = t1 * ss2s + t2 * ss3s
                        t5  = p1s / Z2a

                        do l2 = 1, l1
                           t1   = Q(l2) - r(Nuc(ifunct),l2)
                           t2   = W(l2) - Q(l2)
                           pj1s = t1 * ss1s + t2 * ss2s
                           t4   = pj1s / Z2a
                           ds   = t1 * p1s + t2 * p2s

                           f1 = 1.D0
                           if (l1 .eq. l2) then
                              ds = ds + t3
                              f1 = sq3
                           endif

                           do l3 = 1, 3
                              term = (W(l3) - r(Nucd(kfunct),l3)) * ds
                              if (l1 .eq. l3) then
                                 term = term + t4
                              endif
                              if (l2 .eq. l3) then
                                 term = term + t5
                              endif
                              term = term * ccoef / f1

                              l12 = Ll(l1) + l2
                              if (rexp .lt. rmaxs) then
                                 cool_ind = (l12 + kknumd -7)*Md + kfunct + l3-1
                                 cool(cool_ind) = cool(cool_ind) + term
                              else
                                 cool_ind = (l12 + kknums -7)*Md + kfunct + l3-1
                                 cools(cool_ind) = cools(cool_ind) + real(term)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

               !------(ds|d)
               do kfunct = nsd+npd+1, Md, 6
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     Z2a   = 2.D0 * t0
                     Zc    = 2.D0 * ad(kfunct,nck)
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     uf   = tj * Zij * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ss3s = t2 * FUNCT(3,uf)
                     ss4s = t2 * FUNCT(4,uf)
                     t3   = (sss  - tj * ss1s) / Z2
                     t4   = (ss1s - tj * ss2s) / Z2
                     t5   = (ss2s - tj * ss3s) / Z2
                     t6   = ss2s / Z2a

                     do l1 = 1, 3
                        t1  = Q(l1) - r(Nuc(ifunct),l1)
                        t2  = W(l1) - Q(l1)
                        ps  = t1 * sss + t2 * ss1s
                        p1s = t1 * ss1s + t2 * ss2s
                        p2s = t1 * ss2s + t2 * ss3s
                        p3s = t1 * ss3s + t2 * ss4s
                        t7  = p2s / Z2a

                        do l2 = 1, l1
                           t1   = Q(l2) - r(Nuc(ifunct),l2)
                           t2   = W(l2) - Q(l2)
                           pj1s = t1 * ss1s + t2 * ss2s
                           pj2s = t1 * ss2s + t2 * ss3s
                           ds   = t1 * ps + t2 * p1s
                           d1s  = t1 * p1s + t2 * p2s
                           d2s  = t1 * p2s + t2 * p3s
                           t8   = pj2s / Z2a

                           f1   = 1.0D0
                           if (l1 .eq. l2) then
                              ds  = ds  + t3
                              d1s = d1s + t4
                              d2s = d2s + t5
                              f1  = sq3
                           endif

                           t11 = (ds - ti * d1s) / Zc
                           do l3 = 1, 3
                              t1 = W(l3) - r(Nucd(kfunct),l3)
                              ds1p   = t1 * d2s
                              pis1pk = t1 * p2s
                              pjs1pk = t1 * pj2s

                              if (l1 .eq. l3) then
                                 ds1p   = ds1p   + t8
                                 pis1pk = pis1pk + t6
                              endif
                              if (l2 .eq. l3) then
                                 ds1p   = ds1p   + t7
                                 pjs1pk = pjs1pk + t6
                              endif
                              t9  = pjs1pk / Z2a
                              t10 = pis1pk / Z2a

                              do l4 = 1, l3
                                 term = (W(l4) - r(Nucd(kfunct),l4)) * ds1p
                                 f2   = 1.D0

                                 if (l1 .eq. l4) then
                                    term = term + t9
                                 endif
                                 if (l2 .eq. l4) then
                                    term = term + t10
                                 endif
                                 if (l3 .eq. l4) then
                                    term = term + t11
                                    f2   = sq3
                                 endif
                                 term = term * ccoef / (f1 * f2)

                                 l12 = Ll(l1) + l2
                                 l34 = Ll(l3) + l4
                                 if (rexp .lt. rmaxs) then
                                    cool_ind = (l12 + kknumd -7) * Md + kfunct &
                                             + l34 -1
                                    cool(cool_ind) = cool(cool_ind) + term
                                 else
                                    cool_ind = (l12 + kknums -7) * Md + kfunct &
                                             + l34 -1
                                    cools(cool_ind) = cools(cool_ind) + &
                                                      real(term)
                                 endif
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

   ! (dp|X)
   do ifunct = ns+np+1, M, 6
   do knan   = 1, natomc(Nuc(ifunct))
      jfunct = nnpp(jatc(knan,Nuc(ifunct))) -3

      do kknan = 1, nnp(jatc(knan,Nuc(ifunct))), 3
         jfunct = jfunct +3
         done_sp = .true.
         done_dp = .true.

         do nci = 1, ncont(ifunct)
         do ncj = 1, ncont(jfunct)
            Zij  = a(ifunct,nci) + a(jfunct,ncj)
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            tj   = a(jfunct,ncj) / Zij
            rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

            if (rexp .lt. rmax) then
               if (rexp .lt. rmaxs) then
                  if (done_sp) then
                     do l1 = 1, 6
                     do l2 = 1, 3
                        kknumd = kknumd +1
                        kkind(kknumd) = ifunct + l1 -1 + Jx(jfunct + l2 -1)
                     enddo
                     enddo
                     done_sp = .false.
                  endif
               else
                  if (done_dp) then
                     do l1 = 1, 6
                     do l2 = 1, 3
                        kknums = kknums +1
                        kkinds(kknums) = ifunct + l1 -1 + Jx(jfunct + l2 -1)
                     enddo
                     enddo
                     done_dp = .false.
                  endif
               endif

               Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
               Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
               Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
               sks  = pi52 * exp(-rexp) / Zij

               ! (dp|s)
               do kfunct = 1, nsd
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     uf   = tj * Zij * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ss3s = t2 * FUNCT(3,uf)
                     t3   = (sss  - tj * ss1s) / Z2
                     t4   = (ss1s - tj * ss2s) / Z2

                     lcount = 0
                     do l1 = 1, 3
                        t1  = Q(l1) - r(Nuc(ifunct),l1)
                        t2  = W(l1) - Q(l1)
                        ps  = t1 * sss  + t2 * ss1s
                        p1s = t1 * ss1s + t2 * ss2s
                        p2s = t1 * ss2s + t2 * ss3s
                        t5  = (ps - tj * p1s) / Z2

                        do l2 = 1, l1
                           t1   = Q(l2) - r(Nuc(ifunct),l2)
                           t2   = W(l2) - Q(l2)
                           pjs  = t1 * sss  + t2 * ss1s
                           pj1s = t1 * ss1s + t2 * ss2s
                           ds   = t1 * ps   + t2 * p1s
                           d1s  = t1 * p1s  + t2 * p2s
                           t6   = (pjs - tj * pj1s) / Z2

                           f1 = 1.D0
                           if (l1 .eq. l2) then
                              f1  = sq3
                              ds  = ds  + t3
                              d1s = d1s + t4
                           endif

                           do l3 = 1, 3
                              term = (Q(l3) - r(Nuc(jfunct),l3)) * ds + &
                                     (W(l3) - Q(l3)            ) * d1s
                              if (l1 .eq. l3) then
                                 term = term + t6
                              endif
                              if (l2 .eq. l3) then
                                 term = term + t5
                              endif
                              term = term *ccoef/f1

                              l12    = Ll(l1) + l2
                              lcount = lcount +1
                              if (rexp .lt. rmaxs) then
                                 cool_ind = (lcount + kknumd -19) * Md + kfunct
                                 cool(cool_ind) = cool(cool_ind) + term
                              else
                                 cool_ind = (lcount + kknums -19) * Md + kfunct
                                 cools(cool_ind) = cools(cool_ind) + real(term)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

               ! (dp|p)
               do kfunct = nsd+1, nsd+npd, 3
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     Z2a   = 2.0D0 * t0
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     uf   = tj * Zij * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ss3s = t2 * FUNCT(3,uf)
                     ss4s = t2 * FUNCT(4,uf)
                     t3   = (ss1s - tj * ss2s) / Z2
                     t4   = (ss2s - tj * ss3s) / Z2

                     lcount = 0
                     do l1 = 1, 3
                        t1  = Q(l1) - r(Nuc(ifunct),l1)
                        t2  = W(l1) - Q(l1)
                        ps  = t1 * sss  + t2 * ss1s
                        p1s = t1 * ss1s + t2 * ss2s
                        p2s = t1 * ss2s + t2 * ss3s
                        p3s = t1 * ss3s + t2 * ss4s
                        t5  = (p1s - tj * p2s) / Z2

                        do l2 = 1, l1
                           t1   = Q(l2) - r(Nuc(ifunct),l2)
                           t2   = W(l2) - Q(l2)
                           d1s  = t1 * p1s  + t2 * p2s
                           d2s  = t1 * p2s  + t2 * p3s
                           pj1s = t1 * ss1s + t2 * ss2s
                           pj2s = t1 * ss2s + t2 * ss3s
                           t6   = (pj1s - tj * pj2s) / Z2

                           f1 = 1.D0
                           if (l1 .eq. l2) then
                              d1s = d1s + t3
                              d2s = d2s + t4
                              f1  = sq3
                           endif
                           t9 = d1s / Z2a

                           do l3 = 1, 3
                              t1   = Q(l3) - r(Nuc(jfunct),l3)
                              t2   = W(l3) - Q(l3)
                              d1p  = t1 * d1s  + t2 * d2s
                              pi1p = t1 * p1s  + t2 * p2s
                              pj1p = t1 * pj1s + t2 * pj2s

                              if (l1 .eq. l3) then
                                 d1p  = d1p  + t6
                                 pi1p = pi1p + t3
                              endif
                              if (l2 .eq. l3) then
                                 d1p  = d1p  + t5
                                 pj1p = pj1p + t3
                              endif
                              t7 = pi1p / Z2a
                              t8 = pj1p / Z2a

                              lcount = lcount +1
                              do l4 = 1, 3
                                 term = (W(l4) - r(Nucd(kfunct),l4)) * d1p
                                 if (l1 .eq. l4) then
                                    term = term +t8
                                 endif
                                 if (l2 .eq. l4) then
                                    term = term +t7
                                 endif
                                 if (l3 .eq. l4) then
                                    term = term +t9
                                 endif
                                 term = term * ccoef / f1

                                 l12=Ll(l1)+l2
                                 if (rexp .lt. rmaxs) then
                                    cool_ind = (lcount + kknumd -19) * Md + &
                                               kfunct + l4 -1
                                    cool(cool_ind) = cool(cool_ind) + term
                                 else
                                    cool_ind = (lcount + kknums -19) * Md + &
                                               kfunct + l4 -1
                                    cools(cool_ind)=cools(cool_ind) + real(term)
                                 endif
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

               ! (dp|d)
               do kfunct = nsd+npd+1, Md, 6
                  dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                      + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                      + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                  do nck = 1, ncontd(kfunct)
                     ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                     t0    = ad(kfunct,nck) + Zij
                     Z2a   = 2.D0 * t0
                     Zc    = 2.D0 * ad(kfunct,nck)
                     ti    = Zij / t0
                     tj    = ad(kfunct,nck) / t0

                     W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                     W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                     W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                     t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                     uf   = tj * Zij * dpc
                     sss  = t2 * FUNCT(0,uf)
                     ss1s = t2 * FUNCT(1,uf)
                     ss2s = t2 * FUNCT(2,uf)
                     ss3s = t2 * FUNCT(3,uf)
                     ss4s = t2 * FUNCT(4,uf)
                     ss5s = t2 * FUNCT(5,uf)

                     lcount = 0
                     do l1 = 1, 3
                        t1  = Q(l1) - r(Nuc(ifunct),l1)
                        t2  = W(l1) - Q(l1)
                        ps  = t1 * sss  + t2 * ss1s
                        p1s = t1 * ss1s + t2 * ss2s
                        p2s = t1 * ss2s + t2 * ss3s
                        p3s = t1 * ss3s + t2 * ss4s
                        p4s = t1 * ss4s + t2 * ss5s

                        do l2 = 1, l1
                           t1   = Q(l2) - r(Nuc(ifunct),l2)
                           t2   = W(l2) - Q(l2)
                           ds   = t1 * ps   + t2 * p1s
                           d1s  = t1 * p1s  + t2 * p2s
                           d2s  = t1 * p2s  + t2 * p3s
                           d3s  = t1 * p3s  + t2 * p4s
                           pjs  = t1 * sss  + t2 * ss1s
                           pj1s = t1 * ss1s + t2 * ss2s
                           pj2s = t1 * ss2s + t2 * ss3s
                           pj3s = t1 * ss3s + t2 * ss4s

                           f1 = 1.0D0
                           if (l1 .eq. l2) then
                              ds  = ds  + (sss  - tj * ss1s) / Z2
                              d1s = d1s + (ss1s - tj * ss2s) / Z2
                              d2s = d2s + (ss2s - tj * ss3s) / Z2
                              d3s = d3s + (ss3s - tj * ss4s) / Z2
                              f1  = sq3
                           endif

                           do l3 = 1, 3
                              t1 = Q(l3) - r(Nuc(jfunct),l3)
                              t2 = W(l3) - Q(l3)
                              spks = t1 * ss2s + t2 * ss3s
                              dp   = t1 * ds   + t2 * d1s
                              d1p  = t1 * d1s  + t2 * d2s
                              d2p  = t1 * d2s  + t2 * d3s
                              pi1p = t1 * p1s  + t2 * p2s
                              pi2p = t1 * p2s  + t2 * p3s
                              pj1p = t1 * pj1s + t2 * pj2s
                              pj2p = t1 * pj2s + t2 * pj3s

                              if (l1 .eq. l3) then
                                 dp   = dp   + (pjs  - tj * pj1s) / Z2
                                 d1p  = d1p  + (pj1s - tj * pj2s) / Z2
                                 d2p  = d2p  + (pj2s - tj * pj3s) / Z2
                                 pi1p = pi1p + (ss1s - tj * ss2s) / Z2
                                 pi2p = pi2p + (ss2s - tj * ss3s) / Z2
                              endif
                              if (l2 .eq. l3) then
                                 dp   = dp   + (ps   - tj * p1s ) / Z2
                                 d1p  = d1p  + (p1s  - tj * p2s ) / Z2
                                 d2p  = d2p  + (p2s  - tj * p3s ) / Z2
                                 pj1p = pj1p + (ss1s - tj * ss2s) / Z2
                                 pj2p = pj2p + (ss2s - tj * ss3s) / Z2
                              endif

                              lcount = lcount +1
                              do l4 = 1, 3
                                 t1     = W(l4) - r(Nucd(kfunct),l4)
                                 dp1p   = t1 * d2p
                                 pjpkpl = t1 * pj2p
                                 pipkpl = t1 * pi2p
                                 dspl   = t1 * d2s

                                 if (l1 .eq. l4) then
                                    dp1p   = dp1p   + pj2p / Z2a
                                    pipkpl = pipkpl + spks / Z2a
                                    dspl   = dspl   + pj2s / Z2a
                                 endif
                                 if (l2 .eq. l4) then
                                    dp1p   = dp1p   + pi2p / Z2a
                                    pjpkpl = pjpkpl + spks / Z2a
                                    dspl   = dspl   + p2s  / Z2a
                                 endif
                                 if (l3 .eq. l4) then
                                    dp1p   = dp1p   + d2s  / Z2a
                                    pipkpl = pipkpl + p2s  / Z2a
                                    pjpkpl = pjpkpl + pj2s / Z2a
                                 endif

                                 do l5 = 1, l4
                                    term = (W(l5) - r(Nucd(kfunct),l5)) * dp1p
                                    f2   = 1.0D0
                                    if (l1 .eq. l5) then
                                       term = term + pjpkpl / Z2a
                                    endif
                                    if (l2 .eq. l5) then
                                       term = term + pipkpl / Z2a
                                    endif
                                    if (l3 .eq. l5) then
                                       term = term + dspl  / Z2a
                                    endif
                                    if (l4 .eq. l5) then
                                       term = term + (dp - Zij * d1p / t0) / Zc
                                       f2   = sq3
                                    endif
                                    term = term * ccoef / (f1 * f2)

                                    l45 = Ll(l4) + l5
                                    if (rexp .lt. rmaxs) then
                                       cool_ind = (lcount + kknumd -19) * Md + &
                                                  kfunct + l45 -1
                                       cool(cool_ind) = cool(cool_ind) + term
                                    else
                                       cool_ind = (lcount + kknums -19) * Md + &
                                                   kfunct + l45 -1
                                       cools(cool_ind) = cools(cool_ind) + &
                                                         real(term)
                                    endif
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
   do ifunct = ns+np+1, M, 6
   do knan   = 1, natomc(Nuc(ifunct))
      jfunct = nnpd(jatc(knan,Nuc(ifunct))) -6

      do kknan=1, nnd(jatc(knan,Nuc(ifunct))), 6
         jfunct = jfunct +6

         if (jfunct .le. ifunct) then
            done_sp = .true.
            done_dp = .true.

            do nci = 1, ncont(ifunct)
            do ncj = 1, ncont(jfunct)
               Zij  = a(ifunct,nci) + a(jfunct,ncj)
               Z2   = 2.D0 * Zij
               ti   = a(ifunct,nci) / Zij
               tj   = a(jfunct,ncj) / Zij
               rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

               if (rexp .lt. rmax) then
                  if (rexp .lt. rmaxs) then
                     if (done_sp) then
                        done_sp = .false.
                        if (ifunct .eq. jfunct) then
                           do l1 = 1, 6
                           do l2 = 1, l1
                              kknumd = kknumd +1
                              kkind(kknumd) = ifunct + l1-1 + Jx(jfunct + l2-1)
                           enddo
                           enddo
                        else
                           do l1 = 1, 6
                           do l2 = 1, 6
                              kknumd = kknumd +1
                              kkind(kknumd) = ifunct + l1-1 + Jx(jfunct + l2-1)
                           enddo
                           enddo
                        endif
                     endif
                  else
                     if (done_dp) then
                        done_dp = .false.
                        if (ifunct .eq. jfunct) then
                           do l1 = 1, 6
                           do l2 = 1, l1
                              kknums = kknums +1
                              kkinds(kknums) = ifunct + l1-1 + Jx(jfunct + l2-1)
                           enddo
                           enddo
                        else
                           do l1 = 1, 6
                           do l2 = 1, 6
                              kknums = kknums +1
                              kkinds(kknums) = ifunct + l1-1 + Jx(jfunct + l2-1)
                           enddo
                           enddo
                        endif
                     endif
                  endif

                  Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
                  Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
                  Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
                  sks  = pi52 * exp(-rexp) / Zij

                  do kfunct = 1, nsd
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                        t0    = ad(kfunct,nck) + Zij
                        ti    = Zij / t0
                        tj    = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                        uf   = tj * Zij * dpc
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)
                        ss2s = t2 * FUNCT(2,uf)
                        ss3s = t2 * FUNCT(3,uf)
                        ss4s = t2 * FUNCT(4,uf)
                        t3   = (sss  - tj * ss1s) / Z2
                        t4   = (ss1s - tj * ss2s) / Z2
                        t5   = (ss2s - tj * ss3s) / Z2

                        lcount = 0
                        do  l1=1,3
                           t1  = Q(l1) - r(Nuc(ifunct),l1)
                           t2  = W(l1) - Q(l1)
                           ps  = t1 * sss + t2 * ss1s
                           p1s = t1 * ss1s + t2 * ss2s
                           p2s = t1 * ss2s + t2 * ss3s
                           p3s = t1 * ss3s + t2 * ss4s
                           t6  = (ps  - tj * p1s) / Z2
                           t7  = (p1s - tj * p2s) / Z2

                           do  l2=1,l1
                              t1   = Q(l2) - r(Nuc(ifunct),l2)
                              t2   = W(l2) - Q(l2)
                              pjs  = t1 * sss  + t2 * ss1s
                              pj1s = t1 * ss1s + t2 * ss2s
                              pj2s = t1 * ss2s + t2 * ss3s
                              ds   = t1 * ps   + t2 * p1s
                              d1s  = t1 * p1s  + t2 * p2s
                              d2s  = t1 * p2s  + t2 * p3s
                              t8   = (pjs  - tj * pj1s) / Z2
                              t9   = (pj1s - tj * pj2s) / Z2

                              f1 = 1.D0
                              if (l1 .eq. l2) then
                                 ds  = ds  + t3
                                 d1s = d1s + t4
                                 d2s = d2s + t5
                                 f1  = sq3
                              endif
                              t12 = (ds - tj * d1s) / Z2

                              lij = 3
                              if (ifunct .eq. jfunct) then
                                 lij = l1
                              endif
                              do  l3 = 1, lij
                                 t1   = Q(l3) - r(Nuc(jfunct),l3)
                                 t2   = W(l3) - Q(l3)
                                 pip  = t1 * ps   + t2 * p1s
                                 pi1p = t1 * p1s  + t2 * p2s
                                 pjp  = t1 * pjs  + t2 * pj1s
                                 pj1p = t1 * pj1s + t2 * pj2s
                                 dp   = t1 * ds   + t2 * d1s
                                 d1p  = t1 * d1s  + t2 * d2s

                                 if (l1 .eq. l3) then
                                    pip  = pip  + t3
                                    pi1p = pi1p + t4
                                    dp   = dp   + t8
                                    d1p  = d1p  + t9
                                 endif
                                 if (l2 .eq. l3) then
                                    pjp  = pjp  + t3
                                    pj1p = pj1p + t4
                                    dp   = dp   + t6
                                    d1p  = d1p  + t7
                                 endif
                                 t10 = (pjp - tj * pj1p) / Z2
                                 t11 = (pip - tj * pi1p) / Z2

                                 lk = l3
                                 if (ifunct .eq. jfunct) then
                                    lk = min(l3, Ll(l1)-Ll(l3)+l2)
                                 endif
                                 do l4 = 1, lk
                                    term = (Q(l4) - r(Nuc(jfunct),l4)) * dp &
                                         + (W(l4) - Q(l4)            ) * d1p
                                    f2   = 1.D0

                                    if (l1 .eq. l4) then
                                       term = term +t10
                                    endif
                                    if (l2 .eq. l4) then
                                       term = term +t11
                                    endif
                                    if (l3 .eq. l4) then
                                       term = term +t12
                                       f2   = sq3
                                    endif
                                    term = term * ccoef / (f1 * f2)

                                    lcount = lcount +1
                                    if (rexp .lt. rmaxs) then
                                       if (jfunct .eq. ifunct) then
                                          cool_ind = (lcount + kknumd -22) *Md &
                                                   + kfunct
                                       else
                                          cool_ind = (lcount + kknumd -37) *Md &
                                                   + kfunct
                                       endif
                                       cool(cool_ind) = cool(cool_ind) + term
                                    else
                                       if (jfunct .eq. ifunct) then
                                          cool_ind = (lcount + kknums -22) *Md &
                                                   + kfunct
                                       else
                                          cool_ind = (lcount + kknums -37) *Md &
                                                   + kfunct
                                       endif
                                       cools(cool_ind) = cools(cool_ind) &
                                                       + real(term)
                                    endif
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo

                  !----(dd|p)
                  do kfunct = nsd+1, nsd+npd, 3
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                        t0    = ad(kfunct,nck) + Zij
                        Z2a   = 2.0D0 * t0
                        ti    = Zij / t0
                        tj    = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                        uf   = tj * Zij * dpc
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)
                        ss2s = t2 * FUNCT(2,uf)
                        ss3s = t2 * FUNCT(3,uf)
                        ss4s = t2 * FUNCT(4,uf)
                        ss5s = t2 * FUNCT(5,uf)
                        ta   = (sss  - tj * ss1s) / Z2
                        t3   = (ss1s - tj * ss2s) / Z2
                        t4   = (ss2s - tj * ss3s) / Z2
                        t5   = (ss3s - tj * ss4s) / Z2

                        lcount = 0
                        do l1 = 1, 3
                           t1  = Q(l1) - r(Nuc(ifunct),l1)
                           t2  = W(l1) - Q(l1)
                           ps  = t1 * sss  + t2 * ss1s
                           p1s = t1 * ss1s + t2 * ss2s
                           p2s = t1 * ss2s + t2 * ss3s
                           p3s = t1 * ss3s + t2 * ss4s
                           p4s = t1 * ss4s + t2 * ss5s
                           t6  = (p1s - tj * p2s) / Z2
                           t7  = (p2s - tj * p3s) / Z2

                           do l2 = 1, l1
                              t1   = Q(l2) - r(Nuc(ifunct),l2)
                              t2   = W(l2) - Q(l2)
                              pjs  = t1 * sss  + t2 * ss1s
                              pj1s = t1 * ss1s + t2 * ss2s
                              pj2s = t1 * ss2s + t2 * ss3s
                              pj3s = t1 * ss3s + t2 * ss4s
                              d1s  = t1 * p1s + t2 * p2s
                              d2s  = t1 * p2s + t2 * p3s
                              d3s  = t1 * p3s + t2 * p4s
                              t8   = (pj1s - tj * pj2s) / Z2
                              t9   = (pj2s - tj * pj3s) / Z2

                              f1 = 1.D0
                              if (l1 .eq. l2) then
                                 d1s = d1s + t3
                                 d2s = d2s + t4
                                 d3s = d3s + t5
                                 f1  = sq3
                              endif
                              t18 = (d1s - tj * d2s) / Z2

                              lij = 3
                              if (ifunct .eq. jfunct) then
                                 lij = l1
                              endif
                              do l3 = 1, lij
                                 t1    = Q(l3) - r(Nuc(jfunct),l3)
                                 t2    = W(l3) - Q(l3)
                                 d1pk  = t1 * d1s  + t2 * d2s
                                 d2pk  = t1 * d2s  + t2 * d3s
                                 pjpk  = t1 * pjs  + t2 * pj1s
                                 pj1pk = t1 * pj1s + t2 * pj2s
                                 pj2pk = t1 * pj2s + t2 * pj3s
                                 pipk  = t1 * ps   + t2 * p1s
                                 pi1pk = t1 * p1s  + t2 * p2s
                                 pi2pk = t1 * p2s  + t2 * p3s
                                 spk   = t1 * sss  + t2 * ss1s
                                 s1pk  = t1 * ss1s + t2 * ss2s
                                 s2pk  = t1 * ss2s + t2 * ss3s
                                 t10   = (s1pk - tj * s2pk) / Z2

                                 if (l1 .eq. l3) then
                                    d1pk  = d1pk  + t8
                                    d2pk  = d2pk  + t9
                                    pipk  = pipk  + ta
                                    pi1pk = pi1pk + t3
                                    pi2pk = pi2pk + t4
                                 endif
                                 if (l2 .eq. l3) then
                                    d1pk  = d1pk  + t6
                                    d2pk  = d2pk  + t7
                                    pjpk  = pjpk  + ta
                                    pj1pk = pj1pk + t3
                                    pj2pk = pj2pk + t4
                                 endif
                                 t16 = (pj1pk - tj * pj2pk) / Z2
                                 t17 = (pi1pk - tj * pi2pk) / Z2

                                 lk = l3
                                 if (ifunct .eq. jfunct) then
                                    lk = min(l3, Ll(l1) - Ll(l3) + l2)
                                 endif
                                 do l4 = 1, lk
                                    t1    = Q(l4) - r(Nuc(jfunct),l4)
                                    t2    = W(l4) - Q(l4)
                                    d1d   = t1 * d1pk + t2 * d2pk
                                    pjdkl = t1 * pj1pk + t2 * pj2pk
                                    pidkl = t1 * pi1pk + t2 * pi2pk
                                    d1pl  = t1 * d1s + t2 * d2s
                                    f2    = 1.D0

                                    if (l1 .eq. l4) then
                                       d1d   = d1d   + t16
                                       pidkl = pidkl + t10
                                       d1pl  = d1pl  + t8
                                    endif
                                    if (l2 .eq. l4) then
                                       d1d   = d1d   + t17
                                       pjdkl = pjdkl + t10
                                       d1pl  = d1pl  + t6
                                    endif
                                    if (l3 .eq. l4) then
                                       d1d   = d1d   + t18
                                       pjdkl = pjdkl + t8
                                       pidkl = pidkl + t6
                                       f2    = sq3
                                    endif
                                    t11 = pjdkl / Z2a
                                    t12 = pidkl / Z2a
                                    t13 = d1pl  / Z2a
                                    t14 = d1pk  / Z2a

                                    lcount = lcount +1
                                    do l5 = 1, 3
                                       term = (W(l5) - r(Nucd(kfunct),l5)) * d1d

                                       if (l1 .eq. l5) then
                                          term = term +t11
                                       endif
                                       if (l2 .eq. l5) then
                                          term = term +t12
                                       endif
                                       if (l3 .eq. l5) then
                                          term = term +t13
                                       endif
                                       if (l4 .eq. l5) then
                                          term = term +t14
                                       endif
                                       term = term * ccoef / (f1 * f2)

                                       if (rexp .lt. rmaxs) then
                                          if (jfunct .eq. ifunct) then
                                             cool_ind = (lcount + kknumd-22) * &
                                                        Md + kfunct + l5-1
                                          else
                                             cool_ind = (lcount + kknumd-37) * &
                                                        Md + kfunct + l5-1
                                          endif
                                          cool(cool_ind)=cool(cool_ind) + term
                                       else
                                          if (jfunct .eq. ifunct) then
                                             cool_ind = (lcount + kknums-22) * &
                                                        Md + kfunct + l5-1
                                          else
                                             cool_ind = (lcount + kknums-37) * &
                                                        Md + kfunct + l5-1
                                          endif
                                          cools(cool_ind) = cools(cool_ind) + &
                                                            real(term)
                                       endif
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo

                  ! (dd|d)
                  do kfunct = nsd+npd+1, Md, 6
                     dpc = (Q(1) -r(Nucd(kfunct),1))*(Q(1) -r(Nucd(kfunct),1)) &
                         + (Q(2) -r(Nucd(kfunct),2))*(Q(2) -r(Nucd(kfunct),2)) &
                         + (Q(3) -r(Nucd(kfunct),3))*(Q(3) -r(Nucd(kfunct),3))

                     do nck = 1, ncontd(kfunct)
                        ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                        t0    = ad(kfunct,nck) + Zij
                        Z2a   = 2.D0 * t0
                        Zc    = 2.D0 * ad(kfunct,nck)
                        ti    = Zij / t0
                        tj    = ad(kfunct,nck) / t0

                        W(1) = ti * Q(1) + tj * r(Nucd(kfunct),1)
                        W(2) = ti * Q(2) + tj * r(Nucd(kfunct),2)
                        W(3) = ti * Q(3) + tj * r(Nucd(kfunct),3)

                        t2   = sks / (ad(kfunct,nck) * dsqrt(t0))
                        uf   = tj * Zij*dpc
                        sss  = t2 * FUNCT(0,uf)
                        ss1s = t2 * FUNCT(1,uf)
                        ss2s = t2 * FUNCT(2,uf)
                        ss3s = t2 * FUNCT(3,uf)
                        ss4s = t2 * FUNCT(4,uf)
                        ss5s = t2 * FUNCT(5,uf)
                        ss6s = t2 * FUNCT(6,uf)
                        t3   = (sss  - tj * ss1s) / Z2
                        t4   = (ss1s - tj * ss2s) / Z2
                        t5   = (ss2s - tj * ss3s) / Z2
                        t6   = (ss3s - tj * ss4s) / Z2
                        t6b  = (ss4s - tj * ss5s) / Z2

                        lcount = 0
                        do l1 = 1, 3
                           t1  = Q(l1) - r(Nuc(ifunct),l1)
                           t2  = W(l1) - Q(l1)
                           ps  = t1 * sss  + t2 * ss1s
                           p1s = t1 * ss1s + t2 * ss2s
                           p2s = t1 * ss2s + t2 * ss3s
                           p3s = t1 * ss3s + t2 * ss4s
                           p4s = t1 * ss4s + t2 * ss5s
                           p5s = t1 * ss5s + t2 * ss6s
                           t7  = (ps  - tj * p1s) / Z2
                           t8  = (p1s - tj * p2s) / Z2
                           t9  = (p2s - tj * p3s) / Z2
                           t10 = (p3s - tj * p4s) / Z2

                           do l2 = 1, l1
                              t1   = Q(l2) - r(Nuc(ifunct),l2)
                              t2   = W(l2) - Q(l2)
                              pjs  = t1 * sss  + t2 * ss1s
                              pj1s = t1 * ss1s + t2 * ss2s
                              pj2s = t1 * ss2s + t2 * ss3s
                              pj3s = t1 * ss3s + t2 * ss4s
                              pj4s = t1 * ss4s + t2 * ss5s
                              ds   = t1 * ps   + t2 * p1s
                              d1s  = t1 * p1s  + t2 * p2s
                              d2s  = t1 * p2s  + t2 * p3s
                              d3s  = t1 * p3s  + t2 * p4s
                              d4s  = t1 * p4s  + t2 * p5s
                              t11  = (pjs  - tj * pj1s) / Z2
                              t12  = (pj1s - tj * pj2s) / Z2
                              t13  = (pj2s - tj * pj3s) / Z2
                              t14  = (pj3s - tj * pj4s) / Z2

                              f1 = 1.D0
                              if (l1 .eq. l2) then
                                 ds  = ds  + t3
                                 d1s = d1s + t4
                                 d2s = d2s + t5
                                 d3s = d3s + t6
                                 d4s = d4s + t6b
                                 f1  = sq3
                              endif
                              t16  = (ds  - tj * d1s) / Z2
                              t17  = (d1s - tj * d2s) / Z2
                              t18  = (d2s - tj * d3s) / Z2
                              t22a = d2s / Z2a

                              lij = 3
                              if (ifunct .eq. jfunct) then
                                 lij = l1
                              endif
                              do l3 = 1, lij
                                 t1    = Q(l3) - r(Nuc(jfunct),l3)
                                 t2    = W(l3) - Q(l3)
                                 dpk   = t1 * ds   + t2 * d1s
                                 d1pk  = t1 * d1s  + t2 * d2s
                                 d2pk  = t1 * d2s  + t2 * d3s
                                 d3pk  = t1 * d3s  + t2 * d4s
                                 pjpk  = t1 * pjs  + t2 * pj1s
                                 pj1pk = t1 * pj1s + t2 * pj2s
                                 pj2pk = t1 * pj2s + t2 * pj3s
                                 pj3pk = t1 * pj3s + t2 * pj4s
                                 pipk  = t1 * ps   + t2 * p1s
                                 pi1pk = t1 * p1s  + t2 * p2s
                                 pi2pk = t1 * p2s  + t2 * p3s
                                 pi3pk = t1 * p3s  + t2 * p4s
                                 spk   = t1 * sss  + t2 * ss1s
                                 s1pk  = t1 * ss1s + t2 * ss2s
                                 s2pk  = t1 * ss2s + t2 * ss3s
                                 s3pk  = t1 * ss3s + t2 * ss4s
                                 t15   = (s2pk - tj * s3pk) / Z2

                                 if (l1 .eq. l3) then
                                    dpk   = dpk   + t11
                                    d1pk  = d1pk  + t12
                                    d2pk  = d2pk  + t13
                                    d3pk  = d3pk  + t14
                                    pipk  = pipk  + t3
                                    pi1pk = pi1pk + t4
                                    pi2pk = pi2pk + t5
                                    pi3pk = pi3pk + t6
                                 endif
                                 if (l2 .eq. l3) then
                                    dpk   = dpk   + t7
                                    d1pk  = d1pk  + t8
                                    d2pk  = d2pk  + t9
                                    d3pk  = d3pk  + t10
                                    pjpk  = pjpk  + t3
                                    pj1pk = pj1pk + t4
                                    pj2pk = pj2pk + t5
                                    pj3pk = pj3pk + t6
                                 endif
                                 t20 = pj2pk / Z2a
                                 t21 = pi2pk / Z2a
                                 t22 = d2pk  / Z2a
                                 t24 = (pjpk  - tj * pj1pk) / Z2
                                 t25 = (pj1pk - tj * pj2pk) / Z2
                                 t26 = (pj2pk - tj * pj3pk) / Z2
                                 t27 = (pipk  - tj * pi1pk) / Z2
                                 t28 = (pi1pk - tj * pi2pk) / Z2
                                 t29 = (pi2pk - tj * pi3pk) / Z2

                                 lk = l3
                                 if (ifunct .eq. jfunct) then
                                    lk = min(l3, Ll(l1) - Ll(l3) + l2)
                                 endif
                                 do l4 = 1, lk
                                    t1  = Q(l4) - r(Nuc(jfunct),l4)
                                    t2  = W(l4) - Q(l4)
                                    d0d   = t1 * dpk   + t2 * d1pk
                                    d1d   = t1 * d1pk  + t2 * d2pk
                                    d2d   = t1 * d2pk  + t2 * d3pk
                                    pjdkl = t1 * pj2pk + t2 * pj3pk
                                    pidkl = t1 * pi2pk + t2 * pi3pk
                                    d2pl  = t1 * d2s   + t2 * d3s
                                    sdkl  = t1 * s2pk  + t2 * s3pk
                                    pj2pl = t1 * pj2s  + t2 * pj3s
                                    pi2pl = t1 * p2s   + t2 * p3s

                                    f2 = 1.D0
                                    if (l1 .eq. l4) then
                                       d0d   = d0d   + t24
                                       d1d   = d1d   + t25
                                       d2d   = d2d   + t26
                                       pidkl = pidkl + t15
                                       d2pl  = d2pl  + t13
                                       pi2pl = pi2pl + t5
                                    endif
                                    if (l2 .eq. l4) then
                                       d0d   = d0d   + t27
                                       d1d   = d1d   + t28
                                       d2d   = d2d   + t29
                                       pjdkl = pjdkl + t15
                                       d2pl  = d2pl  + t9
                                       pj2pl = pj2pl + t5
                                    endif
                                    if (l3 .eq. l4) then
                                       sdkl  = sdkl  + t5
                                       d0d   = d0d   + t16
                                       d1d   = d1d   + t17
                                       d2d   = d2d   + t18
                                       pjdkl = pjdkl + t13
                                       pidkl = pidkl + t9
                                       f2    = sq3
                                    endif
                                    t30 = pjdkl / Z2a
                                    t40 = pidkl / Z2a
                                    t50 = sdkl  / Z2a
                                    t60 = pj2pl / Z2a
                                    t70 = pi2pl / Z2a
                                    t80 = d2pl  / Z2a
                                    t23 = (d0d - ti * d1d) / Zc

                                    lcount = lcount +1
                                    do l5 = 1, 3
                                       t1     = W(l5) - r(Nucd(kfunct),l5)
                                       ddp    = t1 * d2d
                                       pjdklp = t1 * pjdkl
                                       pidklp = t1 * pidkl
                                       dijplp = t1 * d2pl
                                       dijpkp = t1 * d2pk

                                       if (l1 .eq. l5) then
                                          ddp    = ddp    + t30
                                          pidklp = pidklp + t50
                                          dijplp = dijplp + t60
                                          dijpkp = dijpkp + t20
                                       endif
                                       if (l2 .eq. l5) then
                                          ddp    = ddp    + t40
                                          pjdklp = pjdklp + t50
                                          dijplp = dijplp + t70
                                          dijpkp = dijpkp + t21
                                       endif
                                       if (l3 .eq. l5) then
                                          ddp    = ddp    + t80
                                          pjdklp = pjdklp + t60
                                          pidklp = pidklp + t70
                                          dijpkp = dijpkp + t22a
                                       endif
                                       if (l4 .eq. l5) then
                                          ddp    = ddp    + t22
                                          pjdklp = pjdklp + t20
                                          pidklp = pidklp + t21
                                          dijplp = dijplp + t22a
                                       endif
                                       t31 = pjdklp / Z2a
                                       t41 = pidklp / Z2a
                                       t51 = dijplp / Z2a
                                       t61 = dijpkp / Z2a

                                       do l6 = 1, l5
                                          term = (W(l6)-r(Nucd(kfunct),l6)) *ddp

                                          f3 = 1.D0
                                          if (l1 .eq. l6) then
                                             term = term + t31
                                          endif
                                          if (l2 .eq. l6) then
                                             term = term + t41
                                          endif
                                          if (l3 .eq. l6) then
                                             term = term + t51
                                          endif
                                          if (l4 .eq. l6) then
                                             term = term + t61
                                          endif
                                          if (l5 .eq. l6) then
                                             term = term + t23
                                             f3   = sq3
                                          endif
                                          term = term * ccoef / (f1 * f2 * f3)

                                          l56 = Ll(l5) + l6
                                          if (rexp .lt. rmaxs) then
                                             if (jfunct .eq. ifunct) then
                                                cool_ind = (lcount + kknumd-22)&
                                                         * Md + kfunct + l56-1
                                             else
                                                cool_ind = (lcount + kknumd-37)&
                                                         * Md + kfunct + l56-1
                                             endif
                                             cool(cool_ind) = cool(cool_ind)   &
                                                             + term
                                          else
                                             if (jfunct .eq. ifunct) then
                                                cool_ind = (lcount + kknums-22)&
                                                         * Md + kfunct + l56-1
                                             else
                                                cool_ind = (lcount + kknums-37)&
                                                         * Md + kfunct + l56-1
                                             endif
                                             cools(cool_ind) = cools(cool_ind) &
                                                             + real(term)
                                          endif
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
