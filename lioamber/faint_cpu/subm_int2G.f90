!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INT2G %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! 2e integral gradients, 2 indexes: density fitting functions.                 !
!                                                                              !
! EXTERNAL INPUT: system information.                                          !
!   · natom: number of QM atoms.                                               !
!   · ntatom: total number of atoms (QM+MM)                                    !
!   · r(ntatom,3): atoms' coordinates.                                         !
!   · d(natom,natom): distances between QM atoms.                              !
!                                                                              !
! INTERNAL INPUT: basis set information.                                       !
!   · M: number of basis functions (without contractions)                      !
!   · Md: number of auxiliary basis functions (without contractions)           !
!   · ncontd(Md): number of contractions per auxiliary function.               !
!   · ad(Md,nl): auxiliary basis function exponents.                           !
!   · cd(Md,nl): auxiliary basis function coefficients.                        !
!   · nshelld(0:3): number of auxiliary basis functions per shell (s,p,d).     !
!   · Nucd(M): atomic index corresponding to auxiliary function i.             !
!   · af(Md): variational coefficient for auxiliary function i.                !
!   · NORM: use custom normalization (now default and deprecated option)       !
!                                                                              !
! OUTPUTS:                                                                     !
!   · f(natom,3): QM gradients (= -forces)                                     !
!                                                                              !
! Original and debugged (or supposed to): Dario Estrin Jul/1992                !
! Refactored:                             Federico Pedron Sep/2018             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_int2G
contains
subroutine int2G(f, natom, ntatom, r, d)

   use liosubs_math , only: FUNCT
   use basis_data   , only: M, Md, NORM, af, nshelld, ad, Nucd, ncontd, cd
   use constants_mod, only: pi5

   implicit none
   ! aux . things
   integer         , intent(in)    :: natom, ntatom
   double precision, intent(in)    :: r(ntatom,3), d(natom,natom)
   double precision, intent(inout) :: f(natom,3)

   double precision  :: f1, f2, cci, ccj, ccoef, ccoef2, SQ3, uf,  Z2, Zij, Zc,&
                        Zc2, Q(3)
   double precision  :: sp, spj, sp1j, s1pk, s2pk, s0s, s1s, s2s, s3s, s4s, s5s
   double precision  :: ps, pp, pd, pis, pip, pid, pjs, pjp, pjd, pi0s, pi0p, &
                        pi0d, pj0s, pj0p, pj0d, p1p, pi2s, pi2p, pi3s, pi4s,  &
                        pj2s, pj2p, pj3s
   double precision  :: ds, dp, dpl, dsd, d0p, d0pl, d1d, d1p, d1s, d2s, d3s, &
                        df
   double precision  :: fs, fp, fd
   double precision  :: ti, tj, t0, t1, t2, t10, t11, t12, t12b, t13, t13a,    &
                        t13b, t14, t14a, t14b, t15, t15a, t15b, t16, t16b, t17,&
                        t17b, t18, t18b, t20, t21, t22, t23, t24, t25, t26,    &
                        t27, t27a, t30, t31, t32, t33, t40, t41, t42, t43

   integer           :: nsd, npd, ndd, nci, ncj, ifunct, jfunct, af_i, af_j,   &
                        lk, lij, l1, l2, l3, l4, l5, Ll(3)

   SQ3 = 1.0D0
   if (NORM) SQ3 = sqrt(3.D0)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) / 2
   enddo
   nsd = nshelld(0); npd = nshelld(1); ndd = nshelld(2)

   ! (s|s)
   do ifunct = 1, nsd
   do jfunct = 1, ifunct
      f1 = 2.0D0
      if (ifunct .eq. jfunct) f1 = 1.0D0

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         t2  = pi5 / (sqrt(Zij) * t0)

         ti   = ad(ifunct,nci) / Zij
         tj   = ad(jfunct,ncj) / Zij
         Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
         Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
         Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

         uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
         s1s = t2 * FUNCT(1,uf)

         ccoef = - cd(ifunct,nci) * cd(jfunct,ncj) * f1 * af(ifunct) *af(jfunct)
         do l1 = 1, 3
            t1 = Q(l1) - r(Nucd(ifunct),l1)
            t2 = Q(l1) - r(Nucd(jfunct),l1)
            ps = t1 * s1s
            sp = t2 * s1s
            f(Nucd(ifunct),l1) = f(Nucd(ifunct),l1) + ad(ifunct,nci) * ccoef *ps
            f(Nucd(jfunct),l1) = f(Nucd(jfunct),l1) + ad(jfunct,ncj) * ccoef *sp
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|s)
   do ifunct = nsd+1, nsd+npd, 3
   do jfunct = 1    , nsd
      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         t2  = pi5 / (sqrt(Zij) * t0)

         ti   = ad(ifunct,nci) / Zij
         tj   = ad(jfunct,ncj) / Zij
         Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
         Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
         Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

         uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
         s0s = t2 * FUNCT(0,uf)
         s1s = t2 * FUNCT(1,uf)
         s2s = t2 * FUNCT(2,uf)
         t10 = 0.5D0 * s1s / Zij
         t20 = 0.5D0 *(s0s - tj * s1s) / ad(ifunct,nci)

         ccoef  = cd(ifunct,nci) * cd(jfunct,ncj)
         do l1 = 1, 3
            ps     = (Q(l1) - r(Nucd(ifunct),l1)) * s2s
            ccoef2 = - ccoef * af(jfunct) * af(ifunct + l1-1)

            do l2 = 1, 3
               pp = (Q(l2) - r(Nucd(jfunct),l2)) * ps
               ds = (Q(l2) - r(Nucd(ifunct),l2)) * ps
               if (l1 .eq. l2) then
                  pp = pp  + t10
                  ds = ds  + t20
                  f(Nucd(ifunct),l2) = f(Nucd(ifunct),l2) - ccoef2 * s0s
               endif

               f(Nucd(jfunct),l2) = f(Nucd(jfunct),l2) + 2.0D0 * ccoef2 * pp * &
                                    ad(jfunct,ncj)
               f(Nucd(ifunct),l2) = f(Nucd(ifunct),l2) + 2.0D0 * ccoef2 * ds * &
                                    ad(ifunct,nci)
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|p)
   do ifunct = nsd+1, nsd+npd, 3
   do jfunct = nsd+1, ifunct , 3
      do nci = 1, ncontd(ifunct)
         do ncj = 1, ncontd(jfunct)
            Zij = ad(ifunct,nci) + ad(jfunct,ncj)
            Z2  = 2.0D0 * Zij
            t0  = ad(ifunct,nci) * ad(jfunct,ncj)
            t2  = pi5 / (sqrt(Zij) * t0)

            ti   = ad(ifunct,nci) / Zij
            tj   = ad(jfunct,ncj) / Zij
            Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
            Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
            Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

            uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
            s1s = t2 * FUNCT(1,uf)
            s2s = t2 * FUNCT(2,uf)
            s3s = t2 * FUNCT(3,uf)
            t25 = s2s  / Z2

            ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
            do l1 = 1, 3
               t1   = Q(l1) - r(Nucd(ifunct),l1)
               pi0s = t1 * s1s
               pis  = t1 * s2s
               pi2s = t1 * s3s
               t15  = 0.5D0 * (pi0s - ti * pis) / ad(jfunct,ncj)
               t17  = pis / Z2
               t20  = pis / Z2

               lij = 3
               if (ifunct .eq. jfunct) lij = l1
               do l2 = 1, lij
                  t2   = Q(l2)-r(Nucd(jfunct),l2)
                  spj  = t2 * s1s
                  sp1j = t2 * s2s
                  p1p  = t2 * pi2s
                  t20  = 0.5D0 * (spj - tj * sp1j) / ad(ifunct,nci)
                  t21  = sp1j / Z2
                  if (l1 .eq. l2) p1p = p1p + t25

                  af_i = ifunct + l1 -1
                  af_j = jfunct + l2 -1
                  f1   = 1.0D0
                  if (af_i .eq. af_j) f1 = 0.5D0
                  ccoef2 = - af(af_i) * af(af_j) * f1 * ccoef

                  do l3 = 1, 3
                     pd = (Q(l3) - r(Nucd(jfunct),l3)) * p1p
                     dp = (Q(l3) - r(Nucd(ifunct),l3)) * p1p

                     if (l1 .eq. l3) then
                        dp = dp + t20
                        pd = pd + t21
                        f(Nucd(ifunct),l3) = f(Nucd(ifunct),l3) - ccoef2 * spj
                     endif
                     if (l2 .eq. l3) then
                        dp = dp + t17
                        pd = pd + t15
                        f(Nucd(jfunct),l3) = f(Nucd(jfunct),l3) - ccoef2 * pi0s
                     endif

                     f(Nucd(ifunct),l3) = f(Nucd(ifunct),l3) + 2.D0 * ccoef2 * &
                                          dp * ad(ifunct,nci)
                     f(Nucd(jfunct),l3) = f(Nucd(jfunct),l3) + 2.D0 * ccoef2 * &
                                          pd *ad(jfunct,ncj)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   enddo

   ! (d|s)
   do ifunct = nsd+npd+1, Md , 6
   do jfunct = 1        , nsd
      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         Z2  = 2.0D0 * Zij
         Zc  = 2.0D0 * ad(ifunct,nci)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         t2  = pi5 / (sqrt(Zij) * t0)

         ti   = ad(ifunct,nci) / Zij
         tj   = ad(jfunct,ncj) / Zij
         Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
         Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
         Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

         uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
         s1s = t2 * FUNCT(1,uf)
         s2s = t2 * FUNCT(2,uf)
         s3s = t2 * FUNCT(3,uf)
         t10 = (s1s - t0 * s2s / (Zij * ad(ifunct,nci))) / Zc

         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
         do l1 = 1, 3
            t1   = Q(l1) - r(Nucd(ifunct),l1)
            pi0s = t1 * s1s
            pis  = t1 * s2s
            pi2s = t1 * s3s
            t12  = pis / Z2
            t22  = (pi0s - tj * pis) / Zc

            do l2 = 1, l1
               t1   = Q(l2) - r(Nucd(ifunct),l2)
               pj0s = t1 * s1s
               pjs  = t1 * s2s
               ds   = t1 * pi2s
               t11  = pjs / Z2
               t21  = (pj0s - tj * pjs) / Zc
               f1   = 1.0D0

               if (l1 .eq. l2) then
                  f1 = SQ3
                  ds = ds + t10
               endif

               ccoef2 = -af(ifunct + Ll(l1) + l2 - 1) * af(jfunct) * ccoef / f1
               cci    = 2.0D0 * ccoef2 * ad(ifunct,nci)
               ccj    = 2.0D0 * ccoef2 * ad(jfunct,ncj)
               do l3 = 1, 3
                  dp = (Q(l3) - r(Nucd(jfunct),l3)) * ds
                  fs = (Q(l3) - r(Nucd(ifunct),l3)) * ds
                  if (l1 .eq. l3) then
                     dp = dp + t11
                     fs = fs + t21
                     f(Nucd(ifunct),l3) = f(Nucd(ifunct),l3) - ccoef2 * pj0s
                  endif
                  if (l2 .eq. l3) then
                     dp = dp + t12
                     fs = fs + t22
                     f(Nucd(ifunct),l3) = f(Nucd(ifunct),l3) - ccoef2 * pi0s
                  endif
                  f(Nucd(ifunct),l3) = f(Nucd(ifunct),l3) + cci * fs
                  f(Nucd(jfunct),l3) = f(Nucd(jfunct),l3) + ccj * dp
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|p)
   do ifunct = nsd+npd+1, Md     , 6
   do jfunct = nsd+1    , nsd+npd, 3
      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         Z2  = 2.0D0 * Zij
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         t2  = pi5 / (sqrt(Zij) * t0)
         Zc  = 2.0D0 * ad(ifunct,nci)
         Zc2 = 2.0D0 * ad(jfunct,ncj)

         ti   = ad(ifunct,nci) / Zij
         tj   = ad(jfunct,ncj) / Zij
         Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
         Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
         Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

         uf   = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
         s0s  = t2 * FUNCT(0,uf)
         s1s  = t2 * FUNCT(1,uf)
         s2s  = t2 * FUNCT(2,uf)
         s3s  = t2 * FUNCT(3,uf)
         s4s  = t2 * FUNCT(4,uf)
         t13a = s1s / Z2
         t13  = s2s / Z2
         t10  = (s0s - tj * s1s) / Zc
         t11  = (s1s - tj * s2s) / Zc
         t12  = (s2s - tj * s3s) / Zc

         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
         do l1 = 1, 3
            t1   = Q(l1) - r(Nucd(ifunct),l1)
            pis  = t1 * s2s
            pi2s = t1 * s3s
            pi3s = t1 * s4s
            t14  = pi2s / Z2

            do l2 = 1, l1
               t1   = Q(l2) - r(Nucd(ifunct),l2)
               pjs  = t1 * s2s
               pj2s = t1 * s3s
               t15  = pj2s / Z2
               ds   = t1 * pis
               d1s  = t1 * pi2s
               d2s  = t1 * pi3s
               f1   = 1.0D0

               if (l1 .eq. l2) then
                  ds  = ds + t10
                  d1s = d1s + t11
                  d2s = d2s + t12
                  f1  = SQ3
               endif
               t18 = (ds - ti * d1s) / Zc2
               t23 = d1s / Z2

               do l3 = 1, 3
                  t0    = Q(l3) - r(Nucd(jfunct),l3)
                  dp    = t0 * d2s
                  pip   = t0 * pi2s
                  pi0p  = t0 * pis
                  pjp   = t0* pj2s
                  pj0p  = t0 * pjs

                  if (l1 .eq. l3) then
                     dp   = dp   + t15
                     pip  = pip  + t13
                     pi0p = pi0p + t13a
                  endif
                  if (l2 .eq. l3) then
                     dp   = dp   + t14
                     pjp  = pjp  + t13
                     pj0p = pj0p + t13a
                  endif
                  t16 = pip / Z2
                  t17 = pjp / Z2
                  t21 = (pj0p - tj * pjp) / Zc
                  t22 = (pi0p - tj * pip) / Zc

                  ccoef2 = - af(ifunct + Ll(l1) + l2 -1) * af(jfunct + l3 -1) *&
                           ccoef/ f1
                  cci    = 2.0D0 * ccoef2 * ad(ifunct,nci)
                  ccj    = 2.0D0 * ccoef2 * ad(jfunct,ncj)

                  do l4 = 1, 3
                     t0  = Q(l4) - r(Nucd(jfunct),l4)
                     t1  = Q(l4) - r(Nucd(ifunct),l4)
                     dsd = t0 * dp
                     fp  = t1 * dp

                     if (l1 .eq. l4) then
                        dsd = dsd + t17
                        fp  = fp  + t21
                        f(Nucd(ifunct),l4) = f(Nucd(ifunct),l4) - ccoef2 * pj0p
                     endif
                     if (l2 .eq. l4) then
                        dsd = dsd + t16
                        fp  = fp  + t22
                        f(Nucd(ifunct),l4) = f(Nucd(ifunct),l4) - ccoef2 * pi0p
                     endif
                     if (l3.eq.l4) then
                        f(Nucd(jfunct),l4) = f(Nucd(jfunct),l4) - ccoef2 * ds
                        dsd = dsd + t18
                        fp  = fp  + t23
                     endif
                     f(Nucd(ifunct),l4) = f(Nucd(ifunct),l4) + cci * fp
                     f(Nucd(jfunct),l4) = f(Nucd(jfunct),l4) + ccj * dsd
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|d)
   do ifunct = nsd+npd+1, Md    , 6
   do jfunct = nsd+npd+1, ifunct, 6
         do nci = 1, ncontd(ifunct)
            do ncj = 1, ncontd(jfunct)
               Zij = ad(ifunct,nci) + ad(jfunct,ncj)
               Z2  = 2.0D0 * Zij
               Zc  = 2.0D0 * ad(ifunct,nci)
               Zc2 = 2.0D0 * ad(jfunct,ncj)
               t0  = ad(ifunct,nci) * ad(jfunct,ncj)
               t2  = pi5 / (sqrt(Zij) * t0)

               ti   = ad(ifunct,nci) / Zij
               tj   = ad(jfunct,ncj) / Zij
               Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
               Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
               Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

               uf   = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
               s0s  = t2 * FUNCT(0,uf)
               s1s  = t2 * FUNCT(1,uf)
               s2s  = t2 * FUNCT(2,uf)
               s3s  = t2 * FUNCT(3,uf)
               s4s  = t2 * FUNCT(4,uf)
               s5s  = t2 * FUNCT(5,uf)
               t13a = s1s / Z2
               t13  = s2s / Z2
               t13b = s3s / Z2
               t10  = (s0s - tj * s1s) / Zc
               t11  = (s1s - tj * s2s) / Zc
               t12  = (s2s - tj * s3s) / Zc
               t12b = (s3s - tj * s4s) / Zc

               ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
               do l1 = 1, 3
                  t1   = Q(l1) - r(Nucd(ifunct),l1)
                  pi0s = t1 * s1s
                  pis  = t1 * s2s
                  pi2s = t1 * s3s
                  pi3s = t1 * s4s
                  pi4s = t1 * s5s
                  t14a = pis  / Z2
                  t14b = pi3s / Z2
                  t14  = pi2s / Z2
                  t25  = (pi0s - ti * pis)  / Zc2
                  t26  = (pis  - ti * pi2s) / Zc2

                  do l2 = 1, l1
                     t1   = Q(l2) - r(Nucd(ifunct),l2)
                     pj0s = t1 * s1s
                     pjs  = t1 * s2s
                     pj2s = t1 * s3s
                     pj3s = t1 * s4s
                     ds   = t1 * pis
                     d1s  = t1 * pi2s
                     d2s  = t1 * pi3s
                     d3s  = t1 * pi4s
                     t15a = pjs / Z2
                     t15b = pj3s / Z2
                     t15  = pj2s / Z2
                     t23  = (pj0s - ti * pjs)  / Zc2
                     t24  = (pjs  - ti * pj2s) / Zc2
                     f1   = 1.0D0

                     if (l1 .eq. l2) then
                        ds  = ds + t10
                        d1s = d1s + t11
                        d2s = d2s + t12
                        d3s = d3s + t12b
                        f1  = SQ3
                     endif
                     t18b = (d1s - ti * d2s) / Zc2

                     lij = 3
                     if (ifunct .eq. jfunct) lij = l1
                     do l3 = 1, lij
                        t0 = Q(l3) - r(Nucd(jfunct),l3)
                        s1pk = t0 * s2s
                        s2pk = t0 * s3s
                        pi2p = t0 * pi3s
                        pip  = t0 * pi2s
                        pi0p = t0 * pis
                        pj2p = t0 * pj3s
                        pjp  = t0 * pj2s
                        pj0p = t0 * pjs
                        d0p  = t0 * d1s
                        dp   = t0 * d2s
                        d1p  = t0 * d3s
                        t27a = s1pk / Z2
                        t27  = s2pk / Z2

                        if (l1 .eq. l3) then
                           d0p  = d0p  + t15a
                           dp   = dp   + t15
                           d1p  = d1p  + t15b
                           pi2p = pi2p + t13b
                           pip  = pip  + t13
                           pi0p = pi0p + t13a
                        endif
                        if (l2 .eq. l3) then
                           d0p  = d0p  + t14a
                           dp   = dp   + t14
                           d1p  = d1p  + t14b
                           pj2p = pj2p + t13b
                           pjp  = pjp  + t13
                           pj0p = pj0p + t13a
                        endif
                        t16b = pi2p / Z2
                        t17b = pj2p / Z2
                        t33  = dp   / Z2
                        t43  = (d0p - ti * dp) / Zc2

                        lk = l3
                        if (ifunct .eq. jfunct) lk = min(l3, Ll(l1) - Ll(l3)+l2)
                        do l4 = 1, lk
                           t0   = Q(l4) - r(Nucd(jfunct),l4)
                           d1d  = t0 * d1p
                           d0pl = t0 * d1s
                           dpl  = t0 * d2s
                           pi0d = t0 * pip
                           pid  = t0 * pi2p
                           pj0d = t0 * pjp
                           pjd  = t0 * pj2p
                           f2   = 1.0D0

                           if (l1 .eq. l4) then
                              d1d  = d1d  + t17b
                              d0pl = d0pl + t15a
                              dpl  = dpl  + t15
                              pi0d = pi0d + t27a
                              pid  = pid  + t27
                           endif
                           if (l2 .eq. l4) then
                              d1d  = d1d  + t16b
                              d0pl = d0pl + t14a
                              dpl  = dpl  + t14
                              pj0d = pj0d + t27a
                              pjd  = pjd  + t27
                           endif
                           if (l3 .eq. l4) then
                              d1d  = d1d  + t18b
                              pi0d = pi0d + t25
                              pid  = pid  + t26
                              pj0d = pj0d + t23
                              pjd  = pjd  + t24
                              f2   = SQ3
                           endif
                           t30 = (pj0d - tj * pjd) / Zc
                           t31 = (pi0d - tj * pid) / Zc
                           t32 = dpl / Z2
                           t40 = pjd / Z2
                           t41 = pid / Z2
                           t42 = (d0pl - ti * dpl) / Zc2

                           af_i = ifunct + Ll(l1) + l2 -1
                           af_j = jfunct + Ll(l3) + l4 -1
                           ccoef2 = 1.0D0
                           if (af_i .eq. af_j) ccoef2 = 0.5D0

                           ccoef2 = - ccoef2 * af(af_i) * af(af_j) * ccoef / &
                                    (f1 * f2)
                           cci = 2.0D0 * ccoef2 * ad(ifunct,nci)
                           ccj = 2.0D0 * ccoef2 * ad(jfunct,ncj)

                           do l5 = 1, 3
                              df = (Q(l5) - r(Nucd(jfunct),l5)) * d1d
                              fd = (Q(l5) - r(Nucd(ifunct),l5)) * d1d

                              if (l1 .eq. l5) then
                                 df = df + t40
                                 fd = fd + t30
                                 f(Nucd(ifunct),l5) = f(Nucd(ifunct),l5) - &
                                                      ccoef2 * pj0d
                              endif
                              if (l2 .eq. l5) then
                                 df = df + t41
                                 fd = fd + t31
                                 f(Nucd(ifunct),l5) = f(Nucd(ifunct),l5) - &
                                                      ccoef2 * pi0d
                              endif
                              if (l3.eq.l5) then
                                 df = df + t42
                                 fd = fd + t32
                                 f(Nucd(jfunct),l5) = f(Nucd(jfunct),l5) - &
                                                      ccoef2 * d0pl
                              endif
                              if (l4.eq.l5) then
                                 df = df + t43
                                 fd = fd + t33
                                 f(Nucd(jfunct),l5) = f(Nucd(jfunct),l5) - &
                                                      ccoef2 * d0p
                              endif
                              f(Nucd(ifunct),l5) = f(Nucd(ifunct),l5) + cci * fd
                              f(Nucd(jfunct),l5) = f(Nucd(jfunct),l5) + ccj * df
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   return
end subroutine
end module subm_int2G
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
