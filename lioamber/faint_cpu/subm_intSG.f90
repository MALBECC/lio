!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates the gradients of overlap using the Obara-Saika recursive method.  !
! This is used in forces calculation.                                          !
!                                                                              !
! Input : basis function information, energy weighted density matrix.          !
! Output: gradients from overlap                                               !
!                                                                              !
! Debugged (or supposed to) 29-7-92 by Dario Estrin.                           !
! Refactored by Federico Pedron 08/2018                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_intSG
contains
subroutine intSG(ff, wgt_rho, r, d, natom, ntatom)
   use garcha_mod   , only: a, c, Nuc, ncont, nshell, M, Md, NORM
   use constants_mod, only: pi32
   implicit none

   integer         , intent(in)    :: natom, ntatom
   double precision, intent(in)    :: r(ntatom,3), d(natom,natom), wgt_rho(:)
   double precision, intent(inout) :: ff(natom,3)

   integer           :: ifunct, jfunct, en_wgt_ind, nci, ncj, lk, lij, l1, l2, &
                        l3, l4, l5, l12, l34, ns, np, nd, M2, M15, ll(3)

   double precision  :: ovlap, fsp, sq3, ccoef, rexp, Zij, Z2, fs, fd, f1, f2, &
                        ti, tj, te, t0, t1, t2, t4, t5, t10, t11, t12, t13,    &
                        t14, t15, t16, t17, ss, spi, spj, spk, ps, pp, pd,     &
                        pidkl, pipk, pis, pjdkl, pjpk, pjs, ds, dp, dd, df,    &
                        dsd, dijpk, dijpl, dijs, Q(3)

   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)
   ns = nshell(0); np = nshell(1); nd = nshell(2)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) / 2
   enddo

   M2 = 2 * M
   ! Energy weighted density matrix
   M15 = 1 + 2 * (M * (M +1)) + Md * (Md +1) + M

   ! (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))

         en_wgt_ind = ifunct + ((M2 - jfunct) * (jfunct -1)) / 2
         te = wgt_rho(en_wgt_ind) * ccoef * 2.0D0 * ss
         t4 = te * a(ifunct,nci)
         t5 = te * a(jfunct,ncj)

         do l1 = 1, 3
            t1 = Q(l1)             - r(Nuc(ifunct),l1)

            ff(Nuc(ifunct),l1) = ff(Nuc(ifunct),l1) + t4 * t1
            ff(Nuc(jfunct),l1) = ff(Nuc(jfunct),l1) + t5 * (t1 + &
                                  (r(Nuc(ifunct),l1) - r(Nuc(jfunct),l1)))
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1   , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            en_wgt_ind = ifunct + l1 -1 + ((M2 - jfunct) * (jfunct -1)) / 2
            t1 = (Q(l1) - r(Nuc(ifunct),l1))
            te = wgt_rho(en_wgt_ind) * ccoef  * ss
            t4 = 2.D0 * te * a(ifunct,nci)
            t5 = 2.D0 * te * a(jfunct,ncj)

            do l2 = 1, 3
               ds = (Q(l2) - r(Nuc(ifunct),l2)) * t1
               pp = (Q(l2) - r(Nuc(jfunct),l2)) * t1
               if (l1 .eq. l2) then
                  ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) - te
                  ds = ds + 1 / Z2
                  pp = pp + 1 / Z2
               endif
               ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) + t4 * ds
               ff(Nuc(jfunct),l2) = ff(Nuc(jfunct),l2) + t5 * pp
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|p)
   do ifunct = ns+1 , ns+np , 3
   do jfunct = ns+1 , ifunct, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            pis = ss * (Q(l1) - r(Nuc(ifunct),l1))
            t11 = pis / Z2

            lij = 3
            if (ifunct .eq. jfunct) then
               lij = l1
            endif
            do l2 = 1, lij
               t2  = Q(l2) - r(Nuc(jfunct),l2)
               spj = t2  * ss
               pp  = t2  * pis
               t13 = spj / Z2

               if (l1 .eq. l2) then
                  pp = pp + t10
               endif

               en_wgt_ind = ifunct + l1-1 + ((M2 - (jfunct + l2-1)) * &
                            (jfunct + l2-2)) / 2
               te = wgt_rho(en_wgt_ind) * ccoef
               t5 = te * 2.D0 * a(jfunct,ncj)
               t4 = te * 2.D0 * a(ifunct,nci)

               do l3 = 1, 3
                  dp = (Q(l3) - r(Nuc(ifunct),l3)) * pp
                  if (l1 .eq. l3) then
                     dp = dp + t13
                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * spj
                  endif
                  if (l2 .eq. l3) then
                     dp=dp+t11
                     ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) - te * pis
                  endif
                  pd = dp + (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) * pp
                  ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * dp
                  ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * pd
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|s)
   do ifunct = ns+np+1, M, 6
   do jfunct = 1      , ns
      dd=d(Nuc(ifunct),Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            pis = ss * (Q(l1) - r(Nuc(ifunct),l1))
            t12 = pis / Z2

            do l2 = 1, l1
               t1 = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = ss * t1
               dijs = t1 * pis
               t11  = pjs / Z2
               f1   = 1.D0

               if (l1 .eq. l2) then
                  f1 = sq3
                  dijs = dijs + t10
               endif

               l12 = l1 * (l1-1) / 2 + l2
               en_wgt_ind = ifunct + l12-1 + ((M2 - jfunct) * (jfunct-1)) / 2

               te = wgt_rho(en_wgt_ind) * ccoef / f1
               t4 = te * 2.D0 * a(ifunct,nci)
               t5 = te * 2.D0 * a(jfunct,ncj)
               do l3 = 1, 3
                  dp = (Q(l3) - r(Nuc(jfunct),l3)) * dijs

                  if (l1 .eq. l3) then
                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * pjs
                     dp = dp + t11
                  endif
                  if (l2 .eq. l3) then
                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * pis
                     dp = dp + t12
                  endif
                  fs = dp - (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) * dijs
                  ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * fs
                  ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * dp
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|p)
   do ifunct = ns+np+1, M    , 6
   do jfunct = ns+1   , ns+np, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            pis = ss * (Q(l1) - r(Nuc(ifunct),l1))
            t11 = pis / Z2

            do l2 = 1, l1
               t1   = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = t1 * ss
               dijs = t1 * pis
               t12  = pjs / Z2
               f1   = 1.D0

               if (l1 .eq. l2) then
                  f1   = sq3
                  dijs = dijs + t0
               endif
               t13 = dijs / Z2

               do l3 = 1, 3
                  t2 = Q(l3) - r(Nuc(jfunct),l3)
                  pipk  = t2 * pis
                  pjpk  = t2 * pjs
                  dijpk = t2 * dijs

                  if (l1 .eq. l3) then
                     pipk  = pipk  + t0
                     dijpk = dijpk + t12
                  endif
                  if (l2 .eq. l3) then
                     pjpk  = pjpk  + t0
                     dijpk = dijpk + t11
                  endif

                  l12 = l1 * (l1-1) / 2 + l2
                  en_wgt_ind = ifunct + l12-1 + ((M2 - (jfunct + l3-1)) * &
                                                (jfunct + l3-2)) / 2
                  te  = wgt_rho(en_wgt_ind) * ccoef / f1
                  t4  = te * 2.D0 * a(ifunct,nci)
                  t5  = te * 2.D0 * a(jfunct,ncj)
                  t14 = pipk / Z2
                  t15 = pjpk / Z2
                  do l4 = 1, 3
                     dsd = (Q(l4)-r(Nuc(jfunct),l4)) * dijpk

                     if (l1 .eq. l4) then
                        dsd = dsd + t15
                        ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) - te * pjpk
                     endif
                     if (l2 .eq. l4) then
                        dsd = dsd + t14
                        ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) - te * pipk
                     endif
                     if (l3 .eq. l4) then
                        dsd = dsd + t13
                        ff(Nuc(jfunct),l4) = ff(Nuc(jfunct),l4) - te * dijs
                     endif
                     fsp = dsd - (r(Nuc(ifunct),l4) - r(Nuc(jfunct),l4)) * dijpk
                     ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) + t4 * fsp
                     ff(Nuc(jfunct),l4) = ff(Nuc(jfunct),l4) + t5 * dsd
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|d)  gradients
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+np+1, ifunct, 6
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            pis = ss * (Q(l1) - r(Nuc(ifunct),l1))
            t17 = pis / Z2

            do l2 = 1, l1
               t1 = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = t1 * ss
               dijs = t1 * pis
               t16  = pjs / Z2
               f1   = 1.D0

               if (l1 .eq. l2) then
                  f1 = sq3
                  dijs = dijs + t0
               endif

               lij = 3
               if (ifunct .eq. jfunct) lij = l1
               do l3 = 1, lij
                  t2    = Q(l3) - r(Nuc(jfunct),l3)
                  spk   = t2 * ss
                  pipk  = t2 * pis
                  pjpk  = t2 * pjs
                  dijpk = t2 * dijs
                  t15   = spk / Z2

                  if (l1 .eq. l3) then
                     pipk  = pipk  + t0
                     dijpk = dijpk + t16
                  endif
                  if (l2 .eq. l3) then
                     pjpk  = pjpk  + t0
                     dijpk = dijpk + t17
                  endif
                  t13 = dijpk / Z2

                  lk = l3
                  if (ifunct .eq. jfunct) lk = min(l3, Ll(l1) - Ll(l3) + l2)
                  do l4 = 1, lk
                     f2 = 1.D0
                     t1 = Q(l4) - r(Nuc(jfunct),l4)
                     ovlap = t1 * dijpk
                     pjdkl = t1 * pjpk
                     pidkl = t1 * pipk
                     dijpl = t1 * dijs

                     if (l1 .eq. l4) then
                        pidkl = pidkl + t15
                        dijpl = dijpl + t16
                        ovlap = ovlap + pjpk / Z2
                     endif
                     if (l2 .eq. l4) then
                        pjdkl = pjdkl + t15
                        dijpl = dijpl + t17
                        ovlap = ovlap + pipk / Z2
                     endif
                     if (l3 .eq. l4) then
                        pjdkl = pjdkl + t16
                        pidkl = pidkl + t17
                        ovlap = ovlap + dijs / Z2
                        f2 = sq3
                     endif
                     t10 = pjdkl / Z2
                     t11 = pidkl / Z2
                     t12 = dijpl / Z2

                     l12 = l1 * (l1-1) / 2 + l2
                     l34 = l3 * (l3-1) / 2 + l4
                     en_wgt_ind = ifunct + l12-1 + ((M2 - (jfunct + l34-1)) * &
                                                    (jfunct + l34 -2)) / 2
                     te = wgt_rho(en_wgt_ind) * ccoef / (f1 * f2)
                     t4 = te * 2.D0 * a(ifunct,nci)
                     t5 = te * 2.D0 * a(jfunct,ncj)

                     do l5 = 1, 3
                        fd = (Q(l5)-r(Nuc(ifunct),l5)) * ovlap

                        if (l1 .eq. l5) then
                           fd = fd + t10
                           ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) - te * pjdkl
                        endif
                        if (l2 .eq. l5) then
                           fd = fd + t11
                           ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) - te * pidkl
                        endif
                        if (l3 .eq. l5) then
                           fd = fd + t12
                           ff(Nuc(jfunct),l5) = ff(Nuc(jfunct),l5) - te * dijpl
                        endif
                        if (l4 .eq. l5) then
                           fd = fd + t13
                           ff(Nuc(jfunct),l5) = ff(Nuc(jfunct),l5) - te * dijpk
                        endif
                        df = fd + (r(Nuc(ifunct),l5) - r(Nuc(jfunct),l5)) *ovlap
                        ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) + t4 * fd
                        ff(Nuc(jfunct),l5) = ff(Nuc(jfunct),l5) + t5 * df
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
end subroutine intSG
end module subm_intSG
