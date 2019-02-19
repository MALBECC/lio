!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INT3G %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! 2e integral gradients, 3 indexes: density fitting functions and wavefunction.!
!                                                                              !
! EXTERNAL INPUT: system information.                                          !
!   · calc_energy: boolean indicating whether to perform XC calculation.       !
!   · natom: number of QM atoms.                                               !
!   · ntatom: total number of atoms (QM+MM)                                    !
!   · r(ntatom,3): atoms' coordinates.                                         !
!   · d(natom,natom): distances between QM atoms.                              !
!   · rho_mat(M,M): density matrix.                                            !
!                                                                              !
! INTERNAL INPUT: basis set information.                                       !
!   · M: number of basis functions (without contractions)                      !
!   · ncont(M): number of contractions per function.                           !
!   · a(M,nl): basis function exponents.                                       !
!   · c(M,nl): basis function coefficients.                                    !
!   · nshell(0:3): number of basis functions per shell (s,p,d).                !
!   · Nuc(M): atomic index corresponding to function i.                        !
!   · Md: number of auxiliary basis functions (without contractions)           !
!   · ncontd(Md): number of contractions per auxiliary function.               !
!   · ad(Md,nl): auxiliary basis function exponents.                           !
!   · cd(Md,nl): auxiliary basis function coefficients.                        !
!   · nshelld(0:3): number of auxiliary basis functions per shell (s,p,d).     !
!   · Nucd(M): atomic index corresponding to auxiliary function i.             !
!   · af(Md): variational coefficient for auxiliary function i.                !
!   · NORM: use custom normalization (now default and deprecated option)       !
!   · rmax: cutoff value (in Ångström) for maximum exponent in the integrals.  !
!                                                                              !
! OUTPUTS:                                                                     !
!   · f(natom,3): QM gradients (= -forces)                                     !
!                                                                              !
! Original and debugged (or supposed to): Dario Estrin Jul/1992                !
! Refactored:                             Federico Pedron Sep/2018             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_int3G
contains
subroutine int3G(frc, calc_energy, rho_mat, r, d, natom, ntatom)

   use subm_int2G   , only: int2G
   use liosubs_math , only: FUNCT
   use basis_data   , only: a, c, nuc, nucd, ad, af, cd, ncont, ncontd, nshell,&
                            nshelld, M, Md, NORM, rmax
   use constants_mod, only: pi52

   implicit none

   logical         , intent(in)    :: calc_energy
   integer         , intent(in)    :: natom, ntatom
   double precision, intent(in)    :: rho_mat(:), r(ntatom,3), d(natom,natom)
   double precision, intent(inout) :: frc(natom,3)

   integer         , allocatable :: Jx(:), Ll(:)
   double precision, allocatable :: Q(:), W(:)

   ! ccoef is the coefficient product between function, Ci*Cj*Ck.
   ! f1, f2 and f3 are the normalization coefficients for d-type functions.
   ! rexp is the function exponent in distance units, while dpc stores distances
   ! between atoms.
   ! uf is the parameter used for Boys Function.
   double precision :: ccoef, rexp, SQ3, uf, dpc, f1, f2, f3

   ! ns, np and nd are the number of basis functions in each shell, while
   ! nsd, npd and ndd are used for the auxiliary basis.
   ! igpu checks whether the gradients are done in GPU or CPU.
   ! rho_ind and af_ind are vector-matrix indeces.
   integer :: ns, nsd, np, npd, nd, ndd, rho_ind, af_ind, igpu

   ! The following thousand variables store temporary results.
   double precision :: s0pk, pss, psf, s2dpm, s2dkl, s1pkpl, s1pk, s1ds, s1dpm,&
                       s2pl, s2pks, s2pkpl, s2pk, s2pjpk, s4pk, s3pl, s3pks,   &
                       s3pk, s3dkl, s2ds, sks, sp0d, sp0js, sp1d, sp1s, sp2js, &
                       sp3js, spd, spjpk, spjs, spk, spp, sps, ss0d, ss0p,     &
                       ss0pj, ss1d, ss1p, ss1pj, ss1pk, ss1s, ss2p, ss2pj,     &
                       ss2pk, ss2s, ss3s, ss4s, ss5s, ss6s, ss7s, ssd, ssf,    &
                       ssp, sspj, sspk, sss, spf
   double precision :: p1s, p2s, p3s, p4s, p5s, p6s, p0pk, p1pk, pi0dd, pi0d,  &
                       pds, pdp, pdd, pi0sd, pi0pp, pi0p, pi0dp, pi0dkl, pi1dp,&
                       pi1dkl, pi1dd, pi0spj, pi1pl, pi1pkpm, pi1pk, pi1d,     &
                       pi1spl, pi1sd, pi1pp, pi1plpm, pi1p, pi2pkpl, pi2pk,    &
                       pi2p, pi2dklp, pi2dkl, pi2spk, pi2spj, pi2pl, pi2pkpm,  &
                       pi3pk, pi3p, pi3dkl, pi2spl, pi2plpm, pij1s, pidklp,    &
                       pidkl, pi4pk, pi3pl, pip0d, pip, pijs, pij3s, pij2s,    &
                       pis2pk, pis1pk, pipkpl, pipk, pip1d, pj0dkl, pj0dd,     &
                       pj0d, pispk, pispj, pj0sd, pj0s, pj0pp, pj0p, pj0dp,    &
                       pjs, pjp0d, pj1p, pj1dp, pj1d, pj1dkl, pj4pk, pjpk,     &
                       pjp1d, pjp, pj1dd, pj1plpm, pj1pl, pj1pkpm, pj1pk,      &
                       pj1spl, pj1sd, pj1s, pj1pp, pj2pk, pj2p, pj2dklp, pj2pl,&
                       pj2pkpm, pj2pkpl, pj2dkl, pj3dkl, pj2spl, pj2s, pj2plpm,&
                       pj3s, pj3pl, pj3p, pj5s, pj4s, pj3pk, pjdklp, pjdkl,    &
                       pjpkpl, pjs1pk, pp0p, pp0d, pjs2pk, pp1p, psd, ps1d,    &
                       pp1s, pp0pl, pp2p, ppd, ppp, ppf, pps, ps, psp, ps0d,   &
                       pp1d, pp1pl
   double precision :: d0d, d0p, d0pk, d0pkd, d0pkp, d0pl, d0pld, d0plp, d0s,  &
                       d1d, d1p, d1pk, d1pkd, d1pkp, d1pl, d1pld, d1plp, d1pp, &
                       d1s, d1spm, d2d, d2p, dds, ddp, ddf, ddd, dij2plp,      &
                       dij2pkp, dfs, dfp, dp0p, dp, dijplp, dfd, dd2p, dijpkp, &
                       dp0pm, dp1d, dp1p, dp1pm, dp1s, dp2p, dpd, dpf, dpk,    &
                       dpp, dps, ds, ds0p, ds1d, ds1p, dsp, dsf, dsd, ds2pl,   &
                       ds2p, ds1s, dss, dspl, dd1pn, dd1s, dd1p, dd1d, dd0pn,  &
                       dd0p, dd, d5s, d4s, d3s, d3pl, d3pk, d3p, d4pk, d3d,    &
                       d2spm, d2s, d2pl, d2pk
   double precision :: fdp, fdd, fss, fsp, fsd, fps, fpp, fpd, fds

   double precision :: ta, tb, ti, tj, te, ty, t0, t1, t2, t2a, t2b, t3, t3a,  &
                       t3b, t4, t4b,  t5, t5a, t5b, t5x, t5y, t6, t6a, t6b,    &
                       t6c, t6d, t7, t7a, t7b, t7c, t7d, t8, t8a, t8b, t9, t9b,&
                       t10, t10a, t10b, t11, t11a, t11b, t12, t12a, t12b, t13, &
                       t13b, t14, t14b, t15, t15a, t15b, t15p, t16, t16a, t16b,&
                       t17, t17a, t17b, t18, t18a, t18b, t20, t20b, t21, t21b, &
                       t22, t22a, t22b, t22c, t22p, t23, t23b, t24, t24b, t25, &
                       t25b, t26, t26b, t27, t27b, t28, t28b, t29, t29b, t30,  &
                       t30a, t30b, t31, t31b, t32, t32b, t33, t33b, t34, t34b, &
                       t35, t35b, t36, t37, t38, t39, t40, t40a, t40b, t41,    &
                       t41b, t50, t50b, t51, t51b, t60, t60b, t61, t61b, t70,  &
                       t70b, t80, t80a, t80b
   double precision :: y2, y2b, y3, y3b, y4, y4b, y6, y6b, y7, y7b, y9, y9b,   &
                       y12, y12b, y13, y13b, y14, y14b, y15, y15b, y16, y16b,  &
                       y17, y17b, y18, y18b, y19, y19b, y20, y21, y22, y23,    &
                       y24, y25, y26, y27, y28, y29, y30, y31
   double precision :: Z2, Z2a, Zc, Zij

   ! Counters for loops.
   integer :: ifunct, jfunct, kfunct, nci, ncj, nck, lk, lij, l1, l2, l3, l4, &
              l5, l6, l7
   allocate(Jx(M), Ll(3), Q(3), W(3))

   ns  = nshell(0)  ; np  = nshell(1) ; nd  = nshell(2)
   nsd = nshelld(0) ; npd = nshelld(1); ndd = nshelld(2)

   SQ3 = 1.0D0
   if (NORM) SQ3 = sqrt(3.D0)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) /2
   enddo
   do ifunct = 1,  M
      Jx(ifunct) = (2*M - ifunct) * (ifunct -1) / 2
   enddo

   call int2G(frc, natom, ntatom, r, d)
   call g2g_timer_sum_start('Exchange-correlation gradients')
   if (calc_energy) then
      call g2g_solve_groups(2, t1, frc)
   else
      call g2g_solve_groups(3, t1, frc)
   endif
   call g2g_timer_sum_stop('Exchange-correlation gradients')

   ! Checks of calculation of Coulomb gradients is done in GPU.
   call g2g_timer_sum_start('Coulomb gradients')
   call aint_query_gpu_level(igpu)
   if (igpu.gt.2) then
      call aint_coulomb_forces(frc)
      call g2g_timer_sum_stop('Coulomb gradients')
      return
   endif

   ! (ss|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      rho_ind = ifunct + Jx(jfunct)

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         tj   = a(jfunct,ncj) / Zij
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = 1, nsd
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))

               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  uf   = ad(kfunct,nck) * ti * dpc
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)

                  te = ccoef * rho_mat(rho_ind) * af(kfunct)
                  ty = 2.D0  * te
                  do l1 = 1, 3
                     pss = (Q(l1) - r(Nuc(ifunct),l1)) * sss + &
                           (W(l1) - Q(l1)            ) * ss1s
                     sps = pss + (r(Nuc(ifunct),l1) - r(Nuc(jfunct),l1)) * sss
                     ssp = (W(l1) - r(Nucd(kfunct),l1)) * ss1s
                     frc(Nuc(ifunct),l1)  = frc(Nuc(ifunct),l1)  + ty * &
                                              a(ifunct,nci) * pss
                     frc(Nuc(jfunct),l1)  = frc(Nuc(jfunct),l1)  + ty * &
                                              a(jfunct,ncj) * sps
                     frc(Nucd(kfunct),l1) = frc(Nucd(kfunct),l1) + ty * &
                                              ad(kfunct,nck) * ssp
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo
   !print*, "ss|s", frc

   ! (ps|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1   , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            ti   = a(ifunct,nci) / Zij
            Z2   = 2.D0 * Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = 1, nsd
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ta   = (sss - tj * ss1s) / Z2
                  tb   = ss1s / Z2a

                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s

                     rho_ind = ifunct + l1 -1 + Jx(jfunct)
                     ty = ccoef * rho_mat(rho_ind) * af(kfunct)
                     te = 2.D0  * ty
                     do l2 = 1, 3
                        t1  = Q(l2) - r(Nuc(ifunct),l2)
                        t2  = W(l2) - Q(l2)

                        dss = t1 * ps + t2 * p1s
                        psp = (W(l2) - r(Nucd(kfunct),l2)) * p1s
                        if (l1 .eq. l2) then
                           dss = dss + ta
                           psp = psp + tb
                           frc(Nuc(ifunct),l2) = frc(Nuc(ifunct),l2) - &
                                                   ty * sss
                        endif
                        pps = dss + (r(Nuc(ifunct),l2) - r(Nuc(jfunct),l2)) * ps

                        frc(Nuc(ifunct),l2)  = frc(Nuc(ifunct),l2) + &
                                                 te * a(ifunct,nci) * dss
                        frc(Nuc(jfunct),l2)  = frc(Nuc(jfunct),l2) + &
                                                 te * a(jfunct,ncj) * pps
                        frc(Nucd(kfunct),l2) = frc(Nucd(kfunct),l2)+ &
                                                 te * ad(kfunct,nck) * psp
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo
   !print*, "ps|s", frc

   ! (pp|s) and grad ients
   do ifunct = ns+1, ns+np , 3
   do jfunct = ns+1, ifunct, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = 1, nsd
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck)*sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)

                  t3 = (sss  - tj * ss1s) / Z2
                  t4 = (ss1s - tj * ss2s) / Z2
                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     t5  = (ps - tj * p1s) / Z2
                     t8  = p1s / Z2a

                     lij = 3
                     if (ifunct .eq. jfunct) lij = l1
                     do l2 = 1, lij
                        t2   = W(l2) - Q(l2)
                        ta   = Q(l2)  - r(Nuc(jfunct),l2)
                        sps  = ta * sss  + t2 * ss1s
                        sp1s = ta * ss1s + t2 * ss2s
                        pps  = ta * ps   + t2 * p1s
                        pp1s = ta * p1s  + t2 * p2s
                        t6   = (sps - tj * sp1s) / Z2
                        t7   = sp1s / Z2a
                        if (l1 .eq. l2) then
                           pps  = pps  + t3
                           pp1s = pp1s + t4
                        endif

                        rho_ind = ifunct + l1 -1 + Jx(jfunct + l2 -1)
                        ty = ccoef * rho_mat(rho_ind) * af(kfunct)
                        te = ty    * 2.0D0

                        do l3 = 1, 3
                           dps = (Q(l3) - r(Nuc(ifunct),l3))  * pps + &
                                 (W(l3) - Q(l3)            )  * pp1s
                           ppp = (W(l3) - r(Nucd(kfunct),l3)) * pp1s

                           if (l1 .eq. l3) then
                              dps = dps + t6
                              ppp = ppp + t7
                              frc(Nuc(ifunct),l3) = frc(Nuc(ifunct),l3) - &
                                                      ty * sps
                           endif
                           if (l2 .eq. l3) then
                              dps = dps + t5
                              ppp = ppp + t8
                              frc(Nuc(jfunct),l3) = frc(Nuc(jfunct),l3) - &
                                                      ty * ps
                           endif
                           pds = (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) * &
                                 pps + dps

                           frc(Nuc(ifunct),l3)  = frc(Nuc(ifunct),l3)  + &
                                                    te * a(ifunct,nci)  * dps
                           frc(Nuc(jfunct),l3)  = frc(Nuc(jfunct),l3)  + &
                                                    te * a(jfunct,ncj)  * pds
                           frc(Nucd(kfunct),l3) = frc(Nucd(kfunct),l3) + &
                                                    te * ad(kfunct,nck) * ppp
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
   !print*, "pp|s", frc

   ! (ds|s) and gradients
   do ifunct = ns+np+1, M , 6
   do jfunct = 1      , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = 1,nsd
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  ti    = Zij  / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ta   = (sss  - tj * ss1s) / Z2
                  tb   = (ss1s - tj * ss2s) / Z2
                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     t12 = (ps - tj * p1s) / Z2
                     t13 = p1s / Z2a
                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        dss  = t1 * ps   + t2 * p1s
                        ds1s = t1 * p1s  + t2 * p2s
                        pj0s = t1 * sss  + t2 * ss1s
                        pj1s = t1 * ss1s + t2 * ss2s
                        t10  = (pj0s - tj * pj1s) / Z2
                        t11  = pj1s / Z2a
                        f1   = 1.D0
                        if (l1 .eq. l2) then
                           dss  = dss  + ta
                           ds1s = ds1s + tb
                           f1   = SQ3
                        endif
                        rho_ind = ifunct + Ll(l1) + l2-1 + Jx(jfunct)
                        te = rho_mat(rho_ind) * af(kfunct) * ccoef / f1
                        ty = 2.D0 * te
                        do l3 = 1, 3
                           dps = (Q(l3) - r(Nuc(jfunct),l3))  * dss + &
                                 (W(l3) - Q(l3)            )  * ds1s
                           dsp = (W(l3) - r(Nucd(kfunct),l3)) * ds1s
                           if (l1 .eq. l3) then
                              dps = dps + t10
                              dsp = dsp + t11
                              frc(Nuc(ifunct),l3) = frc(Nuc(ifunct),l3) - &
                                                      te * pj0s
                           endif
                           if (l2 .eq. l3) then
                              dps = dps + t12
                              dsp = dsp + t13
                              frc(Nuc(ifunct),l3) = frc(Nuc(ifunct),l3) - &
                                                      te * ps
                           endif
                           fss = dps - &
                                 (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) * dss

                           frc(Nuc(ifunct),l3)  = frc(Nuc(ifunct),l3)  + &
                                                    a(ifunct,nci)  * ty * fss
                           frc(Nuc(jfunct),l3)  = frc(Nuc(jfunct),l3)  + &
                                                    a(jfunct,ncj)  * ty * dps
                           frc(Nucd(kfunct),l3) = frc(Nucd(kfunct),l3) + &
                                                    ad(kfunct,nck) * ty * dsp
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
   !print*, "ds|s", frc

   ! (dp|s)
   do ifunct = ns+np+1, M    , 6
   do jfunct = ns+1   , ns+np, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = 1,nsd
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  ti    = Zij  / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf = tj * Zij * dpc
                  t2 = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  t3   = (sss  - tj * ss1s) / Z2
                  t4   = (ss1s - tj * ss2s) / Z2
                  t4b  = (ss2s - tj * ss3s) / Z2

                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     p3s = t1 * ss3s + t2 * ss4s
                     t5  = (ps  - tj * p1s) / Z2
                     t5b = (p1s - tj * p2s) / Z2

                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        pjs  = t1 * sss  + t2 * ss1s
                        pj1s = t1 * ss1s + t2 * ss2s
                        pj2s = t1 * ss2s + t2 * ss3s
                        ds   = t1 * ps   + t2 * p1s
                        d1s  = t1 * p1s  + t2 * p2s
                        d2s  = t1 * p2s  + t2 * p3s
                        t6   = (pjs  - tj * pj1s) / Z2
                        t6b  = (pj1s - tj * pj2s) / Z2
                        f1   = 1.0D0

                        if (l1 .eq. l2) then
                           f1  = SQ3
                           ds  = ds  + t3
                           d1s = d1s + t4
                           d2s = d2s + t4b
                        endif
                        t14 = (ds - tj * d1s) / Z2
                        t15 = d1s / Z2a

                        do l3 = 1, 3
                           t1 = Q(l3) - r(Nuc(jfunct),l3)
                           t2 = W(l3) - Q(l3)
                           dps  = t1 * ds   + t2 * d1s
                           dp1s = t1 * d1s  + t2 * d2s
                           pi0p = t1 * ps   + t2 * p1s
                           pi1p = t1 * p1s  + t2 * p2s
                           pj0p = t1 * pjs  + t2 * pj1s
                           pj1p = t1 * pj1s + t2 * pj2s

                           if (l1 .eq. l3) then
                              dps  = dps  + t6
                              dp1s = dp1s + t6b
                              pi0p = pi0p + t3
                              pi1p = pi1p + t4
                           endif
                           if (l2 .eq. l3) then
                              dps  = dps  + t5
                              dp1s = dp1s + t5b
                              pj0p = pj0p + t3
                              pj1p = pj1p + t4
                           endif
                           t10 = (pj0p - tj * pj1p) / Z2
                           t11 = pj1p / Z2a
                           t12 = (pi0p - tj * pi1p) / Z2
                           t13 = pi1p / Z2a

                           rho_ind = ifunct + Ll(l1) + l2-1 + Jx(jfunct + l3-1)
                           te = af(kfunct) * rho_mat(rho_ind) * ccoef/f1
                           ty = 2.D0 * te
                           do l4 = 1, 3
                              dds = (Q(l4) - r(Nuc(jfunct),l4))  * dps + &
                                    (W(l4) - Q(l4)            )  * dp1s
                              dpp = (W(l4) - r(Nucd(kfunct),l4)) * dp1s

                              if (l1 .eq. l4) then
                                 dds = dds + t10
                                 dpp = dpp + t11
                                 frc(Nuc(ifunct),l4) = frc(Nuc(ifunct),l4) &
                                                         - te * pj0p
                              endif
                              if (l2 .eq. l4) then
                                 dds = dds + t12
                                 dpp = dpp + t13
                                 frc(Nuc(ifunct),l4) = frc(Nuc(ifunct),l4) &
                                                         - te * pi0p
                              endif
                              if (l3 .eq. l4) then
                                 dds = dds + t14
                                 dpp = dpp + t15
                                 frc(Nuc(jfunct),l4) = frc(Nuc(jfunct),l4) &
                                                         - te * ds
                              endif
                              fps = dds - &
                                    (r(Nuc(ifunct),l4) - r(Nuc(jfunct),l4)) *dps
                              frc(Nuc(ifunct),l4)  = frc(Nuc(ifunct),l4)  +&
                                                       a(ifunct,nci)  * ty * fps
                              frc(Nuc(jfunct),l4)  = frc(Nuc(jfunct),l4)  +&
                                                       a(jfunct,ncj)  * ty * dds
                              frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) +&
                                                       ad(kfunct,nck) * ty * dpp
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
   !print*, "dp|s", frc

   ! (dd|s)
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+np+1, ifunct, 6
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = 1,nsd
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  ti    = Zij  / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  ss5s = t2 * FUNCT(5,uf)
                  t3   = (sss  - tj * ss1s) / Z2
                  t4   = (ss1s - tj * ss2s) / Z2
                  t5   = (ss2s - tj * ss3s) / Z2
                  t5b  = (ss3s - tj * ss4s) / Z2

                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     p3s = t1 * ss3s + t2 * ss4s
                     p4s = t1 * ss4s + t2 * ss5s
                     t6  = (ps  - tj * p1s) / Z2
                     t7  = (p1s - tj * p2s) / Z2
                     t7b = (p2s - tj * p3s) / Z2

                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        pjs  = t1 * sss  + t2 * ss1s
                        pj1s = t1 * ss1s + t2 * ss2s
                        pj2s = t1 * ss2s + t2 * ss3s
                        pj3s = t1 * ss3s + t2 * ss4s
                        ds   = t1 * ps   + t2 * p1s
                        d1s  = t1 * p1s  + t2 * p2s
                        d2s  = t1 * p2s  + t2 * p3s
                        d3s  = t1 * p3s  + t2 * p4s
                        t8   = (pjs  - tj * pj1s) / Z2
                        t9   = (pj1s - tj * pj2s) / Z2
                        t9b  = (pj2s - tj * pj3s) / Z2

                        f1   = 1.0D0
                        if (l1 .eq. l2) then
                           ds  = ds  + t3
                           d1s = d1s + t4
                           d2s = d2s + t5
                           d3s = d3s + t5b
                           f1  = SQ3
                        endif
                        t12  = (ds  - tj * d1s) / Z2
                        t12b = (d1s - tj * d2s) / Z2

                        lij = 3
                        if (ifunct .eq. jfunct) lij = l1
                        do l3 = 1, lij
                           t1 = Q(l3) - r(Nuc(jfunct),l3)
                           t2 = W(l3) - Q(l3)
                           s0pk = t1 * sss + t2 * ss1s
                           s1pk = t1 * ss1s + t2 * ss2s
                           s2pk = t1 * ss2s + t2 * ss3s
                           pip  = t1 * ps + t2 * p1s
                           pi1p = t1 * p1s + t2 * p2s
                           pi2p = t1 * p2s + t2 * p3s
                           pjp  = t1 * pjs + t2 * pj1s
                           pj1p = t1 * pj1s + t2 * pj2s
                           pj2p = t1 * pj2s + t2 * pj3s
                           dp   = t1 * ds + t2 * d1s
                           d1p  = t1 * d1s + t2 * d2s
                           d2p  = t1 * d2s + t2 * d3s

                           t13  = (s0pk - tj * s1pk) / Z2
                           t13b = (s1pk - tj * s2pk) / Z2
                           if (l1 .eq. l3) then
                              pip  = pip  + t3
                              pi1p = pi1p + t4
                              pi2p = pi2p + t5
                              dp   = dp   + t8
                              d1p  = d1p  + t9
                              d2p  = d2p  + t9b
                           endif
                           if (l2 .eq. l3) then
                              pjp  = pjp  + t3
                              pj1p = pj1p + t4
                              pj2p = pj2p + t5
                              dp   = dp   + t6
                              d1p  = d1p  + t7
                              d2p  = d2p  + t7b
                           endif
                           t10  = (pjp  - tj * pj1p) / Z2
                           t10b = (pj1p - tj * pj2p) / Z2
                           t11  = (pip  - tj * pi1p) / Z2
                           t11b = (pi1p - tj * pi2p) / Z2
                           t26  = (dp   - tj * d1p)  / Z2
                           t27  = d1p / Z2a

                           lk = l3
                           if (ifunct .eq. jfunct) lk = min(l3,Ll(l1)-Ll(l3)+l2)
                           do l4 = 1, lk
                              t1   = Q(l4) - r(Nuc(jfunct),l4)
                              t2   = W(l4) - Q(l4)
                              dds  = t1 * dp + t2 * d1p
                              dd1s = t1 * d1p + t2 * d2p
                              pi0d = t1 * pip + t2 * pi1p
                              pi1d = t1 * pi1p + t2 * pi2p
                              pj0d = t1 * pjp + t2 * pj1p
                              pj1d = t1 * pj1p + t2 * pj2p
                              d0pl = t1 * ds + t2 * d1s
                              d1pl = t1 * d1s + t2 * d2s

                              f2 = 1.0D0
                              if (l1 .eq. l4) then
                                 dds  = dds  + t10
                                 dd1s = dd1s + t10b
                                 pi0d = pi0d + t13
                                 pi1d = pi1d + t13b
                                 d0pl = d0pl + t8
                                 d1pl = d1pl + t9
                              endif
                              if (l2 .eq. l4) then
                                 dds  = dds  + t11
                                 dd1s = dd1s + t11b
                                 pj0d = pj0d + t13
                                 pj1d = pj1d + t13b
                                 d0pl = d0pl + t6
                                 d1pl = d1pl + t7
                              endif
                              if (l3 .eq. l4) then
                                 dds  = dds + t12
                                 dd1s = dd1s + t12b
                                 pi0d = pi0d + t6
                                 pi1d = pi1d + t7
                                 pj0d = pj0d + t8
                                 pj1d = pj1d + t9
                                 f2   = SQ3
                              endif
                              t20 = (pj0d - tj * pj1d) / Z2
                              t21 = pj1d / Z2a
                              t22 = (pi0d - tj * pi1d) / Z2
                              t23 = pi1d / Z2a
                              t24 = (d0pl - tj * d1pl) / Z2
                              t25 = d1pl / Z2a

                              rho_ind = ifunct + Ll(l1) + l2 -1 + &
                                         Jx(jfunct + Ll(l3) + l4 -1)
                              te = af(kfunct) * rho_mat(rho_ind) * ccoef/(f1*f2)
                              ty = 2.D0 * te
                              do l5 = 1, 3
                                 fds = (Q(l5) - r(Nuc(ifunct),l5) ) * dds + &
                                       (W(l5) - Q(l5)             ) * dd1s
                                 ddp = (W(l5) - r(Nucd(kfunct),l5)) * dd1s
                                 if (l1 .eq. l5) then
                                    fds = fds + t20
                                    ddp = ddp + t21
                                    frc(Nuc(ifunct),l5) = frc(Nuc(ifunct),l5) -&
                                                          te * pj0d
                                 endif
                                 if (l2 .eq. l5) then
                                    fds = fds + t22
                                    ddp = ddp + t23
                                    frc(Nuc(ifunct),l5) = frc(Nuc(ifunct),l5) -&
                                                          te * pi0d
                                 endif
                                 if (l3 .eq. l5) then
                                    fds = fds + t24
                                    ddp = ddp + t25
                                    frc(Nuc(jfunct),l5) = frc(Nuc(jfunct),l5) -&
                                                          te * d0pl
                                 endif
                                 if (l4 .eq. l5) then
                                    fds = fds + t26
                                    ddp = ddp + t27
                                    frc(Nuc(jfunct),l5) = frc(Nuc(jfunct),l5) -&
                                                          te * dp
                                 endif
                                 dfs = fds + dds * &
                                       (r(Nuc(ifunct),l5) - r(Nuc(jfunct),l5))
                                 frc(Nuc(ifunct),l5)  = frc(Nuc(ifunct),l5)  +&
                                                        ty * a(ifunct,nci)  *fds
                                 frc(Nuc(jfunct),l5)  = frc(Nuc(jfunct),l5)  +&
                                                        ty * a(jfunct,ncj)  *dfs
                                 frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5) +&
                                                        ty * ad(kfunct,nck) *ddp
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
   !print*, "dd|s", frc

   ! (ss|p)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+1, nsd+npd, 3
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = ti * ad(kfunct,nck) * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ta   = (sss - ti * ss1s) / Zc
                  tb   = ss1s / Z2a
                  do l1 = 1, 3
                     t1   = W(l1) - r(Nucd(kfunct),l1)
                     ssp  = t1 * ss1s
                     ss1p = t1 * ss2s

                     af_ind   = kfunct + l1-1
                     rho_ind = ifunct + Jx(jfunct)
                     te = ccoef * rho_mat(rho_ind) * af(af_ind)
                     ty = 2.D0*te
                     do l2 = 1, 3
                        ssd = (W(l2) - r(Nucd(kfunct),l2)) * ss1p
                        psp = (Q(l2) - r(Nuc(ifunct),l2) ) * ssp + &
                              (W(l2) - Q(l2))              * ss1p
                        if (l1 .eq. l2) then
                           ssd = ssd + ta
                           psp = psp + tb
                           frc(Nucd(kfunct),l2) = frc(Nucd(kfunct),l2) - te *sss
                        endif
                        spp = psp + (r(Nuc(ifunct),l2) - r(Nuc(jfunct),l2)) *ssp

                        frc(Nuc(ifunct),l2)  = frc(Nuc(ifunct),l2) +&
                                               a(ifunct,nci)  * ty * psp
                        frc(Nuc(jfunct),l2)  = frc(Nuc(jfunct),l2) +&
                                               a(jfunct,ncj)  * ty * spp
                        frc(Nucd(kfunct),l2) = frc(Nucd(kfunct),l2) +&
                                               ad(kfunct,nck) * ty * ssd
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo
   !print*, "ss|p", frc

   ! (ps|p)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1   , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+1,nsd+npd,3
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  Z2a   = 2.0D0 * t0
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = ad(kfunct,nck) * dpc * ti
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  t3   = ss2s / Z2a
                  t3a  = ss1s / Z2a
                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     t5  = (ps - ti * p1s) / Zc
                     t9  = p1s / Z2a

                     do l2 = 1, 3
                        t1     = W(l2) - r(Nucd(kfunct),l2)
                        ss0pj  = t1 * ss1s
                        sspj   = t1 * ss2s
                        pispj  = t1 * p2s
                        pi0spj = t1 * p1s
                        t4     = sspj / Z2a
                        t10    = (ss0pj - tj * sspj) / Z2
                        if (l1 .eq. l2) then
                           pi0spj = pi0spj + t3a
                           pispj  = pispj  + t3
                        endif

                        af_ind=kfunct+l2-1
                        af_ind   = kfunct + l2 -1
                        rho_ind = ifunct + l1 -1 + Jx(jfunct)
                        te = ccoef * rho_mat(rho_ind) * af(af_ind)
                        ty = 2.0D0 * te
                        do l3 = 1, 3
                           psd = (W(l3) - r(Nucd(kfunct),l3)) * pispj
                           dsp = (Q(l3) - r(Nuc(ifunct),l3))  * pi0spj +&
                                 (W(l3) - Q(l3))              * pispj
                           if (l1 .eq. l3) then
                              psd = psd + t4
                              dsp = dsp + t10
                              frc(Nuc(ifunct),l3) = frc(Nuc(ifunct),l3) - &
                                                    te * ss0pj
                           endif
                           if (l2 .eq. l3) then
                              psd = psd + t5
                              dsp = dsp + t9
                              frc(Nucd(kfunct),l3) = frc(Nucd(kfunct),l3) - &
                                                     te * ps
                           endif
                           ppp = dsp + (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) &
                                 * pi0spj
                           frc(Nuc(ifunct),l3)  = frc(Nuc(ifunct),l3)  + &
                                                  a(ifunct,nci)  * ty * dsp
                           frc(Nuc(jfunct),l3)  = frc(Nuc(jfunct),l3)  + &
                                                  a(jfunct,ncj)  * ty * ppp
                           frc(Nucd(kfunct),l3) = frc(Nucd(kfunct),l3) + &
                                                  ad(kfunct,nck) * ty * psd
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
   !print*, "ps|p", frc

   ! (pp|p)
   do ifunct = ns+1, ns+np , 3
   do jfunct = ns+1, ifunct, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij
            do kfunct = nsd+1,nsd+npd,3
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  t3   = (sss  - tj * ss1s) / Z2
                  t4   = (ss1s - tj * ss2s) / Z2
                  t5   = (ss2s - tj * ss3s) / Z2
                  t6a  = ss1s / Z2a
                  t6   = ss2s / Z2a

                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     p3s = t1 * ss3s + t2 * ss4s
                     t8a = p1s / Z2a
                     t8  = p2s / Z2a

                     lij = 3
                     if (ifunct .eq. jfunct) lij = l1
                     do l2 = 1, lij
                        t1    = Q(l2) - r(Nuc(jfunct),l2)
                        t2    = W(l2) - Q(l2)
                        pijs  = t1 * ps   + t2 * p1s
                        pij1s = t1 * p1s  + t2 * p2s
                        pij2s = t1 * p2s  + t2 * p3s
                        spjs  = t1 * ss1s + t2 * ss2s
                        sp2js = t1 * ss2s + t2 * ss3s
                        t7    = sp2js / Z2a
                        t7a   = spjs  / Z2a
                        if (l1 .eq. l2) then
                           pijs  = pijs  + t3
                           pij1s = pij1s + t4
                           pij2s = pij2s + t5
                        endif
                        t11 = (pijs - ti * pij1s) / Zc
                        t12 = pij1s / Z2a

                        do l3 = 1, 3
                           t1    = W(l3) - r(Nucd(kfunct),l3)
                           ppp   = t1 * pij1s
                           spp   = t1 * spjs
                           psp   = t1 * p1s
                           pp1p  = t1 * pij2s
                           spjpk = t1 * sp2js
                           pispk = t1 * p2s
                           if (l1 .eq. l3) then
                              psp   = psp   + t6a
                              ppp   = ppp   + t7a
                              pp1p  = pp1p  + t7
                              pispk = pispk + t6
                           endif
                           if (l2 .eq. l3) then
                              spp   = spp   + t6a
                              ppp   = ppp   + t8a
                              pp1p  = pp1p  + t8
                              spjpk = spjpk + t6
                           endif
                           t9  = spjpk / Z2a
                           t10 = pispk / Z2a
                           t14 = (psp - tj * pispk) / Z2
                           t15 = (spp - tj * spjpk) / Z2

                           af_ind   = kfunct + l3 -1
                           rho_ind = ifunct + l1 -1 + Jx(jfunct + l2 -1)
                           te = ccoef * rho_mat(rho_ind) * af(af_ind)
                           ty = 2.0D0 * te
                           do l4 = 1, 3
                              ppd = (W(l4) - r(Nucd(kfunct),l4)) * pp1p
                              dpp = (Q(l4) - r(Nuc(ifunct),l4) ) * ppp + &
                                    (W(l4) - Q(l4))              * pp1p
                              if (l1 .eq. l4) then
                                 ppd = ppd + t9
                                 dpp = dpp + t15
                                 frc(Nuc(ifunct),l4) = frc(Nuc(ifunct),l4) - &
                                                       te * spp
                              endif
                              if (l2 .eq. l4) then
                                 ppd = ppd + t10
                                 dpp = dpp + t14
                                 frc(Nuc(jfunct),l4) = frc(Nuc(jfunct),l4) - &
                                                       te * psp
                              endif
                              if (l3 .eq. l4) then
                                 ppd = ppd + t11
                                 dpp = dpp + t12
                                 frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) - &
                                                        te * pijs
                              endif
                              pdp = dpp + &
                                    (r(Nuc(ifunct),l4) - r(Nuc(jfunct),l4)) *ppp

                              frc(Nuc(ifunct),l4)  = frc(Nuc(ifunct),l4)  +&
                                                     a(ifunct,nci)  * ty * dpp
                              frc(Nuc(jfunct),l4)  = frc(Nuc(jfunct),l4)  +&
                                                     a(jfunct,ncj)  * ty * pdp
                              frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) +&
                                                     ad(kfunct,nck) * ty * ppd
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
   !print*, "pp|p", frc

   ! (ds|p)
   do ifunct = ns+np+1, M , 6
   do jfunct = 1      , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+1,nsd+npd,3
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  t3a = (sss  - tj * ss1s) / Z2
                  t3  = (ss1s - tj * ss2s) / Z2
                  t3b = (ss2s - tj * ss3s) / Z2
                  t7  = ss1s / Z2a
                  t7b = ss2s / Z2a

                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     p3s = t1 * ss3s + t2 * ss4s
                     t5  = p1s / Z2a
                     t5b = p2s / Z2a

                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        pj1s = t1 * ss1s + t2 * ss2s
                        pj2s = t1 * ss2s + t2 * ss3s
                        d0s  = t1 * ps + t2 * p1s
                        ds   = t1 * p1s + t2 * p2s
                        d2s  = t1 * p2s + t2 * p3s
                        t6   = pj1s / Z2a
                        t6b  = pj2s / Z2a
                        f1   = 1.0D0

                        if (l1 .eq. l2) then
                           d2s = d2s + t3b
                           ds  = ds  + t3
                           d0s = d0s + t3a
                           f1  = SQ3
                        endif
                        t14 = ds / Z2a
                        t15 = (d0s - ti * ds) / Zc

                        do l3 = 1, 3
                           t1   = W(l3) - r(Nucd(kfunct),l3)
                           dsp  = t1 * ds
                           ds1p = t1 * d2s
                           pi0p = t1 * p1s
                           pi1p = t1 * p2s
                           pj0p = t1 * pj1s
                           pj1p = t1 * pj2s

                           if (l1 .eq. l3) then
                              dsp  = dsp + t6
                              ds1p = ds1p + t6b
                              pi0p = pi0p + t7
                              pi1p = pi1p + t7b
                           endif
                           if (l2 .eq. l3) then
                              dsp = dsp + t5
                              ds1p = ds1p + t5b
                              pj0p = pj0p + t7
                              pj1p = pj1p + t7b
                           endif
                           t10 = (pj0p - tj * pj1p) / Z2
                           t11 = pj1p / Z2a
                           t12 = (pi0p - tj * pi1p) / Z2
                           t13 = pi1p / Z2a

                           af_ind   = kfunct + l3 -1
                           rho_ind = ifunct + Ll(l1) + l2 -1 + Jx(jfunct)
                           te = af(af_ind) * rho_mat(rho_ind) * ccoef / f1
                           ty = 2.0D0 * te

                           do l4 = 1, 3
                              dpp = (Q(l4) - r(Nuc(jfunct),l4))  * dsp + &
                                    (W(l4) - Q(l4))              * ds1p
                              dsd = (W(l4) - r(Nucd(kfunct),l4)) * ds1p

                              if (l1 .eq. l4) then
                                 dpp = dpp + t10
                                 dsd = dsd + t11
                                 frc(Nuc(ifunct),l4) = frc(Nuc(ifunct),l4) - &
                                                       te * pj0p
                              endif
                              if (l2 .eq. l4) then
                                 dpp = dpp + t12
                                 dsd = dsd + t13
                                 frc(Nuc(ifunct),l4) = frc(Nuc(ifunct),l4) - &
                                                       te * pi0p
                              endif
                              if (l3 .eq. l4) then
                                 dpp = dpp + t14
                                 dsd = dsd + t15
                                 frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) - &
                                                        te * d0s
                              endif
                              fsp = dpp - &
                                    (r(Nuc(ifunct),l4) - r(Nuc(jfunct),l4)) *dsp

                              frc(Nuc(ifunct),l4)  = frc(Nuc(ifunct),l4)  +&
                                                     ty * a(ifunct,nci)  * fsp
                              frc(Nuc(jfunct),l4)  = frc(Nuc(jfunct),l4)  +&
                                                     ty * a(jfunct,ncj)  * dpp
                              frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) +&
                                                     ty * ad(kfunct,nck) * dsd
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
   !print*, "ds|p", frc

   ! (dp|p)
   do ifunct = ns+np+1, M    , 6
   do jfunct = ns+1   , ns+np, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+1,nsd+npd,3
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  ss5s = t2 * FUNCT(5,uf)
                  t3a  = (sss  - tj * ss1s) / Z2
                  t3   = (ss1s - tj * ss2s) / Z2
                  t4   = (ss2s - tj * ss3s) / Z2
                  t4b  = (ss3s - tj * ss4s) / Z2

                  do l1 = 1, 3
                     t1   = Q(l1) - r(Nuc(ifunct),l1)
                     t2   = W(l1) - Q(l1)
                     ps   = t1 * sss  + t2 * ss1s
                     p1s  = t1 * ss1s + t2 * ss2s
                     p2s  = t1 * ss2s + t2 * ss3s
                     p3s  = t1 * ss3s + t2 * ss4s
                     p4s  = t1 * ss4s + t2 * ss5s
                     t5a  = (ps  - tj * p1s) / Z2
                     t5   = (p1s - tj * p2s) / Z2
                     t5b  = (p2s - tj * p3s) / Z2
                     t12  = p1s / Z2a
                     t12b = p2s / Z2a

                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        ds   = t1 * ps + t2 * p1s
                        d1s  = t1 * p1s + t2 * p2s
                        d2s  = t1 * p2s + t2 * p3s
                        d3s  = t1 * p3s + t2 * p4s
                        pjs  = t1 * sss + t2 * ss1s
                        pj1s = t1 * ss1s + t2 * ss2s
                        pj2s = t1 * ss2s + t2 * ss3s
                        pj3s = t1 * ss3s + t2 * ss4s
                        t6a  = (pjs  - tj * pj1s) / Z2
                        t6   = (pj1s - tj * pj2s) / Z2
                        t6b  = (pj2s - tj * pj3s) / Z2
                        t10  = pj1s / Z2a
                        t10b = pj2s / Z2a
                        f1   = 1.0D0

                        if (l1 .eq. l2) then
                           ds  = ds + t3a
                           d1s = d1s + t3
                           d2s = d2s + t4
                           d3s = d3s + t4b
                           f1  = SQ3
                        endif
                        t9  = d1s / Z2a
                        t9b = d2s / Z2a

                        do l3 = 1, 3
                           t1   = Q(l3) - r(Nuc(jfunct),l3)
                           t2   = W(l3) - Q(l3)
                           s1pk = t1 * ss1s + t2 * ss2s
                           s2pk = t1 * ss2s + t2 * ss3s
                           d0p  = t1 * ds + t2 * d1s
                           d1p  = t1 * d1s + t2 * d2s
                           d2p  = t1 * d2s + t2 * d3s
                           pi1p = t1 * p1s + t2 * p2s
                           pi2p = t1 * p2s + t2 * p3s
                           pj1p = t1 * pj1s + t2 * pj2s
                           pj2p = t1 * pj2s + t2 * pj3s
                           t11  = s1pk / Z2a
                           t11b = s2pk / Z2a

                           if (l1 .eq. l3) then
                              d0p  = d0p + t6a
                              d1p  = d1p + t6
                              d2p  = d2p + t6b
                              pi1p = pi1p + t3
                              pi2p = pi2p + t4
                           endif
                           if (l2 .eq. l3) then
                              d0p  = d0p + t5a
                              d1p  = d1p + t5
                              d2p  = d2p + t5b
                              pj1p = pj1p + t3
                              pj2p = pj2p + t4
                           endif
                           t7  = pi1p / Z2a
                           t7b = pi2p / Z2a
                           t8  = pj1p / Z2a
                           t8b = pj2p / Z2a
                           t26 = d1p  / Z2a
                           t27 = (d0p - ti * d1p) / Zc

                           do l4 = 1, 3
                              t1    = W(l4) - r(Nucd(kfunct),l4)
                              dpp   = t1 * d1p
                              d1pp  = t1 * d2p
                              pj0pp = t1 * pj1p
                              pj1pp = t1 * pj2p
                              pi0pp = t1 * pi1p
                              pi1pp = t1 * pi2p
                              d0pl  = t1 * d1s
                              d1pl  = t1 * d2s

                              if (l1 .eq. l4) then
                                 dpp   = dpp   + t8
                                 d1pp  = d1pp  + t8b
                                 d0pl  = d0pl  + t10
                                 d1pl  = d1pl  + t10b
                                 pi0pp = pi0pp + t11
                                 pi1pp = pi1pp + t11b
                              endif
                              if (l2 .eq. l4) then
                                 dpp   = dpp   + t7
                                 d1pp  = d1pp  + t7b
                                 d0pl  = d0pl  + t12
                                 d1pl  = d1pl  + t12b
                                 pj0pp = pj0pp + t11
                                 pj1pp = pj1pp + t11b
                              endif
                              if (l3 .eq. l4) then
                                 dpp   = dpp   + t9
                                 d1pp  = d1pp  + t9b
                                 pi0pp = pi0pp + t12
                                 pi1pp = pi1pp + t12b
                                 pj0pp = pj0pp + t10
                                 pj1pp = pj1pp + t10b
                              endif
                              t20 = (pj0pp - tj * pj1pp) / Z2
                              t21 = pj1pp / Z2a
                              t22 = (pi0pp - tj * pi1pp) / Z2
                              t23 = pi1pp / Z2a
                              t24 = (d0pl - tj * d1pl) / Z2
                              t25 = d1pl / Z2a

                              af_ind   = kfunct + l4 -1
                              rho_ind = ifunct + Ll(l1) + l2 -1 + &
                                         Jx(jfunct + l3 -1)
                              te = af(af_ind) * rho_mat(rho_ind) * ccoef / f1
                              ty = 2.0D0 * te

                              do l5 = 1, 3
                                 ddp = (Q(l5) - r(Nuc(jfunct),l5))  * dpp  +&
                                       (W(l5) - Q(l5))              * d1pp
                                 dpd = (W(l5) - r(Nucd(kfunct),l5)) * d1pp

                                 if (l1 .eq. l5) then
                                    ddp = ddp + t20
                                    dpd = dpd + t21
                                    frc(Nuc(ifunct),l5) = frc(Nuc(ifunct),l5) -&
                                                          te * pj0pp
                                 endif
                                 if (l2 .eq. l5) then
                                    ddp = ddp + t22
                                    dpd = dpd + t23
                                    frc(Nuc(ifunct),l5) = frc(Nuc(ifunct),l5) -&
                                                          te * pi0pp
                                 endif
                                 if (l3 .eq. l5) then
                                    ddp = ddp + t24
                                    dpd = dpd + t25
                                    frc(Nuc(jfunct),l5) = frc(Nuc(jfunct),l5) -&
                                                          te * d0pl
                                 endif
                                 if (l4 .eq. l5) then
                                    ddp = ddp + t26
                                    dpd = dpd + t27
                                    frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5)&
                                                         - te * d0p
                                 endif
                                 fpp = ddp - dpp * &
                                       (r(Nuc(ifunct),l5) - r(Nuc(jfunct),l5))

                                 frc(Nuc(ifunct),l5)  = frc(Nuc(ifunct),l5)  +&
                                                        a(ifunct,nci)  * ty *fpp
                                 frc(Nuc(jfunct),l5)  = frc(Nuc(jfunct),l5)  +&
                                                        a(jfunct,ncj)  * ty *ddp
                                 frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5) +&
                                                        ad(kfunct,nck) * ty *dpd
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
   !print*, "dp|p", frc

   ! (dd|p)
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+np+1, ifunct, 6
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+1,nsd+npd,3
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  ss5s = t2 * FUNCT(5,uf)
                  ss6s = t2 * FUNCT(6,uf)
                  ta   = (sss  - tj * ss1s) / Z2
                  t3   = (ss1s - tj * ss2s) / Z2
                  t4   = (ss2s - tj * ss3s) / Z2
                  t5   = (ss3s - tj * ss4s) / Z2
                  t5b  = (ss4s - tj * ss5s) / Z2

                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     p3s = t1 * ss3s + t2 * ss4s
                     p4s = t1 * ss4s + t2 * ss5s
                     p5s = t1 * ss5s + t2 * ss6s
                     t6a = (ps  - tj * p1s) / Z2
                     t6  = (p1s - tj * p2s) / Z2
                     t7  = (p2s - tj * p3s) / Z2
                     t7b = (p3s - tj * p4s) / Z2

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
                        t8a  = (pjs  - tj * pj1s) / Z2
                        t8   = (pj1s - tj * pj2s) / Z2
                        t9   = (pj2s - tj * pj3s) / Z2
                        t9b  = (pj3s - tj * pj4s) / Z2
                        f1   = 1.0D0

                        if (l1 .eq. l2) then
                           ds  = ds + ta
                           d1s = d1s + t3
                           d2s = d2s + t4
                           d3s = d3s + t5
                           d4s = d4s + t5b
                           f1  = SQ3
                        endif
                        t18a = (ds  - tj * d1s) / Z2
                        t18  = (d1s - tj * d2s) / Z2
                        t18b = (d2s - tj * d3s) / Z2
                        t35  = d1s / Z2a
                        t35b = d2s / Z2a

                        lij  = 3
                        if (ifunct .eq. jfunct) lij = l1
                        do l3 = 1, lij
                           t1    = Q(l3) - r(Nuc(jfunct),l3)
                           t2    = W(l3) - Q(l3)
                           d0pk  = t1 * ds   + t2 * d1s
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
                           s4pk  = t1 * ss4s + t2 * ss5s
                           t10   = (s1pk - tj * s2pk) / Z2
                           t10b  = (s2pk - tj * s3pk) / Z2

                           if (l1 .eq. l3) then
                              d0pk  = d0pk  + t8a
                              d1pk  = d1pk  + t8
                              d2pk  = d2pk  + t9
                              d3pk  = d3pk  + t9b
                              pipk  = pipk  + ta
                              pi1pk = pi1pk + t3
                              pi2pk = pi2pk + t4
                              pi3pk = pi3pk + t5
                           endif
                           if (l2 .eq. l3) then
                              d0pk = d0pk   + t6a
                              d1pk = d1pk   + t6
                              d2pk = d2pk   + t7
                              d3pk = d3pk   + t7b
                              pjpk = pjpk   + ta
                              pj1pk = pj1pk + t3
                              pj2pk = pj2pk + t4
                              pj3pk = pj3pk + t5
                           endif
                           t16a = (pjpk - tj * pj1pk) / Z2
                           t16  = (pj1pk - tj * pj2pk) / Z2
                           t16b = (pj2pk - tj * pj3pk) / Z2
                           t17a = (pipk - tj * pi1pk) / Z2
                           t17  = (pi1pk - tj * pi2pk) / Z2
                           t17b = (pi2pk - tj * pi3pk) / Z2
                           t30  = pi1pk / Z2a
                           t30b = pi2pk / Z2a
                           t31  = pj1pk / Z2a
                           t31b = pj2pk / Z2a

                           lk = l3
                           if (ifunct .eq. jfunct) lk = min(l3,Ll(l1)-Ll(l3)+l2)
                           do l4 = 1, lk
                              t1     = Q(l4) - r(Nuc(jfunct),l4)
                              t2     = W(l4) - Q(l4)
                              d0d    = t1 * d0pk  + t2 * d1pk
                              d1d    = t1 * d1pk  + t2 * d2pk
                              d2d    = t1 * d2pk  + t2 * d3pk
                              pi1pl  = t1 * p1s   + t2 * p2s
                              pi2pl  = t1 * p2s   + t2 * p3s
                              pj1pl  = t1 * pj1s  + t2 * pj2s
                              pj2pl  = t1 * pj2s  + t2 * pj3s
                              pjdkl  = t1 * pj1pk + t2 * pj2pk
                              pj2dkl = t1 * pj2pk + t2 * pj3pk
                              pidkl  = t1 * pi1pk + t2 * pi2pk
                              pi2dkl = t1 * pi2pk + t2 * pi3pk
                              d1pl   = t1 * d1s   + t2 * d2s
                              d2pl   = t1 * d2s   + t2 * d3s
                              s1ds   = t1 * s1pk  + t2 * s2pk
                              s2ds   = t1 * s2pk  + t2 * s3pk
                              f2     = 1.0D0

                              if (l1 .eq. l4) then
                                 d0d    = d0d    + t16a
                                 d1d    = d1d    + t16
                                 d2d    = d2d    + t16b
                                 pidkl  = pidkl  + t10
                                 pi2dkl = pi2dkl + t10b
                                 d1pl   = d1pl   + t8
                                 d2pl   = d2pl   + t9
                                 pi1pl  = pi1pl  + t3
                                 pi2pl  = pi2pl  + t4
                              endif
                              if (l2 .eq. l4) then
                                 d0d    = d0d    + t17a
                                 d1d    = d1d    + t17
                                 d2d    = d2d    + t17b
                                 pjdkl  = pjdkl  + t10
                                 pj2dkl = pj2dkl + t10b
                                 d1pl   = d1pl   + t6
                                 d2pl   = d2pl   + t7
                                 pj1pl  = pj1pl  + t3
                                 pj2pl  = pj2pl  + t4
                              endif
                              if (l3 .eq. l4) then
                                 d0d    = d0d    + t18a
                                 d1d    = d1d    + t18
                                 d2d    = d2d    + t18b
                                 pjdkl  = pjdkl  + t8
                                 pj2dkl = pj2dkl + t9
                                 pidkl  = pidkl  + t6
                                 pi2dkl = pi2dkl + t7
                                 s1ds   = s1ds   + t3
                                 s2ds   = s2ds   + t4
                                 f2     = SQ3
                              endif
                              t11  = pjdkl  / Z2a
                              t11b = pj2dkl / Z2a
                              t12  = pidkl  / Z2a
                              t12b = pi2dkl / Z2a
                              t13  = d1pl   / Z2a
                              t13b = d2pl   / Z2a
                              t14  = d1pk   / Z2a
                              t14b = d2pk   / Z2a
                              t32  = pi1pl  / Z2a
                              t32b = pi2pl  / Z2a
                              t33  = pj1pl  / Z2a
                              t33b = pj2pl  / Z2a
                              t34  = s1ds   / Z2a
                              t34b = s2ds   / Z2a
                              t28  = d1d    / Z2a
                              t29  = (d0d - ti * d1d) / Zc

                              do l5 = 1, 3
                                 t1    = W(l5) - r(Nucd(kfunct),l5)
                                 ddp   = t1 * d1d
                                 dd1p  = t1 * d2d
                                 pj0dp = t1 * pjdkl
                                 pj1dp = t1 * pj2dkl
                                 pi0dp = t1 * pidkl
                                 pi1dp = t1 * pi2dkl
                                 d0plp = t1 * d1pl
                                 d1plp = t1 * d2pl
                                 d0pkp = t1 * d1pk
                                 d1pkp = t1 * d2pk

                                 if (l1 .eq. l5) then
                                    ddp   = ddp   + t11
                                    dd1p  = dd1p  + t11b
                                    pi0dp = pi0dp + t34
                                    pi1dp = pi1dp + t34b
                                    d0plp = d0plp + t33
                                    d1plp = d1plp + t33b
                                    d0pkp = d0pkp + t31
                                    d1pkp = d1pkp + t31b
                                 endif
                                 if (l2 .eq. l5) then
                                    ddp   = ddp   + t12
                                    dd1p  = dd1p  + t12b
                                    pj0dp = pj0dp + t34
                                    pj1dp = pj1dp + t34b
                                    d0plp = d0plp + t32
                                    d1plp = d1plp + t32b
                                    d0pkp = d0pkp + t30
                                    d1pkp = d1pkp + t30b
                                 endif
                                 if (l3 .eq. l5) then
                                    ddp   = ddp   + t13
                                    dd1p  = dd1p  + t13b
                                    pj0dp = pj0dp + t33
                                    pj1dp = pj1dp + t33b
                                    pi0dp = pi0dp + t32
                                    pi1dp = pi1dp + t32b
                                    d0pkp = d0pkp + t35
                                    d1pkp = d1pkp + t35b
                                 endif
                                 if (l4 .eq. l5) then
                                    ddp   = ddp   + t14
                                    dd1p  = dd1p  + t14b
                                    pj0dp = pj0dp + t31
                                    pj1dp = pj1dp + t31b
                                    pi0dp = pi0dp + t30
                                    pi1dp = pi1dp + t30b
                                    d0plp = d0plp + t35
                                    d1plp = d1plp + t35b
                                 endif
                                 t20 = (pj0dp - tj * pj1dp) / Z2
                                 t21 = pj1dp / Z2a
                                 t22 = (pi0dp - tj * pi1dp) / Z2
                                 t23 = pi1dp / Z2a
                                 t24 = (d0plp - tj * d1plp) / Z2
                                 t25 = d1plp / Z2a
                                 t26 = (d0pkp - tj * d1pkp) / Z2
                                 t27 = d1pkp / Z2a

                                 af_ind   = kfunct + l5 -1
                                 rho_ind = ifunct + Ll(l1) + l2 -1 + &
                                            Jx(jfunct + Ll(l3) + l4 -1)
                                 te = af(af_ind) * rho_mat(rho_ind) * ccoef &
                                     / (f1 * f2)
                                 ty = 2.0D0 * te
                                 do l6 = 1, 3
                                    fdp = (Q(l6) - r(Nuc(ifunct),l6) ) * ddp + &
                                          (W(l6) - Q(l6))              * dd1p
                                    ddd = (W(l6) - r(Nucd(kfunct),l6)) * dd1p
                                    if (l1 .eq. l6) then
                                       fdp = fdp + t20
                                       ddd = ddd + t21
                                       frc(Nuc(ifunct),l6) =frc(Nuc(ifunct),l6)&
                                                           - te * pj0dp
                                    endif
                                    if (l2 .eq. l6) then
                                       fdp = fdp + t22
                                       ddd = ddd + t23
                                       frc(Nuc(ifunct),l6) =frc(Nuc(ifunct),l6)&
                                                           - te * pi0dp
                                    endif
                                    if (l3 .eq. l6) then
                                       fdp = fdp + t24
                                       ddd = ddd + t25
                                       frc(Nuc(jfunct),l6) =frc(Nuc(jfunct),l6)&
                                                           - te * d0plp
                                    endif
                                    if (l4 .eq. l6) then
                                       fdp = fdp + t26
                                       ddd = ddd + t27
                                       frc(Nuc(jfunct),l6) =frc(Nuc(jfunct),l6)&
                                                           - te * d0pkp
                                    endif
                                    if (l5 .eq. l6) then
                                       fdp = fdp + t28
                                       ddd = ddd + t29
                                       frc(Nucd(kfunct),l6)=&
                                                           frc(Nucd(kfunct),l6)&
                                                            - te * d0d
                                    endif
                                    dfp = fdp + ddp * &
                                          (r(Nuc(ifunct),l6) -r(Nuc(jfunct),l6))

                                    frc(Nuc(ifunct),l6)  = frc(Nuc(ifunct),l6) &
                                                         + a(ifunct,nci) *ty*fdp
                                    frc(Nuc(jfunct),l6)  = frc(Nuc(jfunct),l6) &
                                                         + a(jfunct,ncj) *ty*dfp
                                    frc(Nucd(kfunct),l6) = frc(Nucd(kfunct),l6)&
                                                         + ad(kfunct,nck)*ty*ddd
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
   enddo
   enddo
   !print*, "dd|p", frc

   ! (ss|d)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+npd+1, Md, 6
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2    = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = ti * ad(kfunct,nck) * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ta   = (sss  - ti * ss1s) / Zc
                  tb   = (ss1s - ti * ss2s) / Zc

                  do l1 = 1, 3
                     t1   = W(l1) - r(Nucd(kfunct),l1)
                     ss0p = t1 * ss1s
                     ss1p = t1 * ss2s
                     ss2p = t1 * ss3s
                     t13  = ss1p / Z2
                     t11  = (ss0p - ti * ss1p) / Zc

                     do l2 = 1, l1
                        t1    = W(l2) - r(Nucd(kfunct),l2)
                        ss0pj = t1 * ss1s
                        ss1pj = t1 * ss2s
                        ss0d  = t1 * ss1p
                        ss1d  = t1 * ss2p
                        t10   = (ss0pj - ti * ss1pj) / Zc
                        t12   = ss1pj / Z2
                        f1    = 1.0D0

                        if (l1 .eq. l2) then
                           ss0d = ss0d + ta
                           ss1d = ss1d + tb
                           f1   = SQ3
                        endif

                        af_ind   = kfunct + Ll(l1) + l2 -1
                        rho_ind = ifunct + Jx(jfunct)
                        te = af(af_ind) * rho_mat(rho_ind) * ccoef / f1
                        ty = 2.0D0 * te
                        do l3 = 1, 3
                           ssf = (W(l3) - r(Nucd(kfunct),l3)) * ss1d
                           psd = (Q(l3) - r(Nuc(ifunct),l3))  * ss0d + &
                                 (W(l3) - Q(l3))              * ss1d

                           if (l1 .eq. l3) then
                              ssf = ssf + t10
                              psd = psd + t12
                              frc(Nucd(kfunct),l3) = frc(Nucd(kfunct),l3) - &
                                                     te * ss0pj
                           endif
                           if (l2 .eq. l3) then
                              ssf = ssf + t11
                              psd = psd + t13
                              frc(Nucd(kfunct),l3) = frc(Nucd(kfunct),l3) - &
                                                     te * ss0p
                           endif
                           spd = psd + &
                                 (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) * ss0d

                           frc(Nuc(ifunct),l3)  = frc(Nuc(ifunct),l3)  + &
                                                  a(ifunct,nci)  * ty * psd
                           frc(Nuc(jfunct),l3)  = frc(Nuc(jfunct),l3)  + &
                                                  a(jfunct,ncj)  * ty * spd
                           frc(Nucd(kfunct),l3) = frc(Nucd(kfunct),l3) + &
                                                  ad(kfunct,nck) * ty * ssf
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
   !print*, "ss|d", frc

   ! (ps|d)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1   , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+npd+1,Md,6
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  Z2a   = 2.0D0 * t0
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = dpc * ti * ad(kfunct,nck)
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  t3a  = ss1s / Z2a
                  t3   = ss2s / Z2a
                  t3b  = ss3s / Z2a
                  t7   = (sss  - ti * ss1s) / Zc
                  t7b  = (ss1s - ti * ss2s) / Zc
                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     p3s = t1 * ss3s + t2 * ss4s
                     t5  = (ps  - ti * p1s) / Zc
                     t5b = (p1s - ti * p2s) / Zc

                     do l2 = 1, 3
                        t1     = W(l2) - r(Nucd(kfunct),l2)
                        sspj   = t1 * ss2s
                        ss2pj  = t1 * ss3s
                        pi0spj = t1 * p1s
                        pispj  = t1 * p2s
                        pi2spj = t1 * p3s

                        if (l1 .eq. l2) then
                           pi2spj = pi2spj + t3b
                           pispj  = pispj  + t3
                           pi0spj = pi0spj + t3a
                        endif
                        t4  = sspj  / Z2a
                        t4b = ss2pj / Z2a
                        t14 = pispj / Z2a
                        t15 = (pi0spj - ti * pispj) / Zc

                        do l3 = 1, l2
                           t1   = W(l3) - r(Nucd(kfunct),l3)
                           ss0d = t1 * sspj
                           ss1d = t1 * ss2pj
                           ps0d = t1 * pispj
                           ps1d = t1 * pi2spj
                           p0pk = t1 * p1s
                           p1pk = t1 * p2s
                           f1   = 1.0D0

                           if (l1 .eq. l3) then
                              ps0d = ps0d + t4
                              ps1d = ps1d + t4b
                              p0pk = p0pk + t3a
                              p1pk = p1pk + t3
                           endif
                           if (l2 .eq. l3) then
                              ss0d = ss0d + t7
                              ss1d = ss1d + t7b
                              ps0d = ps0d + t5
                              ps1d = ps1d + t5b
                              f1   = SQ3
                           endif
                           t10 = (ss0d - tj * ss1d) / Z2
                           t11 = ss1d / Z2a
                           t12 = p1pk / Z2a
                           t13 = (p0pk - ti * p1pk) / Zc

                           af_ind   = kfunct + Ll(l2) + l3 -1
                           rho_ind = ifunct + l1 -1 + Jx(jfunct)
                           te = af(af_ind) * rho_mat(rho_ind) * ccoef / f1
                           ty = 2.0D0 * te

                           do l4 = 1, 3
                              dsd = (Q(l4) - r(Nuc(ifunct),l4))  * ps0d + &
                                    (W(l4) - Q(l4))              * ps1d
                              psf = (W(l4) - r(Nucd(kfunct),l4)) * ps1d
                              if (l1 .eq. l4) then
                                 dsd = dsd + t10
                                 psf = psf + t11
                                 frc(Nuc(ifunct),l4) = frc(Nuc(ifunct),l4) &
                                                       - te * ss0d
                              endif
                              if (l2 .eq. l4) then
                                 dsd = dsd + t12
                                 psf = psf + t13
                                 frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) &
                                                        - te * p0pk
                              endif
                              if (l3 .eq. l4) then
                                 dsd = dsd + t14
                                 psf = psf + t15
                                 frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) &
                                                      - te * pi0spj
                              endif
                              ppd = dsd + ps0d * &
                                    (r(Nuc(ifunct),l4) - r(Nuc(jfunct),l4))
                              frc(Nuc(ifunct),l4)  = frc(Nuc(ifunct),l4)  + &
                                                     a(ifunct,nci)  * ty * dsd
                              frc(Nuc(jfunct),l4)  = frc(Nuc(jfunct),l4)  + &
                                                     a(jfunct,ncj)  * ty * ppd
                              frc(Nucd(kfunct),l4) = frc(Nucd(kfunct),l4) + &
                                                     ad(kfunct,nck) * ty * psf
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
   !print*, "ps|d", frc

   ! (pp|d)
   do ifunct = ns+1, ns+np  , 3
   do jfunct = ns+1, ifunct , 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+npd+1,Md,6
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  ss5s = t2 * FUNCT(5,uf)
                  t3   = (sss  - tj * ss1s) / Z2
                  t4   = (ss1s - tj * ss2s) / Z2
                  t5   = (ss2s - tj * ss3s) / Z2
                  t5b  = (ss3s - tj * ss4s) / Z2
                  t6   = ss2s / Z2a
                  t6b  = ss3s / Z2a

                  do l1 = 1, 3
                     t1   = Q(l1) - r(Nuc(ifunct),l1)
                     t2   = W(l1) - Q(l1)
                     ps   = t1 * sss  + t2 * ss1s
                     p1s  = t1 * ss1s + t2 * ss2s
                     p2s  = t1 * ss2s + t2 * ss3s
                     p3s  = t1 * ss3s + t2 * ss4s
                     p4s  = t1 * ss4s + t2 * ss5s
                     t8a  = p1s / Z2a
                     t8   = p2s / Z2a
                     t8b  = p3s / Z2a
                     t33  = (ps  - ti * p1s) / Zc
                     t33b = (p1s - ti * p2s) / Zc

                     lij = 3
                     if (ifunct .eq. jfunct) lij = l1
                     do l2 = 1, lij
                        t1    = Q(l2) - r(Nuc(jfunct),l2)
                        t2    = W(l2) - Q(l2)
                        pijs  = t1 * ps   + t2 * p1s
                        pij1s = t1 * p1s  + t2 * p2s
                        pij2s = t1 * p2s  + t2 * p3s
                        pij3s = t1 * p3s  + t2 * p4s
                        sp0js = t1 * sss  + t2 * ss1s
                        spjs  = t1 * ss1s + t2 * ss2s
                        sp2js = t1 * ss2s + t2 * ss3s
                        sp3js = t1 * ss3s + t2 * ss4s
                        t7a   = spjs  / Z2a
                        t7    = sp2js / Z2a
                        t7b   = sp3js / Z2a

                        if (l1 .eq. l2) then
                           pijs  = pijs  + t3
                           pij1s = pij1s + t4
                           pij2s = pij2s + t5
                           pij3s = pij3s + t5b
                        endif
                        t11  = (pijs  - ti * pij1s) / Zc
                        t11b = (pij1s - ti * pij2s) / Zc
                        t31  = (sp0js - ti * spjs)  / Zc
                        t31b = (spjs  - ti * sp2js) / Zc

                        do l3 = 1, 3
                           t1 = W(l3) - r(Nucd(kfunct),l3)
                           pp0p   = t1 * pij1s
                           pp1p   = t1 * pij2s
                           pp2p   = t1 * pij3s
                           spjpk  = t1 * sp2js
                           s2pjpk = t1 * sp3js
                           pispk  = t1 * p2s
                           pi2spk = t1 * p3s
                           sspk   = t1 * ss2s
                           ss2pk  = t1 * ss3s
                           t30    = sspk  / Z2a
                           t30b   = ss2pk / Z2a

                           if (l1 .eq. l3) then
                              pp0p   = pp0p   + t7a
                              pp1p   = pp1p   + t7
                              pp2p   = pp2p   + t7b
                              pispk  = pispk  + t6
                              pi2spk = pi2spk + t6b
                           endif
                           if (l2 .eq. l3) then
                              pp0p   = pp0p   + t8a
                              pp1p   = pp1p   + t8
                              pp2p   = pp2p   + t8b
                              spjpk  = spjpk  + t6
                              s2pjpk = s2pjpk + t6b
                           endif
                           t9   = spjpk  / Z2a
                           t9b  = s2pjpk / Z2a
                           t10  = pispk  / Z2a
                           t10b = pi2spk / Z2a
                           t26  = pp1p   / Z2a
                           t27  = (pp0p - ti * pp1p) / Zc

                           do l4 = 1, l3
                              t1    = W(l4) - r(Nucd(kfunct),l4)
                              pp0d  = t1 * pp1p
                              pp1d  = t1 * pp2p
                              sp0d  = t1 * spjpk
                              sp1d  = t1 * s2pjpk
                              ps0d  = t1 * pispk
                              ps1d  = t1 * pi2spk
                              pp0pl = t1 * pij1s
                              pp1pl = t1 * pij2s
                              f1    = 1.0D0

                              if (l1 .eq. l4) then
                                 pp0d  = pp0d  + t9
                                 pp1d  = pp1d  + t9b
                                 ps0d  = ps0d  + t30
                                 ps1d  = ps1d  + t30b
                                 pp0pl = pp0pl + t7a
                                 pp1pl = pp1pl + t7
                              endif
                              if (l2 .eq. l4) then
                                 pp0d  = pp0d  + t10
                                 pp1d  = pp1d  + t10b
                                 sp0d  = sp0d  + t30
                                 sp1d  = sp1d  + t30b
                                 pp0pl = pp0pl + t8a
                                 pp1pl = pp1pl + t8
                              endif
                              if (l3 .eq. l4) then
                                 pp0d = pp0d + t11
                                 pp1d = pp1d + t11b
                                 sp0d = sp0d + t31
                                 sp1d = sp1d + t31b
                                 ps0d = ps0d + t33
                                 ps1d = ps1d + t33b
                                 f1   = SQ3
                              endif
                              t20 = (sp0d  - tj * sp1d)  / Z2
                              t21 = sp1d / Z2a
                              t22 = (ps0d  - tj * ps1d)  / Z2
                              t23 = ps1d  / Z2a
                              t24 = pp1pl / Z2a
                              t25 = (pp0pl - ti * pp1pl) / Zc

                              af_ind   = kfunct + Ll(l3) + l4 -1
                              rho_ind = ifunct + l1 -1  + Jx(jfunct + l2 -1)
                              te = af(af_ind) * rho_mat(rho_ind) * ccoef / f1
                              ty = 2.0D0 * te

                              do l5 = 1, 3
                                 dpd = (Q(l5) - r(Nuc(ifunct),l5))  * pp0d + &
                                       (W(l5) - Q(l5))              * pp1d
                                 ppf = (W(l5) - r(Nucd(kfunct),l5)) * pp1d

                                 if (l1 .eq. l5) then
                                    dpd = dpd + t20
                                    ppf = ppf + t21
                                    frc(Nuc(ifunct),l5) = frc(Nuc(ifunct),l5)  &
                                                        - te * sp0d
                                 endif
                                 if (l2 .eq. l5) then
                                    dpd = dpd + t22
                                    ppf = ppf + t23
                                    frc(Nuc(jfunct),l5) = frc(Nuc(jfunct),l5)  &
                                                        - te * ps0d
                                 endif
                                 if (l3 .eq. l5) then
                                    dpd = dpd + t24
                                    ppf = ppf + t25
                                    frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5)&
                                                         - te * pp0pl
                                 endif
                                 if (l4 .eq. l5) then
                                    dpd = dpd + t26
                                    ppf = ppf + t27
                                    frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5)&
                                                         - te * pp0p
                                 endif
                                 pdd = dpd + pp0d * &
                                       (r(Nuc(ifunct),l5) - r(Nuc(jfunct),l5))

                                 frc(Nuc(ifunct),l5)  = frc(Nuc(ifunct),l5)  +&
                                                        a(ifunct,nci)  * ty *dpd
                                 frc(Nuc(jfunct),l5)  = frc(Nuc(jfunct),l5)  +&
                                                        a(jfunct,ncj)  * ty *pdd
                                 frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5) +&
                                                        ad(kfunct,nck) * ty *ppf
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
   !print*, "pp|d", frc

   ! (ds|d)
   do ifunct = ns+np+1, M , 6
   do jfunct = 1      , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+npd+1, Md, 6
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  ss5s = t2 * FUNCT(5,uf)
                  t3   = (sss  - tj * ss1s) / Z2
                  t4   = (ss1s - tj * ss2s) / Z2
                  t5   = (ss2s - tj * ss3s) / Z2
                  t5b  = (ss3s - tj * ss4s) / Z2
                  t6   = ss2s / Z2a
                  t6b  = ss3s / Z2a

                  do l1 = 1, 3
                     t1   = Q(l1) - r(Nuc(ifunct),l1)
                     t2   = W(l1) - Q(l1)
                     ps   = t1 * sss  + t2 * ss1s
                     p1s  = t1 * ss1s + t2 * ss2s
                     p2s  = t1 * ss2s + t2 * ss3s
                     p3s  = t1 * ss3s + t2 * ss4s
                     p4s  = t1 * ss4s + t2 * ss5s
                     t7a  = p1s / Z2a
                     t7   = p2s / Z2a
                     t7b  = p3s / Z2a
                     t16  = (ps  - ti * p1s) / Zc
                     t16b = (p1s - ti * p2s) / Zc

                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        pjs  = t1 * sss  + t2 * ss1s
                        pj1s = t1 * ss1s + t2 * ss2s
                        pj2s = t1 * ss2s + t2 * ss3s
                        pj3s = t1 * ss3s + t2 * ss4s
                        ds   = t1 * ps   + t2 * p1s
                        d1s  = t1 * p1s  + t2 * p2s
                        d2s  = t1 * p2s  + t2 * p3s
                        d3s  = t1 * p3s  + t2 * p4s
                        t8a  = pj1s / Z2a
                        t8   = pj2s / Z2a
                        t8b  = pj3s / Z2a
                        t15  = (pjs  - ti * pj1s) / Zc
                        t15b = (pj1s - ti * pj2s) / Zc
                        f1   = 1.0D0

                        if (l1 .eq. l2) then
                           ds  = ds + t3
                           d1s = d1s + t4
                           d2s = d2s + t5
                           d3s = d3s + t5b
                           f1  = SQ3
                        endif
                        t11  = (ds  - ti * d1s) / Zc
                        t11b = (d1s - ti * d2s) / Zc

                        do l3 = 1, 3
                           t1     = W(l3) - r(Nucd(kfunct),l3)
                           ds0p   = t1 * d1s
                           ds1p   = t1 * d2s
                           ds2p   = t1 * d3s
                           pis1pk = t1 * p2s
                           pis2pk = t1 * p3s
                           pjs1pk = t1 * pj2s
                           pjs2pk = t1 * pj3s
                           ss1pk  = t1 * ss2s
                           ss2pk  = t1 * ss3s
                           t12    = ss1pk / Z2a
                           t12b   = ss2pk / Z2a

                           if (l1 .eq. l3) then
                              ds0p   = ds0p   + t8a
                              ds1p   = ds1p   + t8
                              ds2p   = ds2p   + t8b
                              pis1pk = pis1pk + t6
                              pis2pk = pis2pk + t6b
                           endif
                           if (l2 .eq. l3) then
                              ds0p   = ds0p   + t7a
                              ds1p   = ds1p   + t7
                              ds2p   = ds2p   + t7b
                              pjs1pk = pjs1pk + t6
                              pjs2pk = pjs2pk + t6b
                           endif
                           t9   = pjs1pk / Z2a
                           t9b  = pjs2pk / Z2a
                           t10  = pis1pk / Z2a
                           t10b = pis2pk / Z2a
                           t26  = ds1p   / Z2a
                           t27  = (ds0p - ti * ds1p) / Zc

                           do l4 = 1, l3
                              t1    = W(l4) - r(Nucd(kfunct),l4)
                              dsd   = t1 * ds1p
                              ds1d  = t1 * ds2p
                              pj0sd = t1 * pjs1pk
                              pj1sd = t1 * pjs2pk
                              pi0sd = t1 * pis1pk
                              pi1sd = t1 * pis2pk
                              d0pl  = t1 * d1s
                              d1pl  = t1 * d2s
                              f2    = 1.0D0

                              if (l1 .eq. l4) then
                                 dsd   = dsd   + t9
                                 ds1d  = ds1d  + t9b
                                 pi0sd = pi0sd + t12
                                 pi1sd = pi1sd + t12b
                                 d0pl  = d0pl  + t8a
                                 d1pl  = d1pl  + t8
                              endif
                              if (l2 .eq. l4) then
                                 dsd   = dsd   + t10
                                 ds1d  = ds1d  + t10b
                                 pj0sd = pj0sd + t12
                                 pj1sd = pj1sd + t12b
                                 d0pl  = d0pl  + t7a
                                 d1pl  = d1pl  + t7
                              endif
                              if (l3 .eq. l4) then
                                 dsd   = dsd   + t11
                                 ds1d  = ds1d  + t11b
                                 pj0sd = pj0sd + t15
                                 pj1sd = pj1sd + t15b
                                 pi0sd = pi0sd + t16
                                 pi1sd = pi1sd + t16b
                                 f2    = SQ3
                              endif
                              t20 = (pj0sd - tj * pj1sd) / Z2
                              t21 = pj1sd / Z2a
                              t22 = (pi0sd - tj * pi1sd) / Z2
                              t23 = pi1sd / Z2a
                              t24 = d1pl / Z2a
                              t25 = (d0pl  - ti * d1pl) / Zc

                              af_ind   = kfunct + Ll(l3) + l4 -1
                              rho_ind = ifunct + Ll(l1) + l2 -1 + Jx(jfunct)
                              te = af(af_ind) * rho_mat(rho_ind) * ccoef/(f1*f2)
                              ty = 2.0D0 * te

                              do l5 = 1, 3
                                 dpd = (Q(l5) - r(Nuc(jfunct),l5))  * dsd + &
                                       (W(l5) - Q(l5))              * ds1d
                                 dsf = (W(l5) - r(Nucd(kfunct),l5)) * ds1d

                                 if (l1 .eq. l5) then
                                    dpd = dpd + t20
                                    dsf = dsf + t21
                                    frc(Nuc(ifunct),l5) = frc(Nuc(ifunct),l5)  &
                                                        - te * pj0sd
                                 endif
                                 if (l2 .eq. l5) then
                                    dpd = dpd + t22
                                    dsf = dsf + t23
                                    frc(Nuc(ifunct),l5) = frc(Nuc(ifunct),l5)  &
                                                        - te * pi0sd
                                 endif
                                 if (l3 .eq. l5) then
                                    dpd = dpd + t24
                                    dsf = dsf + t25
                                    frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5)&
                                                         - te * d0pl
                                 endif
                                 if (l4 .eq. l5) then
                                    dpd = dpd + t26
                                    dsf = dsf + t27
                                    frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5)&
                                                         - te * ds0p
                                 endif
                                 fsd = dpd - dsd * &
                                       (r(Nuc(ifunct),l5) - r(Nuc(jfunct),l5))

                                 frc(Nuc(ifunct),l5)  = frc(Nuc(ifunct),l5)  +&
                                                        a(ifunct,nci)  * ty *fsd
                                 frc(Nuc(jfunct),l5)  = frc(Nuc(jfunct),l5)  +&
                                                        a(jfunct,ncj)  * ty *dpd
                                 frc(Nucd(kfunct),l5) = frc(Nucd(kfunct),l5) +&
                                                        ad(kfunct,nck) * ty *dsf
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
   !print*, "ds|d", frc

   ! (dp|d) and gradients
   do ifunct = ns+np+1, M    , 6
   do jfunct = ns+1   , ns+np, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij
            do kfunct = nsd+npd+1, Md, 6
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  ss5s = t2 * FUNCT(5,uf)
                  ss6s = t2 * FUNCT(6,uf)

                  ta  = (sss  - tj * ss1s) / Z2
                  t3  = (ss1s - tj * ss2s) / Z2
                  t4  = (ss2s - tj * ss3s) / Z2
                  t5  = (ss3s - tj * ss4s) / Z2
                  t5b = (ss4s - tj * ss5s) / Z2
                  t5x = ss2s / Z2a
                  t5y = ss3s / Z2a
                  do l1 = 1, 3
                     t1  = Q(l1) - r(Nuc(ifunct),l1)
                     t2  = W(l1) - Q(l1)
                     ps  = t1 * sss  + t2 * ss1s
                     p1s = t1 * ss1s + t2 * ss2s
                     p2s = t1 * ss2s + t2 * ss3s
                     p3s = t1 * ss3s + t2 * ss4s
                     p4s = t1 * ss4s + t2 * ss5s
                     p5s = t1 * ss5s + t2 * ss6s
                     t6a = (ps  - tj * p1s) / Z2
                     t6b = (p1s - tj * p2s) / Z2
                     t6c = (p2s - tj * p3s) / Z2
                     t6d = (p3s - tj * p4s) / Z2
                     t9  = p2s / Z2a
                     t9b = p3s / Z2a

                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        ds   = t1 * ps   + t2 * p1s
                        d1s  = t1 * p1s  + t2 * p2s
                        d2s  = t1 * p2s  + t2 * p3s
                        d3s  = t1 * p3s  + t2 * p4s
                        d4s  = t1 * p4s  + t2 * p5s
                        pjs  = t1 * sss  + t2 * ss1s
                        pj1s = t1 * ss1s + t2 * ss2s
                        pj2s = t1 * ss2s + t2 * ss3s
                        pj3s = t1 * ss3s + t2 * ss4s
                        pj4s = t1 * ss4s + t2 * ss5s
                        t7a  = (pjs  - tj * pj1s) / Z2
                        t7b  = (pj1s - tj * pj2s) / Z2
                        t7c  = (pj2s - tj * pj3s) / Z2
                        t7d  = (pj3s - tj * pj4s) / Z2
                        t8   = pj2s / Z2a
                        t8b  = pj3s / Z2a
                        f1   = 1.0D0

                        if (l1 .eq. l2) then
                           ds  = ds  + ta
                           d1s = d1s + t3
                           d2s = d2s + t4
                           d3s = d3s + t5
                           d4s = d4s + t5b
                           f1  = SQ3
                        endif
                        t10a = d1s / Z2a
                        t10  = d2s / Z2a
                        t10b = d3s / Z2a
                        t24  = (ds  - ti * d1s) / Zc
                        t24b = (d1s - ti * d2s) / Zc
                        t28  = d1s / Z2a
                        t28b = d2s / Z2a

                        do l3 = 1, 3
                           t1    = Q(l3) - r(Nuc(jfunct),l3)
                           t2    = W(l3) - Q(l3)
                           s2pks = t1 * ss2s + t2 * ss3s
                           s3pks = t1 * ss3s + t2 * ss4s
                           dp    = t1 * ds   + t2 * d1s
                           d1p   = t1 * d1s  + t2 * d2s
                           d2p   = t1 * d2s  + t2 * d3s
                           d3p   = t1 * d3s  + t2 * d4s
                           pi0p  = t1 * ps   + t2 * p1s
                           pi1p  = t1 * p1s  + t2 * p2s
                           pi2p  = t1 * p2s  + t2 * p3s
                           pi3p  = t1 * p3s  + t2 * p4s
                           pj0p  = t1 * pjs  + t2 * pj1s
                           pj1p  = t1 * pj1s + t2 * pj2s
                           pj2p  = t1 * pj2s + t2 * pj3s
                           pj3p  = t1 * pj3s + t2 * pj4s

                           if (l1 .eq. l3) then
                              dp   = dp   + t7a
                              d1p  = d1p  + t7b
                              d2p  = d2p  + t7c
                              d3p  = d3p  + t7d
                              pi0p = pi0p + ta
                              pi1p = pi1p + t3
                              pi2p = pi2p + t4
                              pi3p = pi3p + t5
                           endif
                           if (l2 .eq. l3) then
                              dp   = dp   + t6a
                              d1p  = d1p  + t6b
                              d2p  = d2p  + t6c
                              d3p  = d3p  + t6d
                              pj0p = pj0p + ta
                              pj1p = pj1p + t3
                              pj2p = pj2p + t4
                              pj3p = pj3p + t5
                           endif
                           t11a = pi1p  / Z2a
                           t11  = pi2p  / Z2a
                           t11b = pi3p  / Z2a
                           t12a = pj1p  / Z2a
                           t12  = pj2p  / Z2a
                           t12b = pj3p  / Z2a
                           t13  = s2pks / Z2a
                           t13b = s3pks / Z2a
                           t22  = (pj0p - ti * pj1p) / Zc
                           t22b = (pj1p - ti * pj2p) / Zc
                           t26  = pj1p  / Z2a
                           t26b = pj2p  / Z2a
                           t25  = (pi0p - ti * pi1p) / Zc
                           t25b = (pi1p - ti * pi2p) / Zc
                           t27  = pi1p  / Z2a
                           t27b = pi2p  / Z2a

                           do l4 = 1, 3
                              t1      = W(l4) - r(Nucd(kfunct),l4)
                              dp0p    = t1 * d1p
                              dp1p    = t1 * d2p
                              dp2p    = t1 * d3p
                              pjpkpl  = t1 * pj2p
                              pj2pkpl = t1 * pj3p
                              pipkpl  = t1 * pi2p
                              pi2pkpl = t1 * pi3p
                              dspl    = t1 * d2s
                              ds2pl   = t1 * d3s
                              pj1spl  = t1 * pj2s
                              pj2spl  = t1 * pj3s
                              pi1spl  = t1 * p2s
                              pi2spl  = t1 * p3s
                              s1pkpl  = t1 * s2pks
                              s2pkpl  = t1 * s3pks

                              if (l1 .eq. l4) then
                                 dp0p    = dp0p    + t12a
                                 dp1p    = dp1p    + t12
                                 dp2p    = dp2p    + t12b
                                 pipkpl  = pipkpl  + t13
                                 pi2pkpl = pi2pkpl + t13b
                                 ds2pl   = ds2pl   + t8b
                                 dspl    = dspl    + t8
                                 pi1spl  = pi1spl  + t5x
                                 pi2spl  = pi2spl  + t5y
                              endif
                              if (l2 .eq. l4) then
                                 dp0p    = dp0p    + t11a
                                 dp1p    = dp1p    + t11
                                 dp2p    = dp2p    + t11b
                                 pjpkpl  = pjpkpl  + t13
                                 pj2pkpl = pj2pkpl + t13b
                                 dspl    = dspl    + t9
                                 ds2pl   = ds2pl   + t9b
                                 pj1spl  = pj1spl  + t5x
                                 pj2spl  = pj2spl  + t5y
                              endif
                              if (l3 .eq. l4) then
                                 dp0p    = dp0p    + t10a
                                 dp1p    = dp1p    + t10
                                 dp2p    = dp2p    + t10b
                                 pipkpl  = pipkpl  + t9
                                 pi2pkpl = pi2pkpl + t9b
                                 pjpkpl  = pjpkpl  + t8
                                 pj2pkpl = pj2pkpl + t8b
                                 s1pkpl  = s1pkpl  + t5x
                                 s2pkpl  = s2pkpl  + t5y
                              endif
                              t14  = pjpkpl  / Z2a
                              t14b = pj2pkpl / Z2a
                              t15  = pipkpl  / Z2a
                              t15b = pi2pkpl / Z2a
                              t16  = dspl    / Z2a
                              t16b = ds2pl   / Z2a
                              t17  = (dp  - ti * d1p)  / Zc
                              t17b = (d1p - ti * d2p) / Zc
                              t20  = s1pkpl / Z2a
                              t20b = s2pkpl / Z2a
                              t21  = pj1spl / Z2a
                              t21b = pj2spl / Z2a
                              t23  = pi1spl / Z2a
                              t23b = pi2spl / Z2a
                              t38  = dp1p   / Z2a
                              t39  = (dp0p - ti * dp1p) / Zc

                              do l5 = 1, l4
                                 t1    = W(l5) - r(Nucd(kfunct),l5)
                                 dpd   = t1 * dp1p
                                 dp1d  = t1 * dp2p
                                 pip0d = t1 * pipkpl
                                 pip1d = t1 * pi2pkpl
                                 pjp0d = t1 * pjpkpl
                                 pjp1d = t1 * pj2pkpl
                                 d0d   = t1 * dspl
                                 d1d   = t1 * ds2pl
                                 dp0pm = t1 * d1p
                                 dp1pm = t1 * d2p
                                 f2    = 1.0D0

                                 if (l1 .eq. l5) then
                                    dpd   = dpd   + t14
                                    dp1d  = dp1d  + t14b
                                    pip0d = pip0d + t20
                                    pip1d = pip1d + t20b
                                    d0d   = d0d   + t21
                                    d1d   = d1d   + t21b
                                    dp0pm = dp0pm + t26
                                    dp1pm = dp1pm + t26b
                                 endif
                                 if (l2 .eq. l5) then
                                    dpd   = dpd   + t15
                                    dp1d  = dp1d  + t15b
                                    pjp0d = pjp0d + t20
                                    pjp1d = pjp1d + t20b
                                    d0d   = d0d   + t23
                                    d1d   = d1d   + t23b
                                    dp0pm = dp0pm + t27
                                    dp1pm = dp1pm + t27b
                                 endif
                                 if (l3 .eq. l5) then
                                    dpd   = dpd   + t16
                                    dp1d  = dp1d  + t16b
                                    pjp0d = pjp0d + t21
                                    pjp1d = pjp1d + t21b
                                    pip0d = pip0d + t23
                                    pip1d = pip1d + t23b
                                    dp0pm = dp0pm + t28
                                    dp1pm = dp1pm + t28b
                                 endif
                                 if (l4 .eq. l5) then
                                    dpd   = dpd   + t17
                                    dp1d  = dp1d  + t17b
                                    pjp0d = pjp0d + t22
                                    pjp1d = pjp1d + t22b
                                    pip0d = pip0d + t25
                                    pip1d = pip1d + t25b
                                    d0d   = d0d   + t24
                                    d1d   = d1d   + t24b
                                    f2    = SQ3
                                 endif
                                 t30 = (pjp0d - tj * pjp1d) / Z2
                                 t31 = pjp1d / Z2a
                                 t32 = (pip0d - tj * pip1d) / Z2
                                 t33 = pip1d / Z2a
                                 t34 = (d0d   - tj * d1d)   / Z2
                                 t35 = d1d   / Z2a
                                 t36 = dp1pm / Z2a
                                 t37 = (dp0pm - ti * dp1pm) / Zc

                                 af_ind   = kfunct + Ll(l4) + l5 -1
                                 rho_ind = ifunct + Ll(l1) + l2 -1 + &
                                            Jx(jfunct + l3 -1)
                                 te = af(af_ind) * rho_mat(rho_ind) * ccoef / &
                                      (f1 * f2)
                                 ty = 2.0D0 * te

                                 do l6 = 1, 3
                                    ddd = (Q(l6) - r(Nuc(jfunct),l6))  * dpd + &
                                          (W(l6) - Q(l6))              * dp1d
                                    dpf = (W(l6) - r(Nucd(kfunct),l6)) * dp1d

                                    if (l1 .eq. l6) then
                                       ddd = ddd + t30
                                       dpf = dpf + t31
                                       frc(Nuc(ifunct),l6)= frc(Nuc(ifunct),l6)&
                                                          - te * pjp0d
                                    endif
                                    if (l2 .eq. l6) then
                                       ddd = ddd + t32
                                       dpf = dpf + t33
                                       frc(Nuc(ifunct),l6)= frc(Nuc(ifunct),l6)&
                                                          - te * pip0d
                                    endif
                                    if (l3 .eq. l6) then
                                       ddd = ddd + t34
                                       dpf = dpf + t35
                                       frc(Nuc(jfunct),l6)= frc(Nuc(jfunct),l6)&
                                                          - te * d0d
                                    endif
                                    if (l4 .eq. l6) then
                                       ddd = ddd + t36
                                       dpf = dpf + t37
                                       frc(Nucd(kfunct),l6)=&
                                                           frc(Nucd(kfunct),l6)&
                                                            - te * dp0pm
                                    endif
                                    if (l5 .eq. l6) then
                                       ddd = ddd + t38
                                       dpf = dpf + t39
                                       frc(Nucd(kfunct),l6)=&
                                                           frc(Nucd(kfunct),l6)&
                                                            - te * dp0p
                                    endif
                                    fpd = ddd - dpd * &
                                          (r(Nuc(ifunct),l6) -r(Nuc(jfunct),l6))

                                    frc(Nuc(ifunct),l6)  = frc(Nuc(ifunct),l6) &
                                                         + a(ifunct,nci)  * ty &
                                                         * fpd
                                    frc(Nuc(jfunct),l6)  = frc(Nuc(jfunct),l6) &
                                                         + a(jfunct,ncj)  * ty &
                                                         * ddd
                                    frc(Nucd(kfunct),l6) = frc(Nucd(kfunct),l6)&
                                                         + ad(kfunct,nck) * ty &
                                                         * dpf
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
   enddo
   enddo
   !print*, "dp|d", frc

   ! (dd|d) and gradients
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+np+1, ifunct, 6
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         tj   = a(jfunct,ncj) / Zij
         rexp = a(ifunct,nci) * tj * d(Nuc(ifunct),Nuc(jfunct))

         if (rexp .lt. rmax) then
            Z2   = 2.D0 * Zij
            ti   = a(ifunct,nci) / Zij
            Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
            Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
            Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
            sks  = pi52 * exp(-rexp) / Zij

            do kfunct = nsd+npd+1, Md, 6
               dpc = (Q(1) - r(Nucd(kfunct),1)) * (Q(1) - r(Nucd(kfunct),1)) + &
                     (Q(2) - r(Nucd(kfunct),2)) * (Q(2) - r(Nucd(kfunct),2)) + &
                     (Q(3) - r(Nucd(kfunct),3)) * (Q(3) - r(Nucd(kfunct),3))
               do nck = 1, ncontd(kfunct)
                  ccoef = c(ifunct,nci) * c(jfunct,ncj) * cd(kfunct,nck)
                  t0    = ad(kfunct,nck) + Zij
                  Z2a   = 2.0D0 * t0
                  Zc    = 2.0D0 * ad(kfunct,nck)
                  ti    = Zij / t0
                  tj    = ad(kfunct,nck) / t0
                  W(1)  = ti * Q(1) + tj * r(Nucd(kfunct),1)
                  W(2)  = ti * Q(2) + tj * r(Nucd(kfunct),2)
                  W(3)  = ti * Q(3) + tj * r(Nucd(kfunct),3)

                  uf   = tj * Zij * dpc
                  t2   = sks / (ad(kfunct,nck) * sqrt(t0))
                  sss  = t2 * FUNCT(0,uf)
                  ss1s = t2 * FUNCT(1,uf)
                  ss2s = t2 * FUNCT(2,uf)
                  ss3s = t2 * FUNCT(3,uf)
                  ss4s = t2 * FUNCT(4,uf)
                  ss5s = t2 * FUNCT(5,uf)
                  ss6s = t2 * FUNCT(6,uf)
                  ss7s = t2 * FUNCT(7,uf)
                  t3   = (sss  - tj * ss1s) / Z2
                  t4   = (ss1s - tj * ss2s) / Z2
                  t5   = (ss2s - tj * ss3s) / Z2
                  t6   = (ss3s - tj * ss4s) / Z2
                  t6b  = (ss4s - tj * ss5s) / Z2
                  t6c  = (ss5s - tj * ss6s) / Z2

                  do l1 = 1, 3
                     t1   = Q(l1) - r(Nuc(ifunct),l1)
                     t2   = W(l1) - Q(l1)
                     ps   = t1 * sss  + t2 * ss1s
                     p1s  = t1 * ss1s + t2 * ss2s
                     p2s  = t1 * ss2s + t2 * ss3s
                     p3s  = t1 * ss3s + t2 * ss4s
                     p4s  = t1 * ss4s + t2 * ss5s
                     p5s  = t1 * ss5s + t2 * ss6s
                     p6s  = t1 * ss6s + t2 * ss7s
                     t7   = (ps  - tj * p1s) / Z2
                     t8   = (p1s - tj * p2s) / Z2
                     t9   = (p2s - tj * p3s) / Z2
                     t10  = (p3s - tj * p4s) / Z2
                     t10b = (p4s - tj * p5s) / Z2
                     y16  = p2s / Z2a
                     y16b = p3s / Z2a

                     do l2 = 1, l1
                        t1   = Q(l2) - r(Nuc(ifunct),l2)
                        t2   = W(l2) - Q(l2)
                        pjs  = t1 * sss  + t2 * ss1s
                        pj1s = t1 * ss1s + t2 * ss2s
                        pj2s = t1 * ss2s + t2 * ss3s
                        pj3s = t1 * ss3s + t2 * ss4s
                        pj4s = t1 * ss4s + t2 * ss5s
                        pj5s = t1 * ss5s + t2 * ss6s
                        ds   = t1 * ps   + t2 * p1s
                        d1s  = t1 * p1s  + t2 * p2s
                        d2s  = t1 * p2s  + t2 * p3s
                        d3s  = t1 * p3s  + t2 * p4s
                        d4s  = t1 * p4s  + t2 * p5s
                        d5s  = t1 * p5s  + t2 * p6s
                        t11  = (pjs  - tj * pj1s) / Z2
                        t12  = (pj1s - tj * pj2s) / Z2
                        t13  = (pj2s - tj * pj3s) / Z2
                        t14  = (pj3s - tj * pj4s) / Z2
                        t14b = (pj4s - tj * pj5s) / Z2
                        y19  = pj2s / Z2a
                        y19b = pj3s / Z2a
                        f1   = 1.0D0

                        if (l1 .eq. l2) then
                           ds  = ds  + t3
                           d1s = d1s + t4
                           d2s = d2s + t5
                           d3s = d3s + t6
                           d4s = d4s + t6b
                           d5s = d5s + t6c
                           f1  = SQ3
                        endif
                        t16  = (ds  - tj * d1s) / Z2
                        t17  = (d1s - tj * d2s) / Z2
                        t18  = (d2s - tj * d3s) / Z2
                        t18b = (d3s - tj * d4s) / Z2
                        t22a = d2s / Z2a
                        t22c = d3s / Z2a

                        lij = 3
                        if (ifunct .eq. jfunct) lij = l1
                        do l3 = 1, lij
                           t1    = Q(l3) - r(Nuc(jfunct),l3)
                           t2    = W(l3) - Q(l3)
                           dpk   = t1 * ds   + t2 * d1s
                           d1pk  = t1 * d1s  + t2 * d2s
                           d2pk  = t1 * d2s  + t2 * d3s
                           d3pk  = t1 * d3s  + t2 * d4s
                           d4pk  = t1 * d4s  + t2 * d5s
                           pjpk  = t1 * pjs  + t2 * pj1s
                           pj1pk = t1 * pj1s + t2 * pj2s
                           pj2pk = t1 * pj2s + t2 * pj3s
                           pj3pk = t1 * pj3s + t2 * pj4s
                           pj4pk = t1 * pj4s + t2 * pj5s
                           pipk  = t1 * ps   + t2 * p1s
                           pi1pk = t1 * p1s  + t2 * p2s
                           pi2pk = t1 * p2s  + t2 * p3s
                           pi3pk = t1 * p3s  + t2 * p4s
                           pi4pk = t1 * p4s  + t2 * p5s
                           spk   = t1 * sss  + t2 * ss1s
                           s1pk  = t1 * ss1s + t2 * ss2s
                           s2pk  = t1 * ss2s + t2 * ss3s
                           s3pk  = t1 * ss3s + t2 * ss4s
                           s4pk  = t1 * ss4s + t2 * ss5s
                           t15p  = (spk  - tj * s1pk) / Z2
                           t15a  = (s1pk - tj * s2pk) / Z2
                           t15   = (s2pk - tj * s3pk) / Z2
                           t15b  = (s3pk - tj * s4pk) / Z2
                           y18   = s2pk / Z2a
                           y18b  = s3pk / Z2a

                           if (l1 .eq. l3) then
                              dpk   = dpk   + t11
                              d1pk  = d1pk  + t12
                              d2pk  = d2pk  + t13
                              d3pk  = d3pk  + t14
                              d4pk  = d4pk  + t14b
                              pipk  = pipk  + t3
                              pi1pk = pi1pk + t4
                              pi2pk = pi2pk + t5
                              pi3pk = pi3pk + t6
                              pi4pk = pi4pk + t6b
                           endif
                           if (l2 .eq. l3) then
                              dpk   = dpk   + t7
                              d1pk  = d1pk  + t8
                              d2pk  = d2pk  + t9
                              d3pk  = d3pk  + t10
                              d4pk  = d4pk  + t10b
                              pjpk  = pjpk  + t3
                              pj1pk = pj1pk + t4
                              pj2pk = pj2pk + t5
                              pj3pk = pj3pk + t6
                              pj4pk = pj4pk + t6b
                           endif
                           t20  = pj2pk / Z2a
                           t20b = pj3pk / Z2a
                           t21  = pi2pk / Z2a
                           t21b = pi3pk / Z2a
                           t22p = d1pk  / Z2a
                           t22  = d2pk  / Z2a
                           t22b = d3pk  / Z2a
                           t24  = (pjpk  - tj * pj1pk) / Z2
                           t25  = (pj1pk - tj * pj2pk) / Z2
                           t26  = (pj2pk - tj * pj3pk) / Z2
                           t26b = (pj3pk - tj * pj4pk) / Z2
                           t27  = (pipk  - tj * pi1pk) / Z2
                           t28  = (pi1pk - tj * pi2pk) / Z2
                           t29  = (pi2pk - tj * pi3pk) / Z2
                           t29b = (pi3pk - tj * pi4pk) / Z2
                           y14  = (dpk  - ti * d1pk) / Zc
                           y14b = (d1pk - ti * d2pk) / Zc

                           lk = l3
                           if (ifunct .eq. jfunct) lk = min(l3,Ll(l1)-Ll(l3)+l2)
                           do l4 = 1, lk
                              t1 = Q(l4) - r(Nuc(jfunct),l4)
                              t2 = W(l4) - Q(l4)
                              dd     = t1 * dpk   + t2 * d1pk
                              d1d    = t1 * d1pk  + t2 * d2pk
                              d2d    = t1 * d2pk  + t2 * d3pk
                              d3d    = t1 * d3pk  + t2 * d4pk
                              pj0dkl = t1 * pjpk  + t2 * pj1pk
                              pj1dkl = t1 * pj1pk + t2 * pj2pk
                              pj2dkl = t1 * pj2pk + t2 * pj3pk
                              pj3dkl = t1 * pj3pk + t2 * pj4pk
                              pi0dkl = t1 * pipk  + t2 * pi1pk
                              pi1dkl = t1 * pi1pk + t2 * pi2pk
                              pi2dkl = t1 * pi2pk + t2 * pi3pk
                              pi3dkl = t1 * pi3pk + t2 * pi4pk
                              d0pl   = t1 * ds    + t2 * d1s
                              d1pl   = t1 * d1s   + t2 * d2s
                              d2pl   = t1 * d2s   + t2 * d3s
                              d3pl   = t1 * d3s   + t2 * d4s
                              s2pl   = t1 * ss2s  + t2 * ss3s
                              s3pl   = t1 * ss3s  + t2 * ss4s
                              s2dkl  = t1 * s2pk  + t2 * s3pk
                              s3dkl  = t1 * s3pk  + t2 * s4pk
                              pj2pl  = t1 * pj2s  + t2 * pj3s
                              pj3pl  = t1 * pj3s  + t2 * pj4s
                              pi2pl  = t1 * p2s   + t2 * p3s
                              pi3pl  = t1 * p3s   + t2 * p4s
                              y17    = s2pl / Z2a
                              y17b   = s3pl / Z2a
                              f2     = 1.0D0

                              if (l1 .eq. l4) then
                                 dd     = dd     + t24
                                 d1d    = d1d    + t25
                                 d2d    = d2d    + t26
                                 d3d    = d3d    + t26b
                                 pi0dkl = pi0dkl + t15p
                                 pi1dkl = pi1dkl + t15a
                                 pi2dkl = pi2dkl + t15
                                 pi3dkl = pi3dkl + t15b
                                 d0pl   = d0pl   + t11
                                 d1pl   = d1pl   + t12
                                 d2pl   = d2pl   + t13
                                 d3pl   = d3pl   + t14
                                 pi2pl  = pi2pl  + t5
                                 pi3pl  = pi3pl  + t6
                              endif
                              if (l2 .eq. l4) then
                                 dd     = dd     + t27
                                 d1d    = d1d    + t28
                                 d2d    = d2d    + t29
                                 d3d    = d3d    + t29b
                                 pj0dkl = pj0dkl + t15p
                                 pj1dkl = pj1dkl + t15a
                                 pj2dkl = pj2dkl + t15
                                 pj3dkl = pj3dkl + t15b
                                 d0pl   = d0pl   + t7
                                 d1pl   = d1pl   + t8
                                 d2pl   = d2pl   + t9
                                 d3pl   = d3pl   + t10
                                 pj2pl  = pj2pl  + t5
                                 pj3pl  = pj3pl  + t6
                              endif
                              if (l3 .eq. l4) then
                                 s2dkl  = s2dkl  + t5
                                 s3dkl  = s3dkl  + t6
                                 dd     = dd     + t16
                                 d1d    = d1d    + t17
                                 d2d    = d2d    + t18
                                 d3d    = d3d    + t18b
                                 pj0dkl = pj0dkl + t11
                                 pj2dkl = pj2dkl + t13
                                 pj1dkl = pj1dkl + t12
                                 pj3dkl = pj3dkl + t14
                                 pi3dkl = pi3dkl + t10
                                 pi2dkl = pi2dkl + t9
                                 pi1dkl = pi1dkl + t8
                                 pi0dkl = pi0dkl + t7
                                 f2     = SQ3
                              endif
                              t30a   = pj1dkl / Z2a
                              t30    = pj2dkl / Z2a
                              t30b   = pj3dkl / Z2a
                              t40a   = pi1dkl / Z2a
                              t40    = pi2dkl / Z2a
                              t40b   = pi3dkl / Z2a
                              t50    = s2dkl  / Z2a
                              t50b   = s3dkl  / Z2a
                              t60    = pj2pl  / Z2a
                              t60b   = pj3pl  / Z2a
                              t70    = pi2pl  / Z2a
                              t70b   = pi3pl  / Z2a
                              t80a   = d1pl   / Z2a
                              t80    = d2pl   / Z2a
                              t80b   = d3pl   / Z2a
                              t23    = (dd  - ti * d1d) / Zc
                              t23b   = (d1d - ti * d2d) / Zc

                              y12  = (pi0dkl - ti * pi1dkl) / Zc
                              y12b = (pi1dkl - ti * pi2dkl) / Zc
                              y13  = (d0pl   - ti * d1pl)   / Zc
                              y13b = (d1pl   - ti * d2pl)   / Zc
                              y15  = (pj0dkl - ti * pj1dkl) / Zc
                              y15b = (pj1dkl - ti * pj2dkl) / Zc

                              do l5 = 1, 3
                                 t1 = W(l5) - r(Nucd(kfunct),l5)
                                 dd0p    = t1 * d1d
                                 ddp     = t1 * d2d
                                 dd2p    = t1 * d3d
                                 pjdklp  = t1 * pj2dkl
                                 pj2dklp = t1 * pj3dkl
                                 pidklp  = t1 * pi2dkl
                                 pi2dklp = t1 * pi3dkl
                                 dijplp  = t1 * d2pl
                                 dij2plp = t1 * d3pl
                                 dijpkp  = t1 * d2pk
                                 dij2pkp = t1 * d3pk
                                 s1dpm   = t1 * s2dkl
                                 s2dpm   = t1 * s3dkl
                                 pj1plpm = t1 * pj2pl
                                 pj2plpm = t1 * pj3pl
                                 pj1pkpm = t1 * pj2pk
                                 pj2pkpm = t1 * pj3pk
                                 pi1plpm = t1 * pi2pl
                                 pi2plpm = t1 * pi3pl
                                 pi1pkpm = t1 * pi2pk
                                 pi2pkpm = t1 * pi3pk
                                 d1spm   = t1 * d2s
                                 d2spm   = t1 * d3s

                                 if (l1 .eq. l5) then
                                    dd0p    = dd0p    + t30a
                                    ddp     = ddp     + t30
                                    dd2p    = dd2p    + t30b
                                    pidklp  = pidklp  + t50
                                    pi2dklp = pi2dklp + t50b
                                    dijplp  = dijplp  + t60
                                    dij2plp = dij2plp + t60b
                                    dijpkp  = dijpkp  + t20
                                    dij2pkp = dij2pkp + t20b
                                    pi1plpm = pi1plpm + y17
                                    pi2plpm = pi2plpm + y17b
                                    pi1pkpm = pi1pkpm + y18
                                    pi2pkpm = pi2pkpm + y18b
                                    d1spm   = d1spm   + y19
                                    d2spm   = d2spm   + y19b
                                 endif
                                 if (l2 .eq. l5) then
                                    dd0p    = dd0p    + t40a
                                    ddp     = ddp     + t40
                                    dd2p    = dd2p    + t40b
                                    pjdklp  = pjdklp  + t50
                                    pj2dklp = pj2dklp + t50b
                                    dijplp  = dijplp  + t70
                                    dij2plp = dij2plp + t70b
                                    dijpkp  = dijpkp  + t21
                                    dij2pkp = dij2pkp + t21b
                                    pj1plpm = pj1plpm + y17
                                    pj2plpm = pj2plpm + y17b
                                    pj1pkpm = pj1pkpm + y18
                                    pj2pkpm = pj2pkpm + y18b
                                    d1spm   = d1spm   + y16
                                    d2spm   = d2spm   + y16b
                                 endif
                                 if (l3 .eq. l5) then
                                    dd0p    = dd0p    + t80a
                                    ddp     = ddp     + t80
                                    dd2p    = dd2p    + t80b
                                    pjdklp  = pjdklp  + t60
                                    pj2dklp = pj2dklp + t60b
                                    pidklp  = pidklp  + t70
                                    pi2dklp = pi2dklp + t70b
                                    dijpkp  = dijpkp  + t22a
                                    dij2pkp = dij2pkp + t22c
                                    s1dpm   = s1dpm   + y17
                                    s2dpm   = s2dpm   + y17b
                                    pj1pkpm = pj1pkpm + y19
                                    pj2pkpm = pj2pkpm + y19b
                                    pi1pkpm = pi1pkpm + y16
                                    pi2pkpm = pi2pkpm + y16b
                                 endif
                                 if (l4 .eq. l5) then
                                    dd0p    = dd0p    + t22p
                                    ddp     = ddp     + t22
                                    dd2p    = dd2p    + t22b
                                    pjdklp  = pjdklp  + t20
                                    pj2dklp = pj2dklp + t20b
                                    pidklp  = pidklp  + t21
                                    pi2dklp = pi2dklp + t21b
                                    dijplp  = dijplp  + t22a
                                    dij2plp = dij2plp + t22c
                                    s1dpm   = s1dpm   + y18
                                    s2dpm   = s2dpm   + y18b
                                    pj1plpm = pj1plpm + y19
                                    pj2plpm = pj2plpm + y19b
                                    pi1plpm = pi1plpm + y16
                                    pi2plpm = pi2plpm + y16b
                                 endif
                                 t31  = pjdklp  / Z2a
                                 t31b = pj2dklp / Z2a
                                 t41  = pidklp  / Z2a
                                 t41b = pi2dklp / Z2a
                                 t51  = dijplp  / Z2a
                                 t51b = dij2plp / Z2a
                                 t61  = dijpkp  / Z2a
                                 t61b = dij2pkp / Z2a
                                 y30  = ddp     / Z2a
                                 y31  = (dd0p - ti * ddp) / Zc
                                 y2   = s1dpm   / Z2a
                                 y2b  = s2dpm   / Z2a
                                 y3   = pj1plpm / Z2a
                                 y3b  = pj2plpm / Z2a
                                 y4   = pj1pkpm / Z2a
                                 y4b  = pj2pkpm / Z2a
                                 y6   = pi1plpm / Z2a
                                 y6b  = pi2plpm / Z2a
                                 y7   = pi1pkpm / Z2a
                                 y7b  = pi2pkpm / Z2a
                                 y9   = d1spm   / Z2a
                                 y9b  = d2spm   / Z2a

                                 do l6 = 1, l5
                                    t1 = W(l6) - r(Nucd(kfunct),l6)
                                    ddd   = t1 * ddp
                                    dd1d  = t1 * dd2p
                                    dd0pn = t1 * d1d
                                    dd1pn = t1 * d2d
                                    pj0dd = t1 * pjdklp
                                    pj1dd = t1 * pj2dklp
                                    pi0dd = t1 * pidklp
                                    pi1dd = t1 * pi2dklp
                                    d0pld = t1 * dijplp
                                    d1pld = t1 * dij2plp
                                    d0pkd = t1 * dijpkp
                                    d1pkd = t1 * dij2pkp
                                    f3    = 1.0D0

                                    if (l1 .eq. l6) then
                                       ddd   = ddd   + t31
                                       dd1d  = dd1d  + t31b
                                       dd0pn = dd0pn + t30a
                                       dd1pn = dd1pn + t30
                                       pi0dd = pi0dd + y2
                                       pi1dd = pi1dd + y2b
                                       d0pld = d0pld + y3
                                       d1pld = d1pld + y3b
                                       d0pkd = d0pkd + y4
                                       d1pkd = d1pkd + y4b
                                    endif
                                    if (l2 .eq. l6) then
                                       ddd   = ddd   + t41
                                       dd1d  = dd1d  + t41b
                                       dd0pn = dd0pn + t40a
                                       dd1pn = dd1pn + t40
                                       d0pld = d0pld + y6
                                       d1pld = d1pld + y6b
                                       d0pkd = d0pkd + y7
                                       d1pkd = d1pkd + y7b
                                       pj0dd = pj0dd + y2
                                       pj1dd = pj1dd + y2b
                                    endif
                                    if (l3 .eq. l6) then
                                       ddd   = ddd   + t51
                                       dd1d  = dd1d  + t51b
                                       dd0pn = dd0pn + t80a
                                       dd1pn = dd1pn + t80
                                       pj0dd = pj0dd + y3
                                       pj1dd = pj1dd + y3b
                                       d0pkd = d0pkd + y9
                                       d1pkd = d1pkd + y9b
                                       pi0dd = pi0dd + y6
                                       pi1dd = pi1dd + y6b
                                    endif
                                    if (l4 .eq. l6) then
                                       ddd   = ddd   + t61
                                       dd1d  = dd1d  + t61b
                                       pi0dd = pi0dd + y7
                                       pi1dd = pi1dd + y7b
                                       d0pld = d0pld + y9
                                       d1pld = d1pld + y9b
                                       dd0pn = dd0pn + t22p
                                       dd1pn = dd1pn + t22
                                       pj0dd = pj0dd + y4
                                       pj1dd = pj1dd + y4b
                                    endif
                                    if (l5 .eq. l6) then
                                       ddd   = ddd   + t23
                                       dd1d  = dd1d  + t23b
                                       pi0dd = pi0dd + y12
                                       pi1dd = pi1dd + y12b
                                       d0pld = d0pld + y13
                                       d1pld = d1pld + y13b
                                       d0pkd = d0pkd + y14
                                       d1pkd = d1pkd + y14b
                                       pj0dd = pj0dd + y15
                                       pj1dd = pj1dd + y15b
                                       f3    = SQ3
                                    endif
                                    y20 = (pj0dd - tj * pj1dd) / Z2
                                    y21 = pj1dd / Z2a
                                    y22 = (pi0dd - tj * pi1dd) / Z2
                                    y23 = pi1dd / Z2a
                                    y24 = (d0pld - tj * d1pld) / Z2
                                    y25 = d1pld / Z2a
                                    y26 = (d0pkd - tj * d1pkd) / Z2
                                    y27 = d1pkd / Z2a
                                    y28 = dd1pn / Z2a
                                    y29 = (dd0pn - ti * dd1pn) / Zc

                                    af_ind   = kfunct + Ll(l5) + l6 -1
                                    rho_ind = ifunct + Ll(l1) + l2 -1 + &
                                               Jx(jfunct + Ll(l3) + l4 -1)
                                    te = af(af_ind) * rho_mat(rho_ind) * ccoef/&
                                         (f1 * f2 * f3)
                                    ty = 2.0D0 * te

                                    do l7 = 1, 3
                                       fdd = (Q(l7) - r(Nuc(ifunct),l7))  *ddd &
                                           + (W(l7) - Q(l7))              *dd1d
                                       ddf = (W(l7) - r(Nucd(kfunct),l7)) *dd1d

                                       if (l1 .eq. l7) then
                                          fdd = fdd + y20
                                          ddf = ddf + y21
                                          frc(Nuc(ifunct),l7) = &
                                                frc(Nuc(ifunct),l7) - te * pj0dd
                                       endif
                                       if (l2 .eq. l7) then
                                          fdd = fdd + y22
                                          ddf = ddf + y23
                                          frc(Nuc(ifunct),l7) = &
                                                frc(Nuc(ifunct),l7) - te * pi0dd
                                       endif
                                       if (l3 .eq. l7) then
                                          fdd = fdd + y24
                                          ddf = ddf + y25
                                          frc(Nuc(jfunct),l7) = &
                                                frc(Nuc(jfunct),l7) - te * d0pld
                                       endif
                                       if (l4 .eq. l7) then
                                          fdd = fdd + y26
                                          ddf = ddf + y27
                                          frc(Nuc(jfunct),l7) = &
                                                frc(Nuc(jfunct),l7) - te * d0pkd
                                       endif
                                       if (l5 .eq. l7) then
                                          fdd = fdd + y28
                                          ddf = ddf + y29
                                          frc(Nucd(kfunct),l7) = &
                                               frc(Nucd(kfunct),l7) - te * dd0pn
                                       endif
                                       if (l6 .eq. l7) then
                                          fdd = fdd + y30
                                          ddf = ddf + y31
                                          frc(Nucd(kfunct),l7) = &
                                                frc(Nucd(kfunct),l7) - te * dd0p
                                       endif
                                       dfd = fdd + ddd * &
                                           (r(Nuc(ifunct),l7)-r(Nuc(jfunct),l7))

                                       frc(Nuc(ifunct),l7)  = &
                                                         frc(Nuc(ifunct),l7)  +&
                                                         a(ifunct,nci) * ty *fdd
                                       frc(Nuc(jfunct),l7)  = &
                                                         frc(Nuc(jfunct),l7)  +&
                                                         a(jfunct,ncj) * ty *dfd
                                       frc(Nucd(kfunct),l7) = &
                                                         frc(Nucd(kfunct),l7) +&
                                                         ad(kfunct,nck)* ty *ddf
                                    enddo
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
   enddo
   enddo
   !print*, "dd|d", frc

   deallocate(Jx, Ll, Q, W)
   call g2g_timer_stop('CoulG')
   call g2g_timer_sum_stop('Coulomb gradients')
   return
end subroutine
end module subm_int3G
