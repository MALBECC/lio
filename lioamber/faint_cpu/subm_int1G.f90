module subm_int1G
contains
subroutine int1G(ff, rho, d, r, Iz, natom, ntatom)
!------------------------------------------------------------------------------!
! Calculates 1e gradients, to be used with MD, using the Obara-Saika recursive !
! method.                                                                      !
! Loops over all basis functions, ordered according to the type (all s, then   !
! all p, then all d). Inside each type, they are ordered in shells.            !
! Inputs : basis function and system information.                              !
! Outputs: forces (gradients) on nuclei.                                       !
!                                                                              !
! Relevant internal variables:                                                 !
! ns: marker for end of s functions.                                           !
! np: marker for end of p functions.                                           !
! nd: marker for end of d functions.                                           !
! r(Nuc(ifunct),j): j component of position of nucleus for basis i (j=1,2,3).  !
! Original and debugged (or supposed to) by Dario Estrin on 28/07/1992         !
! Refactored by Federico Pedron on 25/07/2018                                  !
!------------------------------------------------------------------------------!
   use garcha_mod   , only: a, c, Nuc, ncont, nshell, NORM, M
   use liotemp      , only: FUNCT
   use constants_mod, only: pi, pi32

   implicit none
   ! Inputs and Outputs
   integer         , intent(in)  :: natom, ntatom, Iz(natom)
   double precision, intent(in)  :: d(natom,natom), r(ntatom,3), rho(:)
   double precision, intent(out) :: ff(natom,3)

   ! Auxiliary variables
   integer :: my_natom, igpu, ns, np, nd, M2
   integer :: l1, l2, l3, l4, l5, l12, l34, lk, lij, k_ind, nci, ncj, ifunct, &
              jfunct, iatom, jatom, ll(3)

   double precision :: Q(3), ovlap, alf, alf2, alf3, alf4, temp, temp0, sq3, &
                       ccoef, f1, f2, tn, tna, tn1a, Z2, Zij, q1, q2, q3, uf,&
                       tt, tx, te, ti, tj

   double precision :: ss, sks, spk, spj, skpj, s1p, s0p, &
                       ps, pp, p0s, p1s, p2s, p3s, p4s, pi0p, pi1p, pi2p, pi0d,&
                       pi1d, piks, pikpk, pipk, pis, piNs, pikdkl, pidkl, pj0s,&
                       pj1s, pj2s, pj3s, pj0p, pj1p, pj2p, pj0d, pj1d, pjkpk,  &
                       pjks, pjpk, pjs, pjkdkl, pjdkl, pks, pkp, pkd, pNd, pNp,&
                       pN1p, &
                       ds, dp, dd, d0s, d1s, d2s, d3s, d0p, d1p, d2p, d0pl,  &
                       d1pl, dijs, dijpk, dijks, dijkpk, dijkpl, dijpl, dks, &
                       dkp, dkd, dkf, dN1s, dNs, dNp, dNd, dNf, dsd, &
                       fNs, fNd, fNp, fkd, fkp, fks, fd

   double precision :: t0, t1, t2, t3, t4, t5, t7, t8, t9, t40, t41, &
                       t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, &
                       t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, &
                       t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, &
                       t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, &
                       t60, t61, t62, t63, t64, t65, t66, t67, t68, t69, &
                       t70, t71, t72, t73, t74, t81, t82, t83, t84, t85, &
                       t86, t81b, t82b, t83b, t84b, t85b, t86b, &
                       t90, t91, t92, t93, t94, t95, t96, t97, t98

   double precision, dimension(3) :: dn, dn1, dn2, dn3, dn4, dn5, dn6, dn7,  &
                                     dn8, dn9, dn10, dn2b, dn4b, dn5b, dn7b, &
                                     dn8b, dn9b
   double precision, allocatable  :: s0s(:)  , s1s(:)  , s2s(:)  , s3s(:)  , &
                                     s4s(:)  , s5s(:)  , s6s(:)  , x0x(:,:), &
                                     x1x(:,:), x2x(:,:), x3x(:,:), x4x(:,:)

   allocate(s0s(natom)  , s1s(natom)  , s2s(natom)  , s3s(natom)  , s4s(natom),&
            s5s(natom)  , s6s(natom)  , x0x(natom,3), x1x(natom,3),            &
            x2x(natom,3), x3x(natom,3), x4x(natom,3))

   ! Checks basis set normalization.
   if (NORM) then
      sq3 = sqrt(3.D0)
   else
      sq3 = 1.D0
   endif

   ns = nshell(0); np  = nshell(1); nd = nshell(2); M2 = 2*M

   do l1 = 1,3
      Ll(l1) = l1*(l1-1)/2
   enddo

   ! Overlap ,Kinetic energy and Nuclear Attraction matrix elements evaluation.
   ! Overlap matrix will be kept, while kinetic energy and nuclear attraction
   ! are directly stored in Fock and Energy matrices.
   ff = 0.0D0
   do iatom = 1, natom
      do jatom = 1, iatom-1
         tt = Iz(iatom)*Iz(jatom) / d(iatom,jatom)**1.5D0
         ff(iatom,1) = ff(iatom,1) - tt*(r(iatom,1)-r(jatom,1))
         ff(iatom,2) = ff(iatom,2) - tt*(r(iatom,2)-r(jatom,2))
         ff(iatom,3) = ff(iatom,3) - tt*(r(iatom,3)-r(jatom,3))
      enddo

      do jatom = iatom+1, natom
         tt = Iz(iatom)*Iz(jatom) / d(iatom,jatom)**1.5D0
         ff(iatom,1) = ff(iatom,1) - tt*(r(iatom,1)-r(jatom,1))
         ff(iatom,2) = ff(iatom,2) - tt*(r(iatom,2)-r(jatom,2))
         ff(iatom,3) = ff(iatom,3) - tt*(r(iatom,3)-r(jatom,3))
      enddo
   enddo

   ! Checks if the 1e integrals are done in GPU, doing only the KE terms if so.
   call aint_query_gpu_level(igpu)
   my_natom = natom
   if (igpu.gt.3) my_natom = 0

   ! First loop - (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ! (0|0) calculation
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         alf  = a(ifunct,nci) * a(jfunct,ncj) / Zij
         alf2 = alf * 2.D0
         ti   = a(ifunct,nci) / Zij
         tj   = a(jfunct,ncj) / Zij
         Q(1) = ti*r(Nuc(ifunct),1) + tj*r(Nuc(jfunct),1)
         Q(2) = ti*r(Nuc(ifunct),2) + tj*r(Nuc(jfunct),2)
         Q(3) = ti*r(Nuc(ifunct),3) + tj*r(Nuc(jfunct),3)
         ccoef= c(ifunct,nci) * c(jfunct,ncj)

         ss    = pi32 * exp(-alf*dd) / (Zij*sqrt(Zij))
         ovlap = ss
         sks   = alf * (3.D0 - alf2*dd) * ovlap
         tn    = sks

         ! Loops over nuclei, nuclear attraction matrix elements.
         ! tna accumulates nuclear attraction over all nuclei.
         tna   = 0.D0
         temp0 = 2.D0 * sqrt(Zij/pi) * ss
         do iatom = 1, my_natom
            q1 = Q(1) - r(iatom,1)
            q2 = Q(2) - r(iatom,2)
            q3 = Q(3) - r(iatom,3)
            uf = (q1*q1 + q2*q2 + q3*q3) * Zij

            temp = - Iz(iatom) * temp0
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)

            temp = 2.0D0 * Zij * s1s(iatom)
            x0x(iatom,1) = temp * q1
            x0x(iatom,2) = temp * q2
            x0x(iatom,3) = temp * q3
            tna          = tna + s0s(iatom)
         enddo

         ! L2: different p in the p shell.
         t3 = alf2 * ss
         te = ccoef * rho(ifunct + ((M2-jfunct)*(jfunct-1))/2)
         t4 = te * 2.D0 * a(ifunct, nci)
         t5 = te * 2.D0 * a(jfunct, ncj)

         do l2 = 1, 3
            t1   = Q(l2) - r(Nuc(ifunct),l2)
            piks = t1 * (sks + alf2 * ss)
            tx   = r(Nuc(ifunct),l2) - r(Nuc(jfunct),l2)

            ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) + t4 * piks
            ff(Nuc(jfunct),l2) = ff(Nuc(jfunct),l2) + t5 *(piks + tx*(sks + t3))

            ! Loop over nuclei, specific part
            do iatom = 1, my_natom
               piNs  = t1 * s0s(iatom) - (Q(l2)-r(iatom,l2)) * s1s(iatom)
               ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) + t4 * piNs
               ff(Nuc(jfunct),l2) = ff(Nuc(jfunct),l2) + &
                                    t5 * (piNs + tx * s0s(iatom))
               ff(iatom,l2)       = ff(iatom,l2) + te * x0x(iatom,l2)
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo
   !print*, "s|s", ff

   ! Second loop - (p|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1, ns
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci=1, ncont(ifunct)
      do ncj=1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         Z2   = 2.D0 * Zij
         ti   = a(ifunct,nci) / Zij
         tj   = a(jfunct,ncj) / Zij
         Q(1) = ti*r(Nuc(ifunct),1) + tj*r(Nuc(jfunct),1)
         Q(2) = ti*r(Nuc(ifunct),2) + tj*r(Nuc(jfunct),2)
         Q(3) = ti*r(Nuc(ifunct),3) + tj*r(Nuc(jfunct),3)

         alf   = ti   * a(jfunct,ncj)
         alf2  = 2.D0 * alf
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

         ss   = pi32 * exp(-alf*dd) / (Zij*sqrt(Zij))
         sks  = alf  * (3.D0 - alf2*dd) * ss
         t10 = ss  / Z2
         t15 = sks / Z2
         t20 = t15 - alf * ss / a(ifunct,nci)

         ! Loop over nuclei, common for all shells
         temp0 = 2.D0 * sqrt(Zij/pi) * ss
         do iatom = 1, my_natom
            q1 = Q(1) - r(iatom,1)
            q2 = Q(2) - r(iatom,2)
            q3 = Q(3) - r(iatom,3)
            uf = (q1*q1 + q2*q2 +q3*q3) * Zij

            temp = - temp0 * Iz(iatom)
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)

            temp = Z2 * s1s(iatom)
            x0x(iatom,1) = temp * q1
            x0x(iatom,2) = temp * q2
            x0x(iatom,3) = temp * q3

            temp = Z2 * s2s(iatom)
            x1x(iatom,1) = temp * q1
            x1x(iatom,2) = temp * q2
            x1x(iatom,3) = temp * q3
         enddo

         do l1 = 1, 3
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            ps  = ss  * t1
            pks = sks *t1 + alf2 * ps

            te = ccoef * rho(ifunct + l1 - 1 + ((M2-jfunct)*(jfunct-1))/2)
            t4 = te * 2.D0 * a(ifunct,nci)
            t5 = te * 2.D0 * a(jfunct,ncj)

            do l2 = 1, 3
               t1  = Q(l2) - r(Nuc(ifunct),l2)
               t2  = Q(l2) - r(Nuc(jfunct),l2)
               ds  = t1 * ps
               dks = t1 * pks
               pp  = t2 * ps
               pkp = t2 * pks

               if (l1 .eq. l2) then
                  ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) - te * sks
                  ds  = ds  + t10
                  dks = dks + t20
                  pp  = pp  + t10
                  pkp = pkp + t15
               endif

               dks = dks + alf2 * ds
               pkp = pkp + alf2 * pp

               ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) + t4 * dks
               ff(Nuc(jfunct),l2) = ff(Nuc(jfunct),l2) + t5 * pkp

            enddo
         enddo

         ! Nuclear attraction part
         do iatom = 1, my_natom
            t50 = (s0s(iatom) - s1s(iatom)) / Z2

            do  l1 = 1, 3
               t1  = Q(l1) - r(Nuc(ifunct),l1)
               t2  = Q(l1) - r(iatom,l1)
               p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
               p1s = t1 * s1s(iatom) - t2 * s2s(iatom)

               dn(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
               dn(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
               dn(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
               dn(l1) = dn(l1) + s1s(iatom)

               te = rho(ifunct + l1 - 1 + ((M2-jfunct) * (jfunct-1)) /2) * ccoef
               t4 = te * 2.D0 * a(ifunct,nci)
               t5 = te * 2.D0 * a(jfunct,ncj)

               do l2 = 1, 3
                  dNs = (Q(l2) - r(Nuc(ifunct),l2)) * p0s - &
                        (Q(l2) - r(iatom,l2)      ) * p1s

                  if (l1 .eq. l2) then
                     dNs = dNs + t50
                     ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) - te*s0s(iatom)
                  endif

                  pNp = dNs + (r(Nuc(ifunct),l2) - r(Nuc(jfunct),l2))*p0s
                  ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) + t4 * dNs
                  ff(Nuc(jfunct),l2) = ff(Nuc(jfunct),l2) + t5 * pNp
                  ff(iatom,l2)       = ff(iatom,l2)       + te * dn(l2)
               enddo
            enddo
         enddo
         ! End of nuclear attraction
      enddo
      enddo
   enddo
   enddo
   !print*, "p|s", ff

   ! Third loop - (p|p)
   do ifunct = ns+1, ns+np , 3
   do jfunct = ns+1, ifunct, 3
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         Z2   = 2.D0*Zij
         ti   = a(ifunct,nci) / Zij
         tj   = a(jfunct,ncj) / Zij
         Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)

         alf   = ti   * a(jfunct,ncj)
         alf2  = 2.D0 * alf
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

         ss  = pi32 * exp(-alf*dd) / (Zij*sqrt(Zij))
         sks = alf  * (3.D0 - alf2*dd) * ss
         t10 = ss  / Z2
         t20 = sks / Z2

         ! Loops over nuclei, common for all shells
         temp0 = 2.D0 * sqrt(Zij/pi) * ss
         do iatom = 1, my_natom
            q1 = Q(1) - r(iatom,1)
            q2 = Q(2) - r(iatom,2)
            q3 = Q(3) - r(iatom,3)
            uf = (q1*q1 + q2*q2 + q3*q3) * Zij

            temp = -temp0 * Iz(iatom)
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
            s3s(iatom) = temp * FUNCT(3,uf)

            temp = Z2 * s1s(iatom)
            x0x(iatom,1) = temp * q1
            x0x(iatom,2) = temp * q2
            x0x(iatom,3) = temp * q3

            temp = Z2 * s2s(iatom)
            x1x(iatom,1) = temp * q1
            x1x(iatom,2) = temp * q2
            x1x(iatom,3) = temp * q3

            temp = Z2 * s3s(iatom)
            x2x(iatom,1) = temp * q1
            x2x(iatom,2) = temp * q2
            x2x(iatom,3) = temp * q3
         enddo

         do l1 = 1, 3
            t1   = Q(l1) - r(Nuc(ifunct),l1)
            pis  = ss  * t1
            piks = sks * t1 + alf2 * pis
            t11  = pis  / Z2
            t12  = piks / Z2
            t16  = t12 - pis * alf / a(jfunct,ncj)

            lij = 3
            if (ifunct .eq. jfunct) then
               lij = l1
            endif
            do l2 = 1, lij
               t1   = Q(l2) - r(Nuc(ifunct),l2)
               t2   = Q(l2) - r(Nuc(jfunct),l2)
               spj  = t2 * ss
               skpj = sks * t2 + alf2 * spj
               t13  = spj  / Z2
               t15  = skpj / Z2
               t14  = t15 - alf * spj / a(ifunct,nci)

               pp   = t2*pis
               pkp  = t2*piks
               if (l1 .eq. l2) then
                  pp  = pp  + t10
                  pkp = pkp + t20
               endif
               pkp = pkp + alf2*pp

               k_ind = ifunct+l1-1 + ((M2-(jfunct+l2-1))*(jfunct+l2-2))/2
               te = ccoef*rho(k_ind)
               t5 = te * 2.D0 * a(jfunct,ncj)
               t4 = te * 2.D0 * a(ifunct,nci)

               do l3 = 1, 3
                  t1 = Q(l3) - r(Nuc(ifunct),l3)
                  t2 = Q(l3) - r(Nuc(jfunct),l3)

                  dp  = t1 * pp
                  dkp = t1 * pkp
                  pkd = t2 * pkp
                  if (l1 .eq. l3) then
                     dp  = dp  + t13
                     dkp = dkp + t14
                     pkd = pkd + t15

                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * skpj
                  endif

                  if (l2.eq.l3) then
                     dp  = dp  + t11
                     dkp = dkp + t12
                     pkd = pkd + t16

                     ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) - te * piks
                  endif
                  dkp = dkp + alf2 * dp
                  pkd = pkd + alf2 * (dp + (r(Nuc(ifunct),l3) -&
                                            r(Nuc(jfunct),l3)) * pp)
                  ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * dkp
                  ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * pkd
               enddo
            enddo
         enddo

         ! Nuclear attraction part
         do iatom = 1, my_natom
            t15 = (s0s(iatom)   - s1s(iatom))   / Z2
            t25 = (s1s(iatom)   - s2s(iatom))   / Z2
            t26 = (x0x(iatom,1) - x1x(iatom,1)) / Z2
            t27 = (x0x(iatom,2) - x1x(iatom,2)) / Z2
            t28 = (x0x(iatom,3) - x1x(iatom,3)) / Z2

            do l1 = 1, 3
               t1  = Q(l1) - r(Nuc(ifunct),l1)
               t2  = Q(l1) - r(iatom,l1)
               p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
               p1s = t1 * s1s(iatom) - t2 * s2s(iatom)
               t30 = (p0s - p1s) / Z2
               p2s = t1 * s2s(iatom) - t2 * s3s(iatom)

               ! dn(u) (pi|Au|s)
               dn(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
               dn(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
               dn(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
               dn(l1) = dn(l1) + s1s(iatom)

               dn1(1)  = t1 * x1x(iatom,1) - t2 * x2x(iatom,1)
               dn1(2)  = t1 * x1x(iatom,2) - t2 * x2x(iatom,2)
               dn1(3)  = t1 * x1x(iatom,3) - t2 * x2x(iatom,3)
               dn1(l1) = dn1(l1) + s2s(iatom)

               lij = 3
               if (ifunct .eq. jfunct) then
                  lij = l1
               endif

               do l2 = 1, lij
                  t1   = Q(l2) - r(Nuc(ifunct),l2)
                  t2   = Q(l2) - r(iatom,l2)
                  t3   = Q(l2) - r(Nuc(jfunct),l2)
                  s0p  = t3 * s0s(iatom) - t2 * s1s(iatom)
                  s1p  = t3 * s1s(iatom) - t2 * s2s(iatom)
                  t29  = (s0p - s1p) / Z2

                  dn2(1)  = t3 * dn(1) - t2 * dn1(1)
                  dn2(2)  = t3 * dn(2) - t2 * dn1(2)
                  dn2(3)  = t3 * dn(3) - t2 * dn1(3)
                  dn2(l2) = dn2(l2) + p1s

                  d0s = t1 * p0s - t2 * p1s
                  d1s = t1 * p1s - t2 * p2s

                  pNp  = t3 * p0s - t2 * p1s
                  pN1p = t3 * p1s - t2 * p2s
                  if (l1 .eq. l2) then
                     pNp  = pNp  + t15
                     pN1p = pN1p + t25
                     dn2(1) = dn2(1) + t26
                     dn2(2) = dn2(2) + t27
                     dn2(3) = dn2(3) + t28
                  endif

                  k_ind = ifunct + l1-1 + ((M2-(jfunct+l2-1))*(jfunct+l2-2))/2
                  te = ccoef * rho(k_ind)
                  t4 = te * 2.D0 * a(ifunct,nci)
                  t5 = te * 2.D0 * a(jfunct,ncj)

                  ! Gradients
                  do l3 = 1, 3
                    t1  = Q(l3) - r(Nuc(ifunct),l3)
                    t2  = Q(l3) - r(iatom,l3)
                    dNp = t1 * pNp - t2 * pN1p

                    if (l1.eq.l3) then
                       dNp = dNp + t29
                       ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * s0p
                    endif

                    if (l2.eq.l3) then
                       dNp = dNp + t30
                       ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) - te * p0s
                    endif

                    pNd = dNp + pNp * (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3))
                    ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * dNp
                    ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * pNd
                    ff(iatom,l3)       = ff(iatom,l3)       + te * dn2(l3)
                 enddo
              enddo
           enddo
        enddo
        ! End of nuclear attraction
      enddo
      enddo
   enddo
   enddo
   !print*, "p|p", ff

   ! Fourth loop - (d|s)
   do ifunct = ns+np+1, M, 6
   do jfunct = 1, ns
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci=1, ncont(ifunct)
      do ncj=1, ncont(jfunct)
         Zij = a(ifunct,nci) + a(jfunct,ncj)
         Z2  = 2.D0 * Zij

         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij
         alf  = a(ifunct,nci) * a(jfunct,ncj) / Zij
         alf2 = 2.D0 * alf
         alf3 = alf / a(ifunct,nci)

         ss  = pi32 * exp(-alf*dd) / (Zij*sqrt(Zij))
         sks = alf * (3.D0 - alf2*dd) * ss
         t10 = ss / Z2
         t11 = sks / Z2 - alf3 * ss
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

         ! Loop over nuclei, common for all shells
         temp0 = 2.D0 * sqrt(Zij/pi) * ss
         do iatom = 1, my_natom
            q1 = Q(1) - r(iatom,1)
            q2 = Q(2) - r(iatom,2)
            q3 = Q(3) - r(iatom,3)
            uf = (q1*q1 + q2*q2 + q3*q3) * Zij

            temp = -temp0 * Iz(iatom)
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
            s3s(iatom) = temp * FUNCT(3,uf)

            temp = Z2 * s1s(iatom)
            x0x(iatom,1) = temp * q1
            x0x(iatom,2) = temp * q2
            x0x(iatom,3) = temp * q3

            temp = Z2 * s2s(iatom)
            x1x(iatom,1) = temp * q1
            x1x(iatom,2) = temp * q2
            x1x(iatom,3) = temp * q3

            temp = Z2 * s3s(iatom)
            x2x(iatom,1) = temp * q1
            x2x(iatom,2) = temp * q2
            x2x(iatom,3) = temp * q3
         enddo

         ! The ordering of the d shell goes like: xx,yx,yy,zx,zy,zz
         ! ( 11, 21, 22, 31, 32, 33 )
         do l1 = 1, 3
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            ps  = ss  * t1
            pks = sks * t1 + alf2 * ps
            t12 = ps  / Z2
            t13 = pks / Z2
            t17 = t13 - alf3 * ps

            do l2 = 1, l1
               t1   = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = ss * t1
               pjks = t1 * sks + alf2 * pjs
               t14  = pjs  / Z2
               t15  = pjks / Z2
               t16  = t15 - alf3 * pjs
               tn   = t1 * pks
               ovlap= t1 * ps

               f1 = 1.D0
               if (l1 .eq. l2) then
                  ovlap = ovlap + t10
                  tn    = tn + t11
                  f1    = sq3
               endif
               tn  = tn + alf2 * ovlap

               ! Gradients
               l12   = l1 * (l1-1) / 2 + l2
               k_ind = ifunct + l12-1 + ((M2-jfunct)*(jfunct-1))/2

               te = ccoef * rho(k_ind) / f1
               t4 = te * 2.D0 * a(ifunct,nci)
               t5 = te * 2.D0 * a(jfunct,ncj)

               do l3 = 1, 3
                  t1  = Q(l3) - r(Nuc(jfunct),l3)
                  t2  = Q(l3) - r(Nuc(ifunct),l3)

                  dp  = t1 * ovlap
                  dkp = t1 * tn
                  fks = t2 * tn
                  if (l1.eq.l3) then
                     dp  = dp  + t14
                     dkp = dkp + t15
                     fks = fks + t16
                     ff(Nuc(ifunct),l1) = ff(Nuc(ifunct),l1) - te * pjks
                  endif
                  if (l2.eq.l3) then
                     dp  = dp  + t12
                     dkp = dkp + t13
                     fks = fks + t17
                     ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) - te * pks
                  endif
                  dkp = dkp + alf2 * dp
                  fks = fks + alf2 * (dp - (r(Nuc(ifunct),l3) - &
                                            r(Nuc(jfunct),l3)) * ovlap)

                  ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * fks
                  ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * dkp
               enddo
            enddo
         enddo

         ! Nuclear attraction part
         do iatom = 1, my_natom
            t7  = (s0s(iatom) - s1s(iatom)) / Z2
            t8  = (s1s(iatom) - s2s(iatom)) / Z2
            t26 = (x0x(iatom,1) - x1x(iatom,1)) / Z2
            t27 = (x0x(iatom,2) - x1x(iatom,2)) / Z2
            t28 = (x0x(iatom,3) - x1x(iatom,3)) / Z2

            do l1 = 1, 3
               t1  = Q(l1) - r(Nuc(ifunct),l1)
               t2  = Q(l1) - r(iatom,l1)
               p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
               p1s = t1 * s1s(iatom) - t2 * s2s(iatom)
               p2s = t1 * s2s(iatom) - t2 * s3s(iatom)
               t30 = (p0s - p1s) / Z2

               ! dn(u) (pi|Au|s)
               dn(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
               dn(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
               dn(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
               dn(l1) = dn(l1) + s1s(iatom)

               dn1(1)  = t1 * x1x(iatom,1) - t2 * x2x(iatom,1)
               dn1(2)  = t1 * x1x(iatom,2) - t2 * x2x(iatom,2)
               dn1(3)  = t1 * x1x(iatom,3) - t2 * x2x(iatom,3)
               dn1(l1) = dn1(l1) + s2s(iatom)

               do l2 = 1, l1
                  t1   = Q(l2) - r(Nuc(ifunct),l2)
                  t2   = Q(l2) - r(iatom,l2)
                  tna  = t1 * p0s - t2 * p1s
                  tn1a = t1 * p1s - t2 * p2s
                  pj0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  pj1s = t1 * s1s(iatom) - t2 * s2s(iatom)
                  t29  = (pj0s - pj1s) / Z2

                  dn2(1)  = t1 * dn(1) - t2 * dn1(1)
                  dn2(2)  = t1 * dn(2) - t2 * dn1(2)
                  dn2(3)  = t1 * dn(3) - t2 * dn1(3)
                  dn2(l2) = dn2(l2) + p1s

                  f1 = 1.D0
                  if (l1 .eq. l2) then
                     tna    = tna  + t7
                     tn1a   = tn1a + t8
                     f1     = sq3
                     dn2(1) = dn2(1) + t26
                     dn2(2) = dn2(2) + t27
                     dn2(3) = dn2(3) + t28
                  endif
                  dNs  = tna
                  dN1s = tn1a

                  ! The ordering of the d shell goes like: xx,yx,yy,zx,zy,zz
                  ! ( 11, 21, 22, 31, 32, 33 )
                  ! Gradients
                  l12   = l1 * (l1-1)/2 +l2
                  k_ind = ifunct + l12-1 + ((M2-jfunct)*(jfunct-1))/2
                  te    = ccoef * rho(k_ind) / f1
                  t4    = te * 2.D0 * a(ifunct,nci)
                  t5    = te * 2.D0 * a(jfunct,ncj)

                  do l3 = 1, 3
                     t1 = Q(l3) - r(Nuc(jfunct),l3)
                     t2 = Q(l3) - r(iatom,l3)

                     dNp = t1 * dNs - t2 * dN1s
                     if (l1.eq.l3) then
                        dNp = dNp + t29
                        ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * pj0s
                     endif
                     if (l2.eq.l3) then
                        dNp = dNp + t30
                        ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * p0s
                     endif
                     fNs = dNp - (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) * dNs

                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * fNs
                     ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * dNp
                     ff(iatom,l3)       = ff(iatom,l3)       + te * dn2(l3)
                  enddo
               enddo
            enddo
         enddo
         ! End of nuclear attraction
      enddo
      enddo
   enddo
   enddo
   !print*, "d|s", ff

   ! Fifth loop - (d|p)
   do ifunct = ns+np+1, M, 6
   do jfunct = ns+1, ns+np, 3
      dd = d(Nuc(ifunct),Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)

         Zij = a(ifunct,nci) + a(jfunct,ncj)
         Z2  = 2.D0 * Zij
         Q(1) = (a(ifunct,nci)*r(Nuc(ifunct),1) + &
                 a(jfunct,ncj)*r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci)*r(Nuc(ifunct),2) + &
                 a(jfunct,ncj)*r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci)*r(Nuc(ifunct),3) + &
                 a(jfunct,ncj)*r(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct,nci) * a(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         alf3  = alf / a(ifunct,nci)
         alf4  = alf / a(jfunct,ncj)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

         ss  = pi32 * exp(-alf*dd) / (Zij*sqrt(Zij))
         sks = alf * (3.D0 - alf2*dd) * ss
         t10 = ss  / Z2
         t30 = sks / Z2
         t11 = t30 - alf3 * ss

         ! Loops over nuclei, common for all shells
         temp0 = 2.D0 * sqrt(Zij/pi) * ss
         do iatom = 1, my_natom
            q1 = Q(1) - r(iatom,1)
            q2 = Q(2) - r(iatom,2)
            q3 = Q(3) - r(iatom,3)
            uf = (q1*q1 + q2*q2 + q3*q3) * Zij

            temp = -temp0 * Iz(iatom)
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
            s3s(iatom) = temp * FUNCT(3,uf)
            s4s(iatom) = temp * FUNCT(4,uf)

            temp = Z2 * s1s(iatom)
            x0x(iatom,1) = temp * q1
            x0x(iatom,2) = temp * q2
            x0x(iatom,3) = temp * q3

            temp = Z2 * s2s(iatom)
            x1x(iatom,1) = temp * q1
            x1x(iatom,2) = temp * q2
            x1x(iatom,3) = temp * q3

            temp = Z2 * s3s(iatom)
            x2x(iatom,1) = temp * q1
            x2x(iatom,2) = temp * q2
            x2x(iatom,3) = temp * q3

            temp = Z2 * s4s(iatom)
            x3x(iatom,1) = temp * q1
            x3x(iatom,2) = temp * q2
            x3x(iatom,3) = temp * q3
         enddo

         do l1 = 1, 3
            t1   = Q(l1) - r(Nuc(ifunct),l1)
            pis  = ss * t1
            piks = sks * t1 + alf2 * pis
            t14  = pis / Z2
            t15  = piks / Z2

            do l2 = 1, l1
               t1    = Q(l2) - r(Nuc(ifunct),l2)
               pjs   = ss * t1
               pjks  = sks * t1 + alf2 * pjs
               dijs  = t1 * pis
               dijks = t1 * piks
               f1    = 1.D0

               if (l1 .eq. l2) then
                  f1    = sq3
                  dijs  = dijs  + t10
                  dijks = dijks + t11
               endif

               dijks = dijks + alf2 * dijs
               t12   = pjs   / Z2
               t13   = pjks  / Z2
               t22   = dijs  / Z2
               t23   = dijks / Z2
               t24   = t23 - alf4 * dijs

               do l3 = 1, 3
                  t1    = Q(l3) - r(Nuc(jfunct),l3)
                  ovlap = t1 * dijs
                  tn    = t1 * dijks
                  pipk  = t1 * pis
                  pikpk = t1 * piks
                  pjpk  = t1 * pjs
                  pjkpk = t1 * pjks

                  if (l1.eq.l3) then
                     pipk  = pipk  +  t10
                     pikpk = pikpk + t30
                     ovlap = ovlap + t12
                     tn    = tn + t13
                  endif

                  if (l2.eq.l3) then
                     pjpk  = pjpk  + t10
                     pjkpk = pjkpk + t30
                     ovlap = ovlap + t14
                     tn    = tn + t15
                  endif

                  pjkpk = pjkpk + alf2 * pjpk
                  pikpk = pikpk + alf2 * pipk

                  t16 = pjpk  / Z2
                  t17 = pjkpk / Z2
                  t18 = t17 - alf3 * pjpk
                  t19 = pipk  / Z2
                  t20 = pikpk / Z2
                  t21 = t20 - alf3 * pipk
                  tn  = tn  + alf2 * ovlap

                  ! Gradients
                  l12   = l1 * (l1-1) / 2 + l2
                  k_ind = ifunct + l12-1 + ((M2-(jfunct+l3-1))*(jfunct+l3-2))/2

                  te = rho(k_ind) * ccoef / f1
                  t4 = te * 2.D0 * a(ifunct,nci)
                  t5 = te * 2.D0 * a(jfunct,ncj)

                  do l4 = 1, 3
                     t1 = Q(l4) - r(Nuc(jfunct),l4)
                     t2 = Q(l4) - r(Nuc(ifunct),l4)

                     dsd = t1 * ovlap
                     dkd = t1 * tn
                     fkp = t2 * tn
                     if (l1.eq.l4) then
                        dsd = dsd + t16
                        dkd = dkd + t17
                        fkp = fkp + t18
                        ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) - te * pjkpk
                     endif
                     if (l2.eq.l4) then
                        dsd = dsd + t19
                        dkd = dkd + t20
                        fkp = fkp + t21
                        ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) - te * pikpk
                     endif
                     if (l3.eq.l4) then
                        dsd = dsd + t22
                        dkd = dkd + t24
                        fkp = fkp + t23
                        ff(Nuc(jfunct),l4) = ff(Nuc(jfunct),l4) - te * dijks
                     endif
                     dkd = dkd + alf2 * dsd
                     fkp = fkp + alf2 * (dsd - (r(Nuc(ifunct),l4) - &
                                                r(Nuc(jfunct),l4)) * ovlap)
                     ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) + t4 * fkp
                     ff(Nuc(jfunct),l4) = ff(Nuc(jfunct),l4) + t5 * dkd
                  enddo
               enddo
            enddo
         enddo

         ! Nuclear attraction
         do iatom = 1, my_natom
            t7  = (s0s(iatom) - s1s(iatom)) / Z2
            t8  = (s1s(iatom) - s2s(iatom)) / Z2
            t9  = (s2s(iatom) - s3s(iatom)) / Z2
            t26 = (x0x(iatom,1) - x1x(iatom,1)) / Z2
            t27 = (x0x(iatom,2) - x1x(iatom,2)) / Z2
            t28 = (x0x(iatom,3) - x1x(iatom,3)) / Z2
            t29 = (x1x(iatom,1) - x2x(iatom,1)) / Z2
            t30 = (x1x(iatom,2) - x2x(iatom,2)) / Z2
            t31 = (x1x(iatom,3) - x2x(iatom,3)) / Z2

            do l1 = 1, 3
               t1  = Q(l1) - r(Nuc(ifunct),l1)
               t2  = Q(l1) - r(iatom,l1)
               p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
               p1s = t1 * s1s(iatom) - t2 * s2s(iatom)
               p2s = t1 * s2s(iatom) - t2 * s3s(iatom)
               p3s = t1 * s3s(iatom) - t2 * s4s(iatom)

               ! dn(u) (pi|Au|s)
               dn(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
               dn(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
               dn(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
               dn(l1) = dn(l1) + s1s(iatom)

               dn1(1)  = t1 * x1x(iatom,1) - t2 * x2x(iatom,1)
               dn1(2)  = t1 * x1x(iatom,2) - t2 * x2x(iatom,2)
               dn1(3)  = t1 * x1x(iatom,3) - t2 * x2x(iatom,3)
               dn1(l1) = dn1(l1) + s2s(iatom)

               t51 = (dn(1) - dn1(1)) / Z2
               t52 = (dn(2) - dn1(2)) / Z2
               t53 = (dn(3) - dn1(3)) / Z2

               dn2(1)  = t1 * x2x(iatom,1) - t2 * x3x(iatom,1)
               dn2(2)  = t1 * x2x(iatom,2) - t2 * x3x(iatom,2)
               dn2(3)  = t1 * x2x(iatom,3) - t2 * x3x(iatom,3)
               dn2(l1) = dn2(l1) + s3s(iatom)

               do l2 = 1, l1
                  t1   = Q(l2) - r(Nuc(ifunct),l2)
                  t2   = Q(l2) - r(iatom,l2)
                  pj0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  pj1s = t1 * s1s(iatom) - t2 * s2s(iatom)
                  pj2s = t1 * s2s(iatom) - t2 * s3s(iatom)

                  ! dn3 (dij || s) order 0
                  dn3(1)  = t1 * dn(1) - t2 * dn1(1)
                  dn3(2)  = t1 * dn(2) - t2 * dn1(2)
                  dn3(3)  = t1 * dn(3) - t2 * dn1(3)
                  dn3(l2) = dn3(l2) + p1s

                  ! dn4 (dij || s) order 1
                  dn4(1)  = t1 * dn1(1) - t2 * dn2(1)
                  dn4(2)  = t1 * dn1(2) - t2 * dn2(2)
                  dn4(3)  = t1 * dn1(3) - t2 * dn2(3)
                  dn4(l2) = dn4(l2) + p2s

                  ! dn6 and dn7 used for (pj | s) both orders 0 and 1
                  dn6(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
                  dn6(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
                  dn6(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
                  dn6(l2) = dn6(l2) + s1s(iatom)

                  dn7(1)  = t1 * x1x(iatom,1) - t2 * x2x(iatom,1)
                  dn7(2)  = t1 * x1x(iatom,2) - t2 * x2x(iatom,2)
                  dn7(3)  = t1 * x1x(iatom,3) - t2 * x2x(iatom,3)
                  dn7(l2) = dn7(l2) + s2s(iatom)

                  t54 = (dn6(1) - dn7(1)) / Z2
                  t55 = (dn6(2) - dn7(2)) / Z2
                  t56 = (dn6(3) - dn7(3)) / Z2
                  f1  = 1.D0
                  d0s = t1 * p0s - t2 * p1s
                  d1s = t1 * p1s - t2 * p2s
                  d2s = t1 * p2s - t2 * p3s

                  if (l1 .eq. l2) then
                     f1  = sq3
                     d0s = d0s + t7
                     d1s = d1s + t8
                     d2s = d2s + t9

                     dn3(1) = dn3(1) + t26
                     dn3(2) = dn3(2) + t27
                     dn3(3) = dn3(3) + t28
                     dn4(1) = dn4(1) + t29
                     dn4(2) = dn4(2) + t30
                     dn4(3) = dn4(3) + t31
                  endif

                  do l3 = 1, 3
                     ! dn5 (dij || Pk ) order 0, the one needed for derivatives
                     t1   = Q(l3) - r(Nuc(jfunct),l3)
                     t2   = Q(l3) - r(iatom,l3)
                     tna  = t1 * d0s - t2 * d1s
                     tn1a = t1 * d1s - t2 * d2s

                     pi0p = t1 * p0s - t2 * p1s
                     pi1p = t1 * p1s - t2 * p2s
                     pj0p = t1 * pj0s - t2 * pj1s
                     pj1p = t1 * pj1s - t2 * pj2s

                     dn5(1)  = t1 * dn3(1) - t2 * dn4(1)
                     dn5(2)  = t1 * dn3(2) - t2 * dn4(2)
                     dn5(3)  = t1 * dn3(3) - t2 * dn4(3)
                     dn5(l3) = dn5(l3) + d1s

                     if (l1.eq.l3) then
                        tna  = tna  + (pj0s - pj1s) / Z2
                        tn1a = tn1a + (pj1s - pj2s) / Z2
                        pi0p = pi0p + t7
                        pi1p = pi1p + t8
                        dn5(1) = dn5(1) + t54
                        dn5(2) = dn5(2) + t55
                        dn5(3) = dn5(3) + t56
                     endif

                     if (l2.eq.l3) then
                        tna  = tna  + (p0s - p1s) / Z2
                        tn1a = tn1a + (p1s - p2s) / Z2
                        pj0p = pj0p + t7
                        pj1p = pj1p + t8
                        dn5(1) = dn5(1) + t51
                        dn5(2) = dn5(2) + t52
                        dn5(3) = dn5(3) + t53
                     endif

                     l12   = l1 * (l1-1) / 2 + l2
                     k_ind = ifunct +l12-1+ ((M2-(jfunct+l3-1))*(jfunct+l3-2))/2

                     te = rho(k_ind) * ccoef / f1
                     t4 = te * 2.D0 * a(ifunct,nci)
                     t5 = te * 2.D0 * a(jfunct,ncj)

                     ! Gradients
                     do l4 = 1, 3
                        t1 = Q(l4) - r(Nuc(jfunct),l4)
                        t2 = Q(l4) - r(iatom,l4)

                        dNd = t1 * tna - t2 * tn1a
                        if (l1 .eq. l4) then
                           ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) - te * pj0p
                           dNd = dNd + (pj0p - pj1p) / Z2
                        endif

                        if (l2 .eq. l4) then
                           ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) - te * pi0p
                           dNd = dNd + (pi0p - pi1p) / Z2
                        endif

                        if (l3 .eq. l4) then
                           ff(Nuc(jfunct),l4) = ff(Nuc(jfunct),l4) - te * d0s
                           dNd = dNd + (d0s - d1s) / Z2
                        endif
                        fNp = dNd - (r(Nuc(ifunct),l4) - r(Nuc(jfunct),l4)) *tna

                        ff(Nuc(ifunct),l4) = ff(Nuc(ifunct),l4) + t4 * fNp
                        ff(Nuc(jfunct),l4) = ff(Nuc(jfunct),l4) + t5 * dNd
                        ff(iatom,l4)       = ff(iatom,l4)       + te * dn5(l4)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         ! End of nuclear attraction
      enddo
      enddo
   enddo
   enddo
   !print*, "d|p", ff

   ! Sixth and final loop - (d|d)
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+np+1, ifunct, 6
      dd = d(Nuc(ifunct),Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij = a(ifunct,nci) + a(jfunct,ncj)
         Z2  = 2.D0 * Zij
         Q(1)= (a(ifunct,nci)*r(Nuc(ifunct),1) + &
                a(jfunct,ncj)*r(Nuc(jfunct),1)) / Zij
         Q(2)= (a(ifunct,nci)*r(Nuc(ifunct),2) + &
                a(jfunct,ncj)*r(Nuc(jfunct),2)) / Zij
         Q(3)= (a(ifunct,nci)*r(Nuc(ifunct),3) + &
                a(jfunct,ncj)*r(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct,nci) * a(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         alf3  = alf / a(ifunct,nci)
         alf4  = alf / a(jfunct,ncj)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

         ss  = pi32 * exp(-alf*dd) / (Zij*sqrt(Zij))
         sks = alf * (3.D0 - alf2*dd) * ss
         t0  = ss  / Z2
         t12 = sks / Z2
         t10 = t12 - alf3 * ss

         ! Loops over nuclei, common for all shells
         temp0 = 2.D0 * sqrt(Zij/pi) * ss
         do iatom = 1, my_natom
            q1 = Q(1) - r(iatom,1)
            q2 = Q(2) - r(iatom,2)
            q3 = Q(3) - r(iatom,3)
            uf = (q1*q1 + q2*q2 + q3*q3) * Zij

            temp = -temp0 * Iz(iatom)
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
            s3s(iatom) = temp * FUNCT(3,uf)
            s4s(iatom) = temp * FUNCT(4,uf)
            s5s(iatom) = temp * FUNCT(5,uf)
            s6s(iatom) = temp * FUNCT(6,uf)

            temp = Z2 * s1s(iatom)
            x0x(iatom,1) = temp * q1
            x0x(iatom,2) = temp * q2
            x0x(iatom,3) = temp * q3

            temp = Z2 * s2s(iatom)
            x1x(iatom,1) = temp * q1
            x1x(iatom,2) = temp * q2
            x1x(iatom,3) = temp * q3

            temp = Z2 * s3s(iatom)
            x2x(iatom,1) = temp * q1
            x2x(iatom,2) = temp * q2
            x2x(iatom,3) = temp * q3

            temp = Z2 * s4s(iatom)
            x3x(iatom,1) = temp * q1
            x3x(iatom,2) = temp * q2
            x3x(iatom,3) = temp * q3

            temp = Z2 * s5s(iatom)
            x4x(iatom,1) = temp * q1
            x4x(iatom,2) = temp * q2
            x4x(iatom,3) = temp * q3
         enddo

         do l1 = 1, 3
            t1   = Q(l1) - r(Nuc(ifunct),l1)
            pis  = ss  * t1
            piks = sks * t1 + alf2 * pis
            t14  = pis  / Z2
            t15  = piks / Z2
            t26  = t15 - alf4 * pis

            do l2 = 1, l1
               t1 = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = ss  * t1
               pjks = sks * t1 + alf2 * pjs
               t11  = pjs  / Z2
               t13  = pjks / Z2
               t25  = t13 - alf4 * pjs

               f1 = 1.D0
               t1 = Q(l2) - r(Nuc(ifunct),l2)
               dijs  = t1 * pis
               dijks = t1 * piks

               if (l1 .eq. l2) then
                  f1    = sq3
                  dijs  = dijs  + t0
                  dijks = dijks + t10
               endif
               dijks = dijks + alf2 * dijs
               t20   = dijs  / Z2
               t21   = dijks / Z2 - alf4 * dijs

               lij = 3
               if (ifunct .eq. jfunct) then
                  lij = l1
               endif

               do l3 = 1, lij
                  t2   = Q(l3)-r(Nuc(jfunct),l3)
                  spk  = t2 * ss
                  t22  = spk  / Z2
                  t23  = (t2 * sks + alf2 * spk) / Z2

                  pipk   = t2 * pis
                  pjpk   = t2 * pjs
                  pikpk  = t2 * piks
                  pjkpk  = t2 * pjks
                  dijpk  = t2 * dijs
                  dijkpk = t2 * dijks
                  if (l1.eq.l3) then
                     pipk   = pipk   + t0
                     dijpk  = dijpk  + t11
                     pikpk  = pikpk  + t12
                     dijkpk = dijkpk + t13
                  endif
                  if (l2.eq.l3) then
                     pjpk   = pjpk   + t0
                     dijpk  = dijpk  + t14
                     pjkpk  = pjkpk  + t12
                     dijkpk = dijkpk + t15
                  endif
                  pikpk  = pikpk  + alf2 * pipk
                  pjkpk  = pjkpk  + alf2 * pjpk
                  dijkpk = dijkpk + alf2 * dijpk

                  t16 = pjpk   / Z2
                  t17 = pjkpk  / Z2
                  t18 = pipk   / Z2
                  t19 = pikpk  / Z2
                  t39 = dijpk  / Z2
                  t40 = dijkpk / Z2
                  t41 = t40 - alf4 * dijpk

                  lk = l3
                  if (ifunct .eq. jfunct) then
                    lk = min(l3, Ll(l1) - Ll(l3) + l2)
                  endif
                  do l4 = 1, lk
                     f2=1.D0
                     t1 = Q(l4) - r(Nuc(jfunct),l4)
                     ovlap = t1 * dijpk
                     tn    = t1 * dijkpk

                     pjdkl  = t1 * pjpk
                     pjkdkl = t1 * pjkpk
                     pidkl  = t1 * pipk
                     pikdkl = t1 * pikpk
                     dijpl  = t1 * dijs
                     dijkpl = t1 * dijks
                     if (l1.eq.l4) then
                        ovlap  = ovlap  + t16
                        tn     = tn     + t17
                        pidkl  = pidkl  + t22
                        pikdkl = pikdkl + t23
                        dijpl  = dijpl  + t11
                        dijkpl = dijkpl + t13
                     endif
                     if (l2.eq.l4) then
                        ovlap  = ovlap  + t18
                        tn     = tn     + t19
                        pjdkl  = pjdkl  + t22
                        pjkdkl = pjkdkl + t23
                        dijpl  = dijpl  + t14
                        dijkpl = dijkpl + t15
                     endif
                     if (l3.eq.l4) then
                        ovlap  = ovlap  + t20
                        tn     = tn     + t21
                        pjdkl  = pjdkl  + t11
                        pjkdkl = pjkdkl + t25
                        pidkl  = pidkl  + t14
                        pikdkl = pikdkl + t26
                        f2     = sq3
                     endif
                     pikdkl = pikdkl + alf2 * pidkl
                     pjkdkl = pjkdkl + alf2 * pjdkl
                     dijkpl = dijkpl + alf2 * dijpl

                     ! Gradients
                     ! l12 and l34 go from 1 to 6, spanning the d shell in
                     ! the order xx, xy, yy, zx, zy, zz. The same order should
                     ! be used in ordering the basis set, before calculating
                     ! these matrix elements.
                     l12 = Ll(l1)+l2
                     l34 = Ll(l3)+l4
                     k_ind = ifunct+l12-1+((M2-(jfunct+l34-1))*(jfunct+l34-2))/2

                     te = rho(k_ind) * ccoef / (f1 * f2)
                     t5 = te * 2.D0 * a(jfunct,ncj)
                     t4 = te * 2.D0 * a(ifunct,nci)
                     tn = tn + alf2 * ovlap

                     t30 = pjdkl  / Z2
                     t31 = pjkdkl / Z2
                     t32 = t31 - alf3 * pjdkl
                     t33 = pidkl  / Z2
                     t34 = pikdkl / Z2
                     t35 = t34 - alf3 * pidkl
                     t36 = dijpl  / Z2
                     t37 = dijkpl / Z2
                     t38 = t37 - alf4 * dijpl

                     do l5 = 1, 3
                        t1 = Q(l5) - r(Nuc(ifunct),l5)
                        t2 = Q(l5) - r(Nuc(jfunct),l5)

                        fd  = t1 * ovlap
                        fkd = t1 * tn
                        dkf = t2 * tn
                        if (l1.eq.l5) then
                           fd  = fd  + t30
                           fkd = fkd + t32
                           dkf = dkf + t31
                           ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) - te * pjkdkl
                        endif
                        if (l2.eq.l5) then
                           fd  = fd + t33
                           fkd = fkd + t35
                           dkf = dkf + t34
                           ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) - te * pikdkl
                        endif
                        if (l3.eq.l5) then
                           fd  = fd + t36
                           fkd = fkd + t37
                           dkf = dkf + t38
                           ff(Nuc(jfunct),l5)=ff(Nuc(jfunct),l5) - te * dijkpl
                        endif
                        if (l4.eq.l5) then
                           fd  = fd + t39
                           fkd = fkd + t40
                           dkf = dkf + t41
                           ff(Nuc(jfunct),l5)=ff(Nuc(jfunct),l5) - te * dijkpk
                        endif
                        fkd = fkd + alf2 * fd
                        dkf = dkf + alf2 * (fd + (r(Nuc(ifunct),l5) - &
                                                  r(Nuc(jfunct),l5)) * ovlap)
                        ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) + t4 * fkd
                        ff(Nuc(jfunct),l5) = ff(Nuc(jfunct),l5) + t5 * dkf
                     enddo
                  enddo
               enddo
            enddo
         enddo

         ! Nuclear attraction
         do iatom = 1, my_natom
            t50 = (s0s(iatom) - s1s(iatom)) / Z2
            t51 = (s1s(iatom) - s2s(iatom)) / Z2
            t52 = (s2s(iatom) - s3s(iatom)) / Z2
            t53 = (s3s(iatom) - s4s(iatom)) / Z2
            t26 = (x0x(iatom,1) - x1x(iatom,1)) / Z2
            t27 = (x0x(iatom,2) - x1x(iatom,2)) / Z2
            t28 = (x0x(iatom,3) - x1x(iatom,3)) / Z2
            t29 = (x1x(iatom,1) - x2x(iatom,1)) / Z2
            t30 = (x1x(iatom,2) - x2x(iatom,2)) / Z2
            t31 = (x1x(iatom,3) - x2x(iatom,3)) / Z2
            t32 = (x2x(iatom,1) - x3x(iatom,1)) / Z2
            t33 = (x2x(iatom,2) - x3x(iatom,2)) / Z2
            t34 = (x2x(iatom,3) - x3x(iatom,3)) / Z2

            do l1 = 1, 3
               t1  = Q(l1) - r(Nuc(ifunct),l1)
               t2  = Q(l1) - r(iatom,l1)
               p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
               p1s = t1 * s1s(iatom) - t2 * s2s(iatom)
               p2s = t1 * s2s(iatom) - t2 * s3s(iatom)
               p3s = t1 * s3s(iatom) - t2 * s4s(iatom)
               p4s = t1 * s4s(iatom) - t2 * s5s(iatom)
               t54 = (p0s - p1s) / Z2
               t55 = (p1s - p2s) / Z2
               t56 = (p2s - p3s) / Z2
               t57 = (p3s - p4s) / Z2

               ! dn(u) (pi|Au|s)
               dn(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
               dn(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
               dn(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
               dn(l1) = dn(l1) + s1s(iatom)

               dn1(1)  = t1 * x1x(iatom,1) - t2 * x2x(iatom,1)
               dn1(2)  = t1 * x1x(iatom,2) - t2 * x2x(iatom,2)
               dn1(3)  = t1 * x1x(iatom,3) - t2 * x2x(iatom,3)
               dn1(l1) = dn1(l1) + s2s(iatom)

               dn2(1)  = t1 * x2x(iatom,1) - t2 * x3x(iatom,1)
               dn2(2)  = t1 * x2x(iatom,2) - t2 * x3x(iatom,2)
               dn2(3)  = t1 * x2x(iatom,3) - t2 * x3x(iatom,3)
               dn2(l1) = dn2(l1) + s3s(iatom)

               dn2b(1)  = t1 * x3x(iatom,1) - t2 * x4x(iatom,1)
               dn2b(2)  = t1 * x3x(iatom,2) - t2 * x4x(iatom,2)
               dn2b(3)  = t1 * x3x(iatom,3) - t2 * x4x(iatom,3)
               dn2b(l1) = dn2b(l1) + s4s(iatom)

               t81  = (dn(1) - dn1(1)) / Z2
               t82  = (dn(2) - dn1(2)) / Z2
               t83  = (dn(3) - dn1(3)) / Z2
               t81b = (dn1(1) - dn2(1)) / Z2
               t82b = (dn1(2) - dn2(2)) / Z2
               t83b = (dn1(3) - dn2(3)) / Z2

               do l2 = 1, l1
                  f1   = 1.D0
                  t1   = Q(l2) - r(Nuc(ifunct),l2)
                  t2   = Q(l2) - r(iatom,l2)
                  pj0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  pj1s = t1 * s1s(iatom) - t2 * s2s(iatom)
                  pj2s = t1 * s2s(iatom) - t2 * s3s(iatom)
                  pj3s = t1 * s3s(iatom) - t2 * s4s(iatom)

                  ! dn3 (dij || s) order 0
                  dn3(1)  = t1 * dn(1) - t2 * dn1(1)
                  dn3(2)  = t1 * dn(2) - t2 * dn1(2)
                  dn3(3)  = t1 * dn(3) - t2 * dn1(3)
                  dn3(l2) = dn3(l2) + p1s

                  ! dn4 (dij || s) order 1
                  dn4(1)  = t1 * dn1(1) - t2 * dn2(1)
                  dn4(2)  = t1 * dn1(2) - t2 * dn2(2)
                  dn4(3)  = t1 * dn1(3) - t2 * dn2(3)
                  dn4(l2) = dn4(l2) + p2s

                  ! dn4b (dij || s) order 2
                  dn4b(1)  = t1 * dn2(1) - t2 * dn2b(1)
                  dn4b(2)  = t1 * dn2(2) - t2 * dn2b(2)
                  dn4b(3)  = t1 * dn2(3) - t2 * dn2b(3)
                  dn4b(l2) = dn4b(l2) + p3s

                  ! dn6 and dn7 used for (pj | s) orders 0 and 1
                  dn6(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
                  dn6(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
                  dn6(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
                  dn6(l2) = dn6(l2) + s1s(iatom)

                  dn7(1)  = t1 * x1x(iatom,1) - t2 * x2x(iatom,1)
                  dn7(2)  = t1 * x1x(iatom,2) - t2 * x2x(iatom,2)
                  dn7(3)  = t1 * x1x(iatom,3) - t2 * x2x(iatom,3)
                  dn7(l2) = dn7(l2) + s2s(iatom)

                  dn7b(1)  = t1 * x2x(iatom,1) - t2 * x3x(iatom,1)
                  dn7b(2)  = t1 * x2x(iatom,2) - t2 * x3x(iatom,2)
                  dn7b(3)  = t1 * x2x(iatom,3) - t2 * x3x(iatom,3)
                  dn7b(l2) = dn7b(l2) + s3s(iatom)

                  t84  = (dn6(1) - dn7(1)) / Z2
                  t85  = (dn6(2) - dn7(2)) / Z2
                  t86  = (dn6(3) - dn7(3)) / Z2
                  t84b = (dn7(1) - dn7b(1)) / Z2
                  t85b = (dn7(2) - dn7b(2)) / Z2
                  t86b = (dn7(3) - dn7b(3)) / Z2
                  t58  = (pj0s - pj1s) / Z2
                  t59  = (pj1s - pj2s) / Z2
                  t60  = (pj2s - pj3s) / Z2

                  d0s = t1 * p0s - t2 * p1s
                  d1s = t1 * p1s - t2 * p2s
                  d2s = t1 * p2s - t2 * p3s
                  d3s = t1 * p3s - t2 * p4s
                  if (l1 .eq. l2) then
                     f1  = sq3
                     d0s = d0s + t50
                     d1s = d1s + t51
                     d2s = d2s + t52
                     d3s = d3s + t53
                     dn3(1) = dn3(1) + t26
                     dn3(2) = dn3(2) + t27
                     dn3(3) = dn3(3) + t28
                     dn4(1) = dn4(1) + t29
                     dn4(2) = dn4(2) + t30
                     dn4(3) = dn4(3) + t31
                     dn4b(1) = dn4b(1) + t32
                     dn4b(2) = dn4b(2) + t33
                     dn4b(3) = dn4b(3) + t34
                  endif

                  t61 = (d0s - d1s) / Z2
                  t62 = (d1s - d2s) / Z2
                  t63 = (d2s - d3s) / Z2
                  t96 = (dn3(1) - dn4(1)) / Z2
                  t97 = (dn3(2) - dn4(2)) / Z2
                  t98 = (dn3(3) - dn4(3)) / Z2

                  lij = 3
                  if (ifunct .eq. jfunct) then
                     lij = l1
                  endif
                  do l3 = 1, lij
                     t1 = Q(l3) - r(Nuc(jfunct),l3)
                     t2 = Q(l3) - r(iatom,l3)

                     s0p = t1 * s0s(iatom) - t2 * s1s(iatom)
                     s1p = t1 * s1s(iatom) - t2 * s2s(iatom)
                     t70 = (s0p - s1p) / Z2
                     t71 = (s1p - (t1 * s2s(iatom) - t2 * s3s(iatom))) / Z2

                     d0p = t1 * d0s - t2 * d1s
                     d1p = t1 * d1s - t2 * d2s
                     d2p = t1 * d2s - t2 * d3s

                     pi0p = t1 * p0s  - t2 * p1s
                     pi1p = t1 * p1s  - t2 * p2s
                     pi2p = t1 * p2s  - t2 * p3s
                     pj0p = t1 * pj0s - t2 * pj1s
                     pj1p = t1 * pj1s - t2 * pj2s
                     pj2p = t1 * pj2s - t2 * pj3s

                     ! dn8 and dn8b (pi||pk),   dn9 and dn9b (pj||pk)
                     dn8(1)  = t1 * dn(1) - t2 * dn1(1)
                     dn8(2)  = t1 * dn(2) - t2 * dn1(2)
                     dn8(3)  = t1 * dn(3) - t2 * dn1(3)
                     dn8(l3) = dn8(l3)+p1s
                     dn8b(1)  = t1 * dn1(1) - t2 * dn2(1)
                     dn8b(2)  = t1 * dn1(2) - t2 * dn2(2)
                     dn8b(3)  = t1 * dn1(3) - t2 * dn2(3)
                     dn8b(l3) = dn8b(l3)+p2s

                     dn9(1)  = t1 * dn6(1) - t2 * dn7(1)
                     dn9(2)  = t1 * dn6(2) - t2 * dn7(2)
                     dn9(3)  = t1 * dn6(3) - t2 * dn7(3)
                     dn9(l3) = dn9(l3) + pj1s
                     dn9b(1)  = t1 * dn7(1) - t2 * dn7b(1)
                     dn9b(2)  = t1 * dn7(2) - t2 * dn7b(2)
                     dn9b(3)  = t1 * dn7(3) - t2 * dn7b(3)
                     dn9b(l3) = dn9b(l3) + pj2s

                     ! dn5 (dij || pk) dn5b (dij ||pk) order 1
                     dn5(1) = t1 * dn3(1) - t2 * dn4(1)
                     dn5(2) = t1 * dn3(2) - t2 * dn4(2)
                     dn5(3) = t1 * dn3(3) - t2 * dn4(3)
                     dn5(l3) = dn5(l3) + d1s
                     dn5b(1)  = t1 * dn4(1) - t2 * dn4b(1)
                     dn5b(2)  = t1 * dn4(2) - t2 * dn4b(2)
                     dn5b(3)  = t1 * dn4(3) - t2 * dn4b(3)
                     dn5b(l3) = dn5b(l3) + d2s

                     if (l1.eq.l3) then
                        d0p  = d0p + t58
                        d1p  = d1p + t59
                        d2p  = d2p + t60
                        pi0p = pi0p + t50
                        pi1p = pi1p + t51
                        pi2p = pi2p + t52
                        dn5(1)  = dn5(1) + t84
                        dn5(2)  = dn5(2) + t85
                        dn5(3)  = dn5(3) + t86
                        dn5b(1) = dn5b(1) + t84b
                        dn5b(2) = dn5b(2) + t85b
                        dn5b(3) = dn5b(3) + t86b
                        dn8(1)  = dn8(1) + t26
                        dn8(2)  = dn8(2) + t27
                        dn8(3)  = dn8(3) + t28
                        dn8b(1) = dn8b(1) + t29
                        dn8b(2) = dn8b(2) + t30
                        dn8b(3) = dn8b(3) + t31
                     endif
                     if (l2.eq.l3) then
                        d0p  = d0p + t54
                        d1p  = d1p + t55
                        d2p  = d2p + t56
                        pj0p = pj0p + t50
                        pj1p = pj1p + t51
                        pj2p = pj2p + t52
                        dn5(1)  = dn5(1) + t81
                        dn5(2)  = dn5(2) + t82
                        dn5(3)  = dn5(3) + t83
                        dn5b(1) = dn5b(1) + t81b
                        dn5b(2) = dn5b(2) + t82b
                        dn5b(3) = dn5b(3) + t83b
                        dn9(1)  = dn9(1) + t26
                        dn9(2)  = dn9(2) + t27
                        dn9(3)  = dn9(3) + t28
                        dn9b(1) = dn9b(1) + t29
                        dn9b(2) = dn9b(2) + t30
                        dn9b(3) = dn9b(3) + t31
                     endif

                     t64 = (d0p - d1p) / Z2
                     t65 = (d1p - d2p) / Z2
                     t66 = (pi0p - pi1p) / Z2
                     t67 = (pi1p - pi2p) / Z2
                     t68 = (pj0p - pj1p) / Z2
                     t69 = (pj1p - pj2p) / Z2
                     t90 = (dn9(1) - dn9b(1)) / Z2
                     t91 = (dn9(2) - dn9b(2)) / Z2
                     t92 = (dn9(3) - dn9b(3)) / Z2
                     t93 = (dn8(1) - dn8b(1)) / Z2
                     t94 = (dn8(2) - dn8b(2)) / Z2
                     t95 = (dn8(3) - dn8b(3)) / Z2

                     lk = l3
                     if (ifunct .eq. jfunct) then
                        lk = min(l3, Ll(l1) - Ll(l3) + l2)
                     endif
                     do l4 = 1, lk
                        f2 = 1.D0
                        t1 = Q(l4) - r(Nuc(jfunct),l4)
                        t2 = Q(l4) - r(iatom,l4)
                        tna  = t1 * d0p - t2 * d1p
                        tn1a = t1 * d1p - t2 * d2p

                        ! dn10 : (dij || dkl) nuclear derivative needed
                        dn10(1)  = t1 * dn5(1) - t2 * dn5b(1)
                        dn10(2)  = t1 * dn5(2) - t2 * dn5b(2)
                        dn10(3)  = t1 * dn5(3) - t2 * dn5b(3)
                        dn10(l4) = dn10(l4)+d1p

                        d0pl = t1 * d0s - t2 * d1s
                        d1pl = t1 * d1s - t2 * d2s
                        pj0d = t1 * pj0p - t2 * pj1p
                        pj1d = t1 * pj1p - t2 * pj2p
                        pi0d = t1 * pi0p - t2 * pi1p
                        pi1d = t1 * pi1p - t2 * pi2p

                        if (l4.eq.l1) then
                           tna  = tna  + t68
                           tn1a = tn1a + t69
                           d0pl = d0pl + t58
                           d1pl = d1pl + t59
                           pi0d = pi0d + t70
                           pi1d = pi1d + t71
                           dn10(1) = dn10(1) + t90
                           dn10(2) = dn10(2) + t91
                           dn10(3) = dn10(3) + t92
                        endif
                        if (l4.eq.l2) then
                           tna  = tna  + t66
                           tn1a = tn1a + t67
                           d0pl = d0pl + t54
                           d1pl = d1pl + t55
                           pj0d = pj0d + t70
                           pj1d = pj1d + t71
                           dn10(1) = dn10(1) + t93
                           dn10(2) = dn10(2) + t94
                           dn10(3) = dn10(3) + t95
                        endif
                        if (l4.eq.l3) then
                           f2 = sq3
                           tna  = tna  + t61
                           tn1a = tn1a  + t62
                           pj0d = pj0d + t58
                           pj1d = pj1d + t59
                           pi0d = pi0d + t54
                           pi1d = pi1d + t55
                           dn10(1) = dn10(1) + t96
                           dn10(2) = dn10(2) + t97
                           dn10(3) = dn10(3) + t98
                        endif

                        t72 = (pj0d - pj1d) / Z2
                        t73 = (pi0d - pi1d) / Z2
                        t74 = (d0pl - d1pl) / Z2

                        ! Gradients
                        l12   = Ll(l1) + l2
                        l34   = Ll(l3) + l4
                        k_ind = ifunct + l12-1 + &
                                ((M2-(jfunct+l34-1))*(jfunct+l34-2))/2

                        te = rho(k_ind)* ccoef / (f1 * f2)
                        t4 = te * 2.D0 * a(ifunct,nci)
                        t5 = te * 2.D0 * a(jfunct,ncj)

                        do l5 = 1, 3
                           t1 = Q(l5) - r(Nuc(ifunct),l5)
                           t2 = Q(l5) - r(iatom,l5)

                           fNd = t1 * tna - t2 * tn1a
                           if (l1.eq.l5) then
                              fNd = fNd + t72
                              ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) - te *pj0d
                           endif
                           if (l2.eq.l5) then
                              fNd = fNd + t73
                              ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) - te *pi0d
                           endif
                           if (l3.eq.l5) then
                              fNd = fNd + t74
                              ff(Nuc(jfunct),l5) = ff(Nuc(jfunct),l5) - te *d0pl
                           endif
                           if (l4.eq.l5) then
                              fNd = fNd + t64
                              ff(Nuc(jfunct),l5) = ff(Nuc(jfunct),l5) - te * d0p
                           endif
                           dNf = (r(Nuc(ifunct),l5) - r(Nuc(jfunct),l5)) *tna +&
                                 fNd
                           ff(Nuc(ifunct),l5) = ff(Nuc(ifunct),l5) + t4 * fNd
                           ff(Nuc(jfunct),l5) = ff(Nuc(jfunct),l5) + t5 * dNf
                           ff(iatom,l5)       = ff(iatom,l5) + te * dn10(l5)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         ! End of nuclear attraction
      enddo
      enddo
   enddo
   enddo
   !print*, "d|d", ff

   deallocate(s0s, s1s, s2s, s3s, s4s, s5s, s6s, x0x, x1x, x2x, x3x, x4x)
   return
end subroutine int1G
end module subm_int1G
