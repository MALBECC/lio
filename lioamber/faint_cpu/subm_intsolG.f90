
! GRADIENT VERSION
! 1 e integrals
! using the Obara - saika recursive method.
! debugged ( or supposed to) 28-7-92
! Dario Estrin

module subm_intsolG; contains
subroutine intsolG(frc_qm, frc_mm, natom, ntatom, RMM, d, r, pc, Iz)

   use garcha_mod   , only: M, Md, a, c, Nuc, Ncont, nshell, rmax, NORM
   use liotemp      , only: FUNCT
   use constants_mod, only: pi, pi32

   implicit none
   integer         , intent(in)    :: natom, ntatom, Iz(natom)
   double precision, intent(in)    :: RMM(:), r(ntatom,3), d(natom,natom), &
                                      pc(ntatom)
   double precision, intent(inout) :: frc_qm(natom,3), frc_mm(ntatom,3)

   integer          :: ns, np, nd, iatom, jatom, ifunct, jfunct, nci, ncj, MM, &
                       M2, rho_ind, lk, lij, l1, l2, l3, l4, l5, l12, l34, Ll(3)
   double precision :: SQ3, f1, f2, Q(3), q1, q2, q3, rexp, term0, term, Zij,  &
                       Z2, uf, ccoef, te
   double precision :: s0p, s1p, s2p, sNpi
   double precision :: p0s, p1s, p2s, p3s, p4s, pi0p, pi0d, pi1p, pi1d, pi2p,  &
                       pj0s, pj0p, pj0d, pj1s, pj1p, pj1d, pj2s, pj2p, pj3s,   &
                       pNp, pNd, pN1p, piNs
   double precision :: d0s, d0p, d1s, d1p, d2s, d3s, d2p, d0pl, d1pl, dNs, dNp,&
                       dNd, dNf, dN1s, dN1p
   double precision :: fNs, fNp, fNd
   double precision :: t1, t2, t3, t4, t5, t7, t8, t9, t15, t25, t26, t27, t28,&
                       t29, t30, t31, t32, t33, t34, t50, t51, t52, t53, t54,  &
                       t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65,  &
                       t66, t67, t68, t69, t70, t71, t72, t73, t74, t81, t81b, &
                       t82, t82b, t83, t83b, t84, t84b, t85, t85b, t86, t86b,  &
                       t90, t91, t92, t93, t94, t95, t96, t97, t98
   double precision :: dn(3)  , dn1(3) , dn2(3) , dn3(3) , dn4(3) , dn5(3) , &
                       dn6(3) , dn7(3) , dn8(3) , dn9(3) , dn10(3), dn2b(3), &
                       dn4b(3), dn5b(3), dn7b(3), dn8b(3), dn9b(3), dn11(3), &
                       dn12(3)
   double precision, allocatable :: s0s(:), s1s(:), s2s(:), s3s(:), s4s(:), &
                                    s5s(:), s6s(:), x0x(:,:), x1x(:,:),     &
                                    x2x(:,:), x3x(:,:), x4x(:,:)

   allocate(s0s(ntatom), s1s(ntatom), s2s(ntatom), s3s(ntatom), s4s(ntatom), &
            s5s(ntatom), s6s(ntatom))
   allocate(x0x(ntatom,3), x1x(ntatom,3), x2x(ntatom,3), x3x(ntatom,3), &
            x4x(ntatom,3))

   SQ3 = 1.D0
   if (NORM) SQ3 = sqrt(3.D0)

   ns = nshell(0); np = nshell(1); nd = nshell(2)
   MM  = M * (M +1) / 2
   M2  = 2 * M

   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) / 2
   enddo

   ! Nuclear repulsion
   do iatom = 1, natom
   do jatom = natom+1, ntatom
      t1   = r(iatom,1) - r(jatom,1)
      t2   = r(iatom,2) - r(jatom,2)
      t3   = r(iatom,3) - r(jatom,3)
      term = t1 * t1 + t2 * t2 + t3 * t3
      term = - Iz(iatom) * pc(jatom) / (term * dsqrt(term))

      frc_qm(iatom,1) = frc_qm(iatom,1) + term * t1
      frc_qm(iatom,2) = frc_qm(iatom,2) + term * t2
      frc_qm(iatom,3) = frc_qm(iatom,3) + term * t3

      frc_mm(jatom,1) = frc_mm(jatom,1) - term * t1
      frc_mm(jatom,2) = frc_mm(jatom,2) - term * t2
      frc_mm(jatom,3) = frc_mm(jatom,3) - term * t3
   enddo
   enddo
   !print*, frc_qm
   !print*, frc_mm

   ! (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         Z2   = 2.D0 * Zij
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            term0 = 2.D0 * PI * exp( - rexp) / Zij

            Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
            Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
            Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

            do iatom = natom+1, ntatom
               uf  = ((Q(1) - r(iatom,1)) * (Q(1) - r(iatom,1)) + &
                      (Q(2) - r(iatom,2)) * (Q(2) - r(iatom,2)) + &
                      (Q(3) - r(iatom,3)) * (Q(3) - r(iatom,3))) * Zij

               term = - pc(iatom) * term0
               s0s(iatom) = term * FUNCT(0,uf)
               s1s(iatom) = term * FUNCT(1,uf)

               term = Z2 * s1s(iatom)
               x0x(iatom,1) = term * (Q(1) - r(iatom,1))
               x0x(iatom,2) = term * (Q(2) - r(iatom,2))
               x0x(iatom,3) = term * (Q(3) - r(iatom,3))
            enddo

            rho_ind = ifunct + ((M2 - jfunct) * (jfunct -1)) / 2
            te = RMM(rho_ind) * ccoef
            t4 = 2.D0 * te * a(ifunct,nci)
            t5 = 2.D0 * te * a(jfunct,ncj)

            do l2 = 1, 3
               t1 = Q(l2) - r(Nuc(ifunct),l2)
               t2 = r(Nuc(ifunct),l2) - r(Nuc(jfunct),l2)

               do iatom = natom+1, ntatom
                  piNs = t1 * s0s(iatom) - (Q(l2) - r(iatom,l2)) * s1s(iatom)
                  sNpi = piNs + t2 * s0s(iatom)
                  frc_qm(Nuc(ifunct),l2) = frc_qm(Nuc(ifunct),l2) + t4 * piNs
                  frc_qm(Nuc(jfunct),l2) = frc_qm(Nuc(jfunct),l2) + t5 * sNpi
                  frc_mm(iatom,l2)       = frc_mm(iatom,l2) + te * x0x(iatom,l2)
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo
   !print*, "ss", frc_qm
   !print*, "ss", frc_mm

   ! (p|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1, ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            term0 = 2.D0 * PI * exp( - rexp) / Zij
            Z2    = 2.D0 * Zij
            Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
            Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
            Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

            do iatom = natom+1, ntatom
               q1   = Q(1) - r(iatom,1)
               q2   = Q(2) - r(iatom,2)
               q3   = Q(3) - r(iatom,3)
               uf   = (q1 * q1 + q2 * q2 + q3 * q3) * Zij

               term = - pc(iatom) * term0
               s0s(iatom) = term * FUNCT(0,uf)
               s1s(iatom) = term * FUNCT(1,uf)
               s2s(iatom) = term * FUNCT(2,uf)

               term = Z2 * s1s(iatom)
               x0x(iatom,1) = term * q1
               x0x(iatom,2) = term * q2
               x0x(iatom,3) = term * q3

               term = Z2 * s2s(iatom)
               x1x(iatom,1) = term * q1
               x1x(iatom,2) = term * q2
               x1x(iatom,3) = term * q3
            enddo

            do iatom = natom+1, ntatom
               t50 = (s0s(iatom) - s1s(iatom)) / Z2

               do l1 = 1, 3
                  t1  = Q(l1) - r(Nuc(ifunct),l1)
                  t2  = Q(l1) - r(iatom,l1)
                  p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  p1s = t1 * s1s(iatom) - t2 * s2s(iatom)

                  dn(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
                  dn(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
                  dn(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
                  dn(l1) = dn(l1) + s1s(iatom)

                  rho_ind = ifunct + l1-1 + ((M2 - jfunct) * (jfunct -1)) / 2
                  te = RMM(rho_ind) * ccoef
                  t4 = 2.0D0 * te * a(ifunct, nci)
                  t5 = 2.0D0 * te * a(jfunct, ncj)

                  do l2 = 1, 3
                     t1  = Q(l2) - r(Nuc(ifunct),l2)
                     t2  = Q(l2) - r(iatom,l2)
                     dNs = t1 * p0s - t2 * p1s

                     if (l1 .eq. l2) then
                        dNs = dNs + t50
                        frc_qm(Nuc(ifunct),l2) = frc_qm(Nuc(ifunct),l2) - &
                                                 te * s0s(iatom)
                     endif
                     pNp = dNs + (r(Nuc(ifunct),l2) - r(Nuc(jfunct),l2)) * p0s

                     frc_qm(Nuc(ifunct),l2) = frc_qm(Nuc(ifunct),l2) + t4 * dNs
                     frc_qm(Nuc(jfunct),l2) = frc_qm(Nuc(jfunct),l2) + t5 * pNp
                     frc_mm(iatom,l2)       = frc_mm(iatom,l2)       + te*dn(l2)
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo
   !print*, "ps", frc_qm
   !print*, "ps", frc_mm

   ! (p|p)
   do ifunct = ns+1, ns+np , 3
   do jfunct = ns+1, ifunct, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            Z2    = 2.D0 * Zij
            term0 = 2.D0 * PI * exp( - rexp) / Zij
            Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
            Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
            Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

            do iatom = natom+1, ntatom
               q1   = Q(1) - r(iatom,1)
               q2   = Q(2) - r(iatom,2)
               q3   = Q(3) - r(iatom,3)
               uf   = (q1 * q1 + q2 * q2 + q3 * q3) * Zij

               term = - pc(iatom) * term0
               s0s(iatom) = term * FUNCT(0,uf)
               s1s(iatom) = term * FUNCT(1,uf)
               s2s(iatom) = term * FUNCT(2,uf)
               s3s(iatom) = term * FUNCT(3,uf)

               term = Z2 * s1s(iatom)
               x0x(iatom,1) = term * q1
               x0x(iatom,2) = term * q2
               x0x(iatom,3) = term * q3

               term = Z2 * s2s(iatom)
               x1x(iatom,1) = term * q1
               x1x(iatom,2) = term * q2
               x1x(iatom,3) = term * q3

               term = Z2 * s3s(iatom)
               x2x(iatom,1) = term * q1
               x2x(iatom,2) = term * q2
               x2x(iatom,3) = term * q3
            enddo

            do iatom = natom+1, ntatom
               t15 = (s0s(iatom)   - s1s(iatom))   / Z2
               t25 = (s1s(iatom)   - s2s(iatom))   / Z2
               t26 = (x0x(iatom,1) - x1x(iatom,1)) / Z2
               t27 = (x0x(iatom,2) - x1x(iatom,2)) / Z2
               t28 = (x0x(iatom,3) - x1x(iatom,3)) / Z2

               do l1 = 1, 3
                  t1 = Q(l1) - r(Nuc(ifunct),l1)
                  t2 = Q(l1) - r(iatom,l1)
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

                  lij = 3
                  if (ifunct .eq. jfunct) lij = l1
                  do l2 = 1, lij
                     t1  = Q(l2) - r(Nuc(ifunct),l2)
                     t2  = Q(l2) - r(iatom,l2)
                     t3  = Q(l2) - r(Nuc(jfunct),l2)

                     s0p  = t3 * s0s(iatom) - t2 * s1s(iatom)
                     s1p  = t3 * s1s(iatom) - t2 * s2s(iatom)
                     pNp  = t3 * p0s - t2 * p1s
                     pN1p = t3 * p1s - t2 * p2s
                     t29  = (s0p - s1p) / Z2

                     dn2(1) = t3 * dn(1) - t2 * dn1(1)
                     dn2(2) = t3 * dn(2) - t2 * dn1(2)
                     dn2(3) = t3 * dn(3) - t2 * dn1(3)
                     dn2(l2) = dn2(l2) + p1s

                     d0s = t1 * p0s - t2 * p1s
                     d1s = t1 * p1s - t2 * p2s

                     if (l1 .eq. l2) then
                        pNp    = pNp  + t15
                        pN1p   = pN1p + t25
                        dn2(1) = dn2(1) + t26
                        dn2(2) = dn2(2) + t27
                        dn2(3) = dn2(3) + t28
                     endif

                     rho_ind = ifunct + l1-1 + ((M2 - (jfunct + l2-1)) * &
                               (jfunct + l2-2)) / 2
                     te = RMM(rho_ind) * ccoef
                     t4 = 2.0D0 * te * a(ifunct, nci)
                     t5 = 2.0D0 * te * a(jfunct, ncj)

                     do l3 = 1, 3
                        dNp = (Q(l3) - r(Nuc(ifunct),l3)) * pNp - &
                              (Q(l3) - r(iatom,l3))       * pN1p

                        if (l1 .eq. l3) then
                           dNp = dNp + t29
                           frc_qm(Nuc(ifunct),l3) = frc_qm(Nuc(ifunct),l3) - &
                                                    te * s0p
                        endif
                        if (l2 .eq. l3) then
                           dNp = dNp + t30
                           frc_qm(Nuc(jfunct),l3) = frc_qm(Nuc(jfunct),l3) - &
                                                    te * p0s
                        endif
                        pNd = dNp + (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) *pNp

                        frc_qm(Nuc(ifunct),l3) = frc_qm(Nuc(ifunct),l3) + t4*dNp
                        frc_qm(Nuc(jfunct),l3) = frc_qm(Nuc(jfunct),l3) + t5*pNd
                        frc_mm(iatom,l3) = frc_mm(iatom,l3) + te * dn2(l3)
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo
   !print*, "pp", frc_qm
   !print*, "pp", frc_mm

   ! (d|s)
   do ifunct = ns+np+1, M, 6
   do jfunct = 1      , ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            Z2    = 2.D0 * Zij
            term0 = 2.D0 * PI * exp( - rexp) / Zij
            Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
            Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
            Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

            do iatom = natom+1, ntatom
               q1   = Q(1) - r(iatom,1)
               q2   = Q(2) - r(iatom,2)
               q3   = Q(3) - r(iatom,3)
               uf   = (q1 * q1 + q2 * q2 + q3 * q3) * Zij

               term = - pc(iatom) * term0
               s0s(iatom) = term * FUNCT(0,uf)
               s1s(iatom) = term * FUNCT(1,uf)
               s2s(iatom) = term * FUNCT(2,uf)
               s3s(iatom) = term * FUNCT(3,uf)

               term = Z2 * s1s(iatom)
               x0x(iatom,1) = term * q1
               x0x(iatom,2) = term * q2
               x0x(iatom,3) = term * q3

               term = Z2 * s2s(iatom)
               x1x(iatom,1) = term * q1
               x1x(iatom,2) = term * q2
               x1x(iatom,3) = term * q3

               term = Z2 * s3s(iatom)
               x2x(iatom,1) = term * q1
               x2x(iatom,2) = term * q2
               x2x(iatom,3) = term * q3
            enddo

            do iatom = natom+1, ntatom
               t7  = (s0s(iatom)   - s1s(iatom)  ) / Z2
               t8  = (s1s(iatom)   - s2s(iatom)  ) / Z2
               t26 = (x0x(iatom,1) - x1x(iatom,1)) / Z2
               t27 = (x0x(iatom,2) - x1x(iatom,2)) / Z2
               t28 = (x0x(iatom,3) - x1x(iatom,3)) / Z2

               do l1 = 1, 3
                  t1 = Q(l1) - r(Nuc(ifunct),l1)
                  t2 = Q(l1) - r(iatom,l1)
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
                     dNs  = t1 * p0s - t2 * p1s
                     dN1s = t1 * p1s - t2 * p2s
                     pj0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                     pj1s = t1 * s1s(iatom) - t2 * s2s(iatom)
                     t29  = (pj0s - pj1s) / Z2

                     dn2(1)  = t1 * dn(1) - t2 * dn1(1)
                     dn2(2)  = t1 * dn(2) - t2 * dn1(2)
                     dn2(3)  = t1 * dn(3) - t2 * dn1(3)
                     dn2(l2) = dn2(l2) + p1s

                     f1 = 1.0D0
                     if (l1 .eq. l2) then
                        dNs    = dNs  + t7
                        dN1s   = dN1s + t8
                        f1     = SQ3
                        dn2(1) = dn2(1) + t26
                        dn2(2) = dn2(2) + t27
                        dn2(3) = dn2(3) + t28
                     endif

                     l12 = Ll(l1) + l2
                     rho_ind = ifunct + l12 -1 + ((M2 - jfunct) * (jfunct -1))/2

                     te = RMM(rho_ind)* ccoef / f1
                     t4 = 2.0D0 * te * a(ifunct,nci)
                     t5 = 2.0D0 * te * a(jfunct,ncj)

                     do l3 = 1, 3
                        dNp = (Q(l3) - r(Nuc(jfunct),l3)) * dNs - &
                              (Q(l3) - r(iatom,l3)      ) * dN1s
                        if (l1 .eq. l3) then
                           dNp = dNp + t29
                           frc_qm(Nuc(ifunct),l3) = frc_qm(Nuc(ifunct),l3) - &
                                                    te * pj0s
                        endif
                        if (l2 .eq. l3) then
                           dNp = dNp + t30
                           frc_qm(Nuc(ifunct),l3) = frc_qm(Nuc(ifunct),l3) - &
                                                    te * p0s
                        endif
                        fNs = dNp - (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) *dNs

                        frc_qm(Nuc(ifunct),l3) = frc_qm(Nuc(ifunct),l3) + t4*fNs
                        frc_qm(Nuc(jfunct),l3) = frc_qm(Nuc(jfunct),l3) + t5*dNp
                        frc_mm(iatom,l3)       = frc_mm(iatom,l3) + te * dn2(l3)
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo
   !print*, "ds", frc_qm
   !print*, "ds", frc_mm

   ! (d|p)
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+1   , ns+np , 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            Z2    = 2.D0 * Zij
            term0 = 2.D0 * PI * exp( - rexp) / Zij
            Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
            Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
            Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

            do iatom = natom+1, ntatom
               q1   = Q(1) - r(iatom,1)
               q2   = Q(2) - r(iatom,2)
               q3   = Q(3) - r(iatom,3)
               uf   = (q1 * q1 + q2 * q2 + q3 * q3) * Zij

               term = - pc(iatom) * term0
               s0s(iatom) = term * FUNCT(0,uf)
               s1s(iatom) = term * FUNCT(1,uf)
               s2s(iatom) = term * FUNCT(2,uf)
               s3s(iatom) = term * FUNCT(3,uf)
               s4s(iatom) = term * FUNCT(4,uf)

               term = Z2 * s1s(iatom)
               x0x(iatom,1) = term * q1
               x0x(iatom,2) = term * q2
               x0x(iatom,3) = term * q3

               term = Z2 * s2s(iatom)
               x1x(iatom,1) = term * q1
               x1x(iatom,2) = term * q2
               x1x(iatom,3) = term * q3

               term = Z2 * s3s(iatom)
               x2x(iatom,1) = term * q1
               x2x(iatom,2) = term * q2
               x2x(iatom,3) = term * q3

               term = Z2 * s4s(iatom)
               x3x(iatom,1) = term * q1
               x3x(iatom,2) = term * q2
               x3x(iatom,3) = term * q3
            enddo

            do iatom = natom+1, ntatom
               t7  = (s0s(iatom)   - s1s(iatom)  ) / Z2
               t8  = (s1s(iatom)   - s2s(iatom)  ) / Z2
               t9  = (s2s(iatom)   - s3s(iatom)  ) / Z2
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

                  dn2(1)  = t1 * x2x(iatom,1) - t2 * x3x(iatom,1)
                  dn2(2)  = t1 * x2x(iatom,2) - t2 * x3x(iatom,2)
                  dn2(3)  = t1 * x2x(iatom,3) - t2 * x3x(iatom,3)
                  dn2(l1) = dn2(l1) + s3s(iatom)

                  t51 = (dn(1) - dn1(1)) / Z2
                  t52 = (dn(2) - dn1(2)) / Z2
                  t53 = (dn(3) - dn1(3)) / Z2

                  do l2 = 1, l1
                     t1 = Q(l2) - r(Nuc(ifunct),l2)
                     t2 = Q(l2) - r(iatom,l2)
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

                     ! dn6 (pj || s) order 0
                     dn6(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
                     dn6(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
                     dn6(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
                     dn6(l2) = dn6(l2) + s1s(iatom)

                     ! dn7 (pj || s) order 1
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
                        f1  = SQ3
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
                        t1 = Q(l3) - r(Nuc(jfunct),l3)
                        t2 = Q(l3) - r(iatom,l3)
                        dNp  = t1 * d0s - t2 * d1s
                        dN1p = t1 * d1s - t2 * d2s

                        pi0p = t1 * p0s - t2 * p1s
                        pi1p = t1 * p1s - t2 * p2s
                        pj0p = t1 * pj0s - t2 * pj1s
                        pj1p = t1 * pj1s - t2 * pj2s

                        ! dn5 (dij || Pk ) order 0
                        dn5(1)  = t1 * dn3(1) - t2 * dn4(1)
                        dn5(2)  = t1 * dn3(2) - t2 * dn4(2)
                        dn5(3)  = t1 * dn3(3) - t2 * dn4(3)
                        dn5(l3) = dn5(l3) + d1s

                        if (l1 .eq. l3) then
                           dNp    = dNp  + (pj0s - pj1s) / Z2
                           dN1p   = dN1p + (pj1s - pj2s) / Z2
                           pi0p   = pi0p + t7
                           pi1p   = pi1p + t8
                           dn5(1) = dn5(1) + t54
                           dn5(2) = dn5(2) + t55
                           dn5(3) = dn5(3) + t56
                        endif
                        if (l2 .eq. l3) then
                           dNp    = dNp  + (p0s - p1s) / Z2
                           dN1p   = dN1p + (p1s - p2s) / Z2
                           pj0p   = pj0p + t7
                           pj1p   = pj1p + t8
                           dn5(1) = dn5(1) + t51
                           dn5(2) = dn5(2) + t52
                           dn5(3) = dn5(3) + t53
                        endif

                        l12 = Ll(l1) + l2
                        rho_ind = ifunct + l12 -1 + (M2 - jfunct - l3 +1) * &
                                                    (jfunct + l3 -2) / 2
                        te = RMM(rho_ind) * ccoef / f1
                        t4 = 2.0D0 * te * a(ifunct,nci)
                        t5 = 2.0D0 * te * a(jfunct,ncj)

                        do l4 = 1, 3
                           dNd = (Q(l4) - r(Nuc(jfunct),l4)) * dNp - &
                                 (Q(l4) - r(iatom,l4)      ) * dN1p

                           if (l1 .eq. l4) then
                              frc_qm(Nuc(ifunct),l4) = frc_qm(Nuc(ifunct),l4) &
                                                     - te * pj0p
                              dNd = dNd + (pj0p - pj1p) / Z2
                           endif
                           if (l2 .eq. l4) then
                              frc_qm(Nuc(ifunct),l4) = frc_qm(Nuc(ifunct),l4) &
                                                     - te * pi0p
                              dNd = dNd + (pi0p - pi1p) / Z2
                           endif
                           if (l3.eq.l4) then
                              frc_qm(Nuc(jfunct),l4) = frc_qm(Nuc(jfunct),l4) &
                                                     - te * d0s
                              dNd = dNd + (d0s - d1s) / Z2
                           endif
                           fNp = dNd - dNp * &
                                       (r(Nuc(ifunct),l4) - r(Nuc(jfunct),l4))

                           frc_qm(Nuc(ifunct),l4) = frc_qm(Nuc(ifunct),l4) + &
                                                    t4 * fNp
                           frc_qm(Nuc(jfunct),l4) = frc_qm(Nuc(jfunct),l4) + &
                                                    t5 * dNd
                           frc_mm(iatom,l4)       = frc_mm(iatom,l4)       + &
                                                    te * dn5(l4)
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
   !print*, "dp", frc_qm
   !print*, "dp", frc_mm

   ! (d|d)
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+np+1, ifunct, 6
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            Z2    = 2.D0 * Zij
            term0 = 2.D0 * PI * exp( - rexp) / Zij
            Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
            Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
            Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                    a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

            do iatom = natom+1, ntatom
               q1   = Q(1) - r(iatom,1)
               q2   = Q(2) - r(iatom,2)
               q3   = Q(3) - r(iatom,3)
               uf   = (q1 * q1 + q2 * q2 + q3 * q3) * Zij

               term = - pc(iatom) * term0
               s0s(iatom) = term * FUNCT(0,uf)
               s1s(iatom) = term * FUNCT(1,uf)
               s2s(iatom) = term * FUNCT(2,uf)
               s3s(iatom) = term * FUNCT(3,uf)
               s4s(iatom) = term * FUNCT(4,uf)
               s5s(iatom) = term * FUNCT(5,uf)
               s6s(iatom) = term * FUNCT(6,uf)

               term = Z2 * s1s(iatom)
               x0x(iatom,1) = term * q1
               x0x(iatom,2) = term * q2
               x0x(iatom,3) = term * q3

               term = Z2 * s2s(iatom)
               x1x(iatom,1) = term * q1
               x1x(iatom,2) = term * q2
               x1x(iatom,3) = term * q3

               term = Z2 * s3s(iatom)
               x2x(iatom,1) = term * q1
               x2x(iatom,2) = term * q2
               x2x(iatom,3) = term * q3

               term = Z2 * s4s(iatom)
               x3x(iatom,1) = term * q1
               x3x(iatom,2) = term * q2
               x3x(iatom,3) = term * q3

               term = Z2 * s5s(iatom)
               x4x(iatom,1) = term * q1
               x4x(iatom,2) = term * q2
               x4x(iatom,3) = term * q3
            enddo

            do iatom = natom +1, ntatom
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
                  t1 = Q(l1) - r(Nuc(ifunct),l1)
                  t2 = Q(l1) - r(iatom,l1)
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

                  t81  = (dn(1)  - dn1(1)) / Z2
                  t82  = (dn(2)  - dn1(2)) / Z2
                  t83  = (dn(3)  - dn1(3)) / Z2
                  t81b = (dn1(1) - dn2(1)) / Z2
                  t82b = (dn1(2) - dn2(2)) / Z2
                  t83b = (dn1(3) - dn2(3)) / Z2

                  do l2 = 1, l1
                     f1 = 1.D0
                     t1 = Q(l2) - r(Nuc(ifunct),l2)
                     t2 = Q(l2) - r(iatom,l2)
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

                     ! dn6 (pj || s) order 0
                     dn6(1)  = t1 * x0x(iatom,1) - t2 * x1x(iatom,1)
                     dn6(2)  = t1 * x0x(iatom,2) - t2 * x1x(iatom,2)
                     dn6(3)  = t1 * x0x(iatom,3) - t2 * x1x(iatom,3)
                     dn6(l2) = dn6(l2) + s1s(iatom)

                     ! dn7 (pj || s) order 1
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
                     t58  = (pj0s   - pj1s)    / Z2
                     t59  = (pj1s   - pj2s)    / Z2
                     t60  = (pj2s   - pj3s)    / Z2

                     d0s = t1 * p0s - t2 * p1s
                     d1s = t1 * p1s - t2 * p2s
                     d2s = t1 * p2s - t2 * p3s
                     d3s = t1 * p3s - t2 * p4s

                     if (l1 .eq. l2) then
                        f1  = SQ3
                        d0s     = d0s + t50
                        d1s     = d1s + t51
                        d2s     = d2s + t52
                        d3s     = d3s + t53
                        dn3(1)  = dn3(1) + t26
                        dn3(2)  = dn3(2) + t27
                        dn3(3)  = dn3(3) + t28
                        dn4(1)  = dn4(1) + t29
                        dn4(2)  = dn4(2) + t30
                        dn4(3)  = dn4(3) + t31
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
                     if (ifunct .eq. jfunct) lij = l1

                     do l3 = 1, lij
                        t1 = Q(l3) - r(Nuc(jfunct),l3)
                        t2 = Q(l3) - r(iatom,l3)
                        s0p = t1 * s0s(iatom) - t2 * s1s(iatom)
                        s1p = t1 * s1s(iatom) - t2 * s2s(iatom)
                        s2p = t1 * s2s(iatom) - t2 * s3s(iatom)
                        t70 = (s0p - s1p) / Z2
                        t71 = (s1p - s2p) / Z2

                        d0p  = t1 * d0s  - t2 * d1s
                        d1p  = t1 * d1s  - t2 * d2s
                        d2p  = t1 * d2s  - t2 * d3s
                        pi0p = t1 * p0s  - t2 * p1s
                        pi1p = t1 * p1s  - t2 * p2s
                        pi2p = t1 * p2s  - t2 * p3s
                        pj0p = t1 * pj0s - t2 * pj1s
                        pj1p = t1 * pj1s - t2 * pj2s
                        pj2p = t1 * pj2s - t2 * pj3s

                        ! dn8 and dn8b (pi||pk)
                        dn8(1)   = t1 * dn(1) - t2 * dn1(1)
                        dn8(2)   = t1 * dn(2) - t2 * dn1(2)
                        dn8(3)   = t1 * dn(3) - t2 * dn1(3)
                        dn8(l3)  = dn8(l3) + p1s
                        dn8b(1)  = t1 * dn1(1) - t2 * dn2(1)
                        dn8b(2)  = t1 * dn1(2) - t2 * dn2(2)
                        dn8b(3)  = t1 * dn1(3) - t2 * dn2(3)
                        dn8b(l3) = dn8b(l3) + p2s

                        ! dn9 and dn9b (pj||pk)
                        dn9(1)  = t1 * dn6(1) - t2 * dn7(1)
                        dn9(2)  = t1 * dn6(2) - t2 * dn7(2)
                        dn9(3)  = t1 * dn6(3) - t2 * dn7(3)
                        dn9(l3) = dn9(l3) + pj1s

                        dn9b(1)  = t1 * dn7(1) - t2 * dn7b(1)
                        dn9b(2)  = t1 * dn7(2) - t2 * dn7b(2)
                        dn9b(3)  = t1 * dn7(3) - t2 * dn7b(3)
                        dn9b(l3) = dn9b(l3) + pj2s

                        ! dn5 (dij || pk) dn5b (dij ||pk) order 1
                        dn5(1)  = t1 * dn3(1) - t2 * dn4(1)
                        dn5(2)  = t1 * dn3(2) - t2 * dn4(2)
                        dn5(3)  = t1 * dn3(3) - t2 * dn4(3)
                        dn5(l3) = dn5(l3) + d1s

                        dn5b(1)  = t1 * dn4(1) - t2 * dn4b(1)
                        dn5b(2)  = t1 * dn4(2) - t2 * dn4b(2)
                        dn5b(3)  = t1 * dn4(3) - t2 * dn4b(3)
                        dn5b(l3) = dn5b(l3) + d2s

                        if (l1 .eq. l3) then
                           d0p     = d0p + t58
                           d1p     = d1p + t59
                           d2p     = d2p + t60
                           pi0p    = pi0p + t50
                           pi1p    = pi1p + t51
                           pi2p    = pi2p + t52
                           dn5(1)  = dn5(1)  + t84
                           dn5(2)  = dn5(2)  + t85
                           dn5(3)  = dn5(3)  + t86
                           dn5b(1) = dn5b(1) + t84b
                           dn5b(2) = dn5b(2) + t85b
                           dn5b(3) = dn5b(3) + t86b
                           dn8(1)  = dn8(1)  + t26
                           dn8(2)  = dn8(2)  + t27
                           dn8(3)  = dn8(3)  + t28
                           dn8b(1) = dn8b(1) + t29
                           dn8b(2) = dn8b(2) + t30
                           dn8b(3) = dn8b(3) + t31
                        endif
                        if (l2 .eq. l3) then
                           d0p     = d0p     + t54
                           d1p     = d1p     + t55
                           d2p     = d2p     + t56
                           pj0p    = pj0p    + t50
                           pj1p    = pj1p    + t51
                           pj2p    = pj2p    + t52
                           dn5(1)  = dn5(1)  + t81
                           dn5(2)  = dn5(2)  + t82
                           dn5(3)  = dn5(3)  + t83
                           dn5b(1) = dn5b(1) + t81b
                           dn5b(2) = dn5b(2) + t82b
                           dn5b(3) = dn5b(3) + t83b
                           dn9(1)  = dn9(1)  + t26
                           dn9(2)  = dn9(2)  + t27
                           dn9(3)  = dn9(3)  + t28
                           dn9b(1) = dn9b(1) + t29
                           dn9b(2) = dn9b(2) + t30
                           dn9b(3) = dn9b(3) + t31
                        endif

                        t64 = (d0p  - d1p ) / Z2
                        t65 = (d1p  - d2p ) / Z2
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
                        if (ifunct .eq. jfunct) lk = min(l3, Ll(l1) - Ll(l3)+l2)
                        do l4 = 1, lk
                           f2 = 1.D0
                           t1 = Q(l4) - r(Nuc(jfunct),l4)
                           t2 = Q(l4) - r(iatom,l4)

                           ! dn10 (dij || dkl) nuclear derivative
                           dn10(1)  = t1 * dn5(1) - t2 * dn5b(1)
                           dn10(2)  = t1 * dn5(2) - t2 * dn5b(2)
                           dn10(3)  = t1 * dn5(3) - t2 * dn5b(3)
                           dn10(l4) = dn10(l4)+d1p

                           dNp  = t1 * d0p - t2 * d1p
                           dN1p = t1 * d1p - t2 * d2p
                           d0pl = t1 * d0s - t2 * d1s
                           d1pl = t1 * d1s - t2 * d2s
                           pj0d = t1 * pj0p - t2 * pj1p
                           pj1d = t1 * pj1p - t2 * pj2p
                           pi0d = t1 * pi0p - t2 * pi1p
                           pi1d = t1 * pi1p - t2 * pi2p

                           if (l4 .eq. l1) then
                              dNp  = dNp + t68
                              dN1p = dN1p + t69
                              d0pl = d0pl + t58
                              d1pl = d1pl + t59
                              pi0d = pi0d + t70
                              pi1d = pi1d + t71
                              dn10(1) = dn10(1) + t90
                              dn10(2) = dn10(2) + t91
                              dn10(3) = dn10(3) + t92
                           endif
                           if (l4 .eq. l2) then
                              dNp  = dNp + t66
                              dN1p = dN1p + t67
                              d0pl = d0pl + t54
                              d1pl = d1pl + t55
                              pj0d = pj0d + t70
                              pj1d = pj1d + t71
                              dn10(1) = dn10(1) + t93
                              dn10(2) = dn10(2) + t94
                              dn10(3) = dn10(3) + t95
                           endif
                           if (l4 .eq. l3) then
                              f2 = SQ3
                              dNp  = dNp  + t61
                              dN1p = dN1p + t62
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

                           l12 = Ll(l1) + l2
                           l34 = Ll(l3) + l4
                           rho_ind = (M2 - jfunct - l34 +1) * &
                                      (jfunct + l34 -2) / 2 + ifunct + l12 -1

                           te = RMM(rho_ind) * ccoef / (f1 * f2)
                           t4 = 2.0D0 * te * a(ifunct,nci)
                           t5 = 2.0D0 * te * a(jfunct,ncj)

                           do l5 = 1, 3
                              fNd = (Q(l5) - r(Nuc(ifunct),l5)) * dNp - &
                                    (Q(l5) - r(iatom,l5)      ) * dN1p

                              if (l1 .eq. l5) then
                                 fNd = fNd + t72
                                 frc_qm(Nuc(ifunct),l5) =frc_qm(Nuc(ifunct),l5)&
                                                        - te * pj0d
                              endif
                              if (l2 .eq. l5) then
                                 fNd = fNd + t73
                                 frc_qm(Nuc(ifunct),l5) =frc_qm(Nuc(ifunct),l5)&
                                                        - te * pi0d
                              endif
                              if (l3 .eq. l5) then
                                 fNd = fNd + t74
                                 frc_qm(Nuc(jfunct),l5) =frc_qm(Nuc(jfunct),l5)&
                                                        - te * d0pl
                              endif
                              if (l4 .eq. l5) then
                                 fNd = fNd + t64
                                 frc_qm(Nuc(jfunct),l5) =frc_qm(Nuc(jfunct),l5)&
                                                        - te * d0p
                              endif
                              dNf = fNd + dNp * &
                                   (r(Nuc(ifunct),l5) - r(Nuc(jfunct),l5))

                              frc_qm(Nuc(ifunct),l5) = frc_qm(Nuc(ifunct),l5) +&
                                                       t4 * fNd
                              frc_qm(Nuc(jfunct),l5) = frc_qm(Nuc(jfunct),l5) +&
                                                       t5 * dNf
                              frc_mm(iatom,l5)       = frc_mm(iatom,l5)       +&
                                                       te * dn10(l5)
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
   !print*, "dd", frc_qm
   !print*, "dd", frc_mm

   deallocate(s0s, s1s, s2s, s3s, s4s, s5s, s6s, x0x, x1x, x2x, x3x, x4x)
   return

end subroutine
end module subm_intsolG
