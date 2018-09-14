!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! 1e integrals for the solvent point charges interaction with electronic       !
! density. Calculates the Fock matrix elements vía the Obara-Saika recursive   !
! methods. This is the same as nuclear attraction elements.                    !
!                                                                              !
! Written by Dario Estrin, Buenos Aires, August 1994.                          !
! Refactored by Federico Pedron, Buenos Aires, August 2018.                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_intsol
contains
subroutine intsol(Rho, Hmat, Iz, pc, r, d, natom, ntatom, E1s, Ens, elec)

   use liotemp      , only: FUNCT
   use garcha_mod   , only: a, c, nuc, ncont, rmax, nshell, M, NORM, Md
   use constants_mod, only: pi

   implicit none
   integer         , intent(in)    :: natom, ntatom, Iz(natom)
   logical         , intent(in)    :: elec
   double precision, intent(in)    :: pc(ntatom), r(ntatom,3), d(natom,natom)
   double precision, intent(out)   :: E1s, Ens
   double precision, intent(inout) :: Rho(:), Hmat(:)

   integer           :: M11, ns, np, nd, iatom, jatom, ifunct, jfunct, nci, &
                        ncj, Ll(3), lk, lij, l1, l2, l3, l4, M2, vecmat_ind
   double precision  :: sq3, rexp, ccoef, term, tna, uf, Z2, Zij, t1, t2, f1, &
                        f2, p3s, p2s, p1s, p0s, pj2s, pj1s, pj1p, pj0s, pj0p, &
                        pi1p, pi0p, dd2, d1s, d2s, d1p, d0s, d0p, Q(3)
   double precision, allocatable :: s0s(:), s1s(:), s2s(:), s3s(:), s4s(:)

   allocate(s0s(ntatom), s1s(ntatom), s2s(ntatom), s3s(ntatom), s4s(ntatom))

   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)
   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) / 2
   enddo

   ns = nshell(0); np = nshell(1); nd = nshell(2)
   M2 = 2 * M; M11 = 1 + 3 * M * (M +1) / 2 + Md * (Md +1)

   ! Nuclear attraction / repulsion between nuclei and classical atomic
   ! charges.
   Ens = 0.0D0
   do iatom = 1, natom
   do jatom = natom+1, ntatom
      term = dsqrt((r(iatom,1) - r(jatom,1)) * (r(iatom,1) - r(jatom,1)) + &
                   (r(iatom,2) - r(jatom,2)) * (r(iatom,2) - r(jatom,2)) + &
                   (r(iatom,3) - r(jatom,3)) * (r(iatom,3) - r(jatom,3)))
      Ens  = Ens + Iz(iatom) * pc(jatom) / term
   enddo
   enddo

   E1s = 0.0D0
   if (.not. elec) return

   ! (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            term  = 2.D0 * PI * exp(-rexp) / Zij
            tna   = 0.D0

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
               tna = tna - pc(iatom) * term * FUNCT(0,uf)
            enddo

            term = ccoef * tna
            vecmat_ind = ifunct + ((M2 - jfunct) * (jfunct -1)) / 2
            Hmat(vecmat_ind) = Hmat(vecmat_ind) + term
            E1s = E1s + Rho(vecmat_ind) * term
         endif
      enddo
      enddo
   enddo
   enddo

   ! (p|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1, ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            term  = 2.D0 * PI * exp(-rexp) / Zij
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
               s1s(iatom) = term * FUNCT(1,uf)
               s0s(iatom) = term * FUNCT(0,uf)
            enddo

            do l1 = 1, 3
               t1  = Q(l1) - r(Nuc(ifunct),l1)
               tna = 0.D0
               do iatom = natom+1, ntatom
                  tna  = tna - pc(iatom) * &
                         (t1 * s0s(iatom) - (Q(l1) - r(iatom,l1)) * s1s(iatom))
               enddo
               term = ccoef * tna

               vecmat_ind = ifunct + l1 -1 + ((M2 - jfunct) * (jfunct -1)) / 2
               Hmat(vecmat_ind) = Hmat(vecmat_ind) + term
               E1s = E1s + Rho(vecmat_ind) * term
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

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
            term  = 2.D0 * PI * exp(-rexp) / Zij
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
               s1s(iatom) = term * FUNCT(1,uf)
               s0s(iatom) = term * FUNCT(0,uf)
               s2s(iatom) = term * FUNCT(2,uf)
            enddo

            do iatom = natom+1, ntatom
               do l1 = 1, 3
                  t1 = Q(l1) - r(Nuc(ifunct),l1)
                  t2 = Q(l1) - r(iatom,l1)
                  p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  p1s = t1 * s1s(iatom) - t2 * s2s(iatom)

                  lij = 3
                  if (ifunct .eq. jfunct) lij = l1
                  do l2 = 1, lij
                     t1  = Q(l2) - r(Nuc(jfunct),l2)
                     t2  = Q(l2) - r(iatom,l2)
                     tna = t1 * p0s - t2 * p1s
                     if (l1 .eq. l2) tna = tna + (s0s(iatom) - s1s(iatom)) / Z2
                     tna = tna * pc(iatom)

                     vecmat_ind = ifunct + l1-1 + &
                                ((M2 - (jfunct + l2 -1)) * (jfunct + l2 -2)) / 2
                     term     = - tna * ccoef
                     Hmat(vecmat_ind) = Hmat(vecmat_ind) + term
                     E1s = E1s + Rho(vecmat_ind) * term
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! (d|s)
   do ifunct = ns+np+1, M, 6
   do jfunct = 1, ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            Z2    = 2.D0 * Zij
            term  = 2.D0 * PI * exp(-rexp) / Zij
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
               s1s(iatom) = term * FUNCT(1,uf)
               s0s(iatom) = term * FUNCT(0,uf)
               s2s(iatom) = term * FUNCT(2,uf)
            enddo

            do iatom = natom+1, ntatom
               do l1 = 1, 3
                  t1  = Q(l1) - r(Nuc(ifunct),l1)
                  t2  = Q(l1) - r(iatom,l1)
                  p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  p1s = t1 * s1s(iatom) - t2 * s2s(iatom)

                  do l2 = 1, l1
                     tna = (Q(l2) - r(Nuc(ifunct),l2)) * p0s &
                         - (Q(l2) - r(iatom,l2)      ) * p1s

                     f1 = 1.D0
                     if (l1 .eq. l2) then
                        tna = tna + (s0s(iatom) - s1s(iatom)) / Z2
                        f1  = sq3
                     endif
                     term = - tna * pc(iatom) * ccoef / f1

                     vecmat_ind = ifunct + Ll(l1) + l2 -1 + &
                                ((M2 - jfunct) * (jfunct -1)) / 2
                     Hmat(vecmat_ind) = Hmat(vecmat_ind) + term
                     E1s = E1s + Rho(vecmat_ind) * term
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! (d|p)
   do ifunct = ns+np+1, M    , 6
   do jfunct = ns+1   , ns+np, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rexp = a(ifunct,nci) * a(jfunct,ncj) * d(Nuc(ifunct),Nuc(jfunct)) / Zij

         if (rexp .lt. rmax) then
            ccoef = c(ifunct,nci) * c(jfunct,ncj)
            Z2    = 2.D0 * Zij
            term  = 2.D0 * PI * exp(-rexp) / Zij
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
               s1s(iatom) = term * FUNCT(1,uf)
               s0s(iatom) = term * FUNCT(0,uf)
               s2s(iatom) = term * FUNCT(2,uf)
               s3s(iatom) = term * FUNCT(3,uf)
            enddo

            do iatom = natom+1, ntatom
               do l1 = 1, 3
                  t1  = Q(l1) - r(Nuc(ifunct),l1)
                  t2  = Q(l1) - r(iatom,l1)
                  p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  p1s = t1 * s1s(iatom) - t2 * s2s(iatom)
                  p2s = t1 * s2s(iatom) - t2 * s3s(iatom)

                  do l2 = 1, l1
                     t1 = Q(l2) - r(Nuc(ifunct),l2)
                     t2 = Q(l2) - r(iatom,l2)
                     pj0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                     pj1s = t1 * s1s(iatom) - t2 * s2s(iatom)

                     f1  = 1.D0
                     d0s = t1 * p0s - t2 * p1s
                     d1s = t1 * p1s - t2 * p2s
                     if (l1 .eq. l2) then
                        f1  = sq3
                        d0s = d0s + (s0s(iatom) - s1s(iatom)) / Z2
                        d1s = d1s + (s1s(iatom) - s2s(iatom)) / Z2
                     endif

                     do l3 = 1, 3
                        tna = (Q(l3) - r(Nuc(jfunct),l3)) * d0s - &
                              (Q(l3) - r(iatom,l3)      ) * d1s

                        if (l1 .eq. l3) tna = tna + (pj0s - pj1s) / Z2
                        if (l2 .eq. l3) tna = tna + (p0s  - p1s ) / Z2
                        term = - tna * pc(iatom) * ccoef / f1

                        vecmat_ind = ifunct + Ll(l1) + l2 -1 + &
                                   ((M2 - (jfunct + l3-1)) * (jfunct + l3 -2))/2
                        Hmat(vecmat_ind) = Hmat(vecmat_ind) + term
                        E1s = E1s + Rho(vecmat_ind) * term
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

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
            term  = 2.D0 * PI * exp(-rexp) / Zij
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
               s1s(iatom) = term * FUNCT(1,uf)
               s0s(iatom) = term * FUNCT(0,uf)
               s2s(iatom) = term * FUNCT(2,uf)
               s3s(iatom) = term * FUNCT(3,uf)
               s4s(iatom) = term * FUNCT(4,uf)
            enddo

            do iatom = natom+1, ntatom
               do l1 = 1, 3
                  t1  = Q(l1) - r(Nuc(ifunct),l1)
                  t2  = Q(l1) - r(iatom,l1)
                  p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                  p1s = t1 * s1s(iatom) - t2 * s2s(iatom)
                  p2s = t1 * s2s(iatom) - t2 * s3s(iatom)
                  p3s = t1 * s3s(iatom) - t2 * s4s(iatom)

                  do l2 = 1, l1
                     f1 = 1.D0
                     t1 = Q(l2) - r(Nuc(ifunct),l2)
                     t2 = Q(l2) - r(iatom,l2)
                     pj0s = t1 * s0s(iatom) - t2 * s1s(iatom)
                     pj1s = t1 * s1s(iatom) - t2 * s2s(iatom)
                     pj2s = t1 * s2s(iatom) - t2 * s3s(iatom)
                     d0s = t1 * p0s - t2 * p1s
                     d1s = t1 * p1s - t2 * p2s
                     d2s = t1 * p2s - t2 * p3s

                     if (l1 .eq. l2) then
                        f1  = sq3
                        d0s = d0s + (s0s(iatom) - s1s(iatom)) / Z2
                        d1s = d1s + (s1s(iatom) - s2s(iatom)) / Z2
                        d2s = d2s + (s2s(iatom) - s3s(iatom)) / Z2
                     endif

                     lij = 3
                     if (ifunct .eq. jfunct) lij = l1
                     do l3 = 1, lij
                        t1 = Q(l3) - r(Nuc(jfunct),l3)
                        t2 = Q(l3) - r(iatom,l3)
                        d0p  = t1 * d0s  - t2 * d1s
                        d1p  = t1 * d1s  - t2 * d2s
                        pi0p = t1 * p0s  - t2 * p1s
                        pi1p = t1 * p1s  - t2 * p2s
                        pj0p = t1 * pj0s - t2 * pj1s
                        pj1p = t1 * pj1s - t2 * pj2s

                        if (l1 .eq. l3) then
                           d0p  = d0p  + (pj0s - pj1s) / Z2
                           d1p  = d1p  + (pj1s - pj2s) / Z2
                           pi0p = pi0p + (s0s(iatom) - s1s(iatom)) / Z2
                           pi1p = pi1p + (s1s(iatom) - s2s(iatom)) / Z2
                        endif
                        if (l2 .eq. l3) then
                           d0p  = d0p  + (p0s - p1s) / Z2
                           d1p  = d1p  + (p1s - p2s) / Z2
                           pj0p = pj0p + (s0s(iatom) - s1s(iatom)) / Z2
                           pj1p = pj1p + (s1s(iatom) - s2s(iatom)) / Z2
                        endif

                        lk = l3
                        if (ifunct .eq. jfunct) lk = min(l3, Ll(l1) -Ll(l3) +l2)
                        do l4 = 1, lk
                           f2  = 1.D0
                           t1  = Q(l4) - r(Nuc(jfunct),l4)
                           t2  = Q(l4) - r(iatom,l4)
                           tna = t1 * d0p - t2 * d1p

                           if (l4 .eq. l1) tna = tna + (pj0p - pj1p) / Z2
                           if (l4 .eq. l2) tna = tna + (pi0p - pi1p) / Z2
                           if (l4 .eq. l3) then
                              tna = tna + (d0s - d1s) / Z2
                              f2  = sq3
                           endif
                           term = - pc(iatom) * tna * ccoef / (f1 * f2)

                           vecmat_ind = ifunct + Ll(l1) + l2 -1 + &
                                       ((M2 - (jfunct + Ll(l3) + l4-1)) * &
                                       (jfunct + Ll(l3) + l4 -2)) / 2
                           Hmat(vecmat_ind) = Hmat(vecmat_ind) + term
                           E1s = E1s + Rho(vecmat_ind) * term
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

   deallocate(s0s, s1s, s2s, s3s, s4s)
   return
end subroutine
end module subm_intsol
