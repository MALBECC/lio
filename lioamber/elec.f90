!##############################################################################!
!## ELEC ######################################################################!
! Calculation of electrical potential for an arbitrary point of space, making  !
! 1e integrals via the Obara-Saika recursive method.                           !
!                                                                              !
! Inputs:                                                                      !
!   · NX, NY, NZ      : Number of x, y and z points.                           !
!   · xMin, yMin, zMin: Minimum value for x, y and z coordinates.              !
!   · deltax          : Increments for x, y and z coordinates (same for all).  !
!                                                                              !
! Outputs:                                                                     !
!   · A file containing the electrostatic potential surface.                   !
!                                                                              !
!##############################################################################!
subroutine elec(NX, NY, NZ, deltax, xMin, yMin, zMin)
   use garcha_mod   , only: r, d, natom, cube_elec_file, RMM, Iz
   use constants_mod, only: PI, PI32
   use basis_data   , only: M, norm, nShell, nCont, nuc, a, c
   use liosubs_math , only: funct

   implicit none
   double precision, intent(in) :: xMin, yMin, zMin, deltax
   integer         , intent(in) :: NX, NY, NZ

   integer          :: iatom, ifunct, i_ind, jatom, jfunct, j_ind, l1, l2, l3, &
                       l4, lij, lk, M2, MM, nnx, nny, nnz, ntotal, ns, np, nd, &
                       nci, ncj, rho_ind

   double precision :: Q(3), xi(3)
   double precision :: V_nuc, tna, temp, SQ3, Z2, Zij, ccoef, rExp, uf, t1, t2,&
                       f1, f2, s0s, s1s, s2s, s3s, s4s, pj0s, pj0p, pj1s, pj1p,&
                       p0s, p1s, p2s, p3s, PI0p, PI1p, pj2s, d0s, d0p, d1s,    &
                       d1p, d2s

   double precision, allocatable :: pote(:)
   double precision, parameter   :: RMAX = 3.0D0
   integer         , parameter   :: LL(3) = (/0, 1, 3/) ! LL is l * (l -1) / 2

   ! Variable initialization.
   SQ3 = 1.0D0
   if (norm) SQ3 = dsqrt(3.D0)

   ns = nShell(0); np = nShell(1); nd = nShell(2)
   MM = M * (M +1) / 2 ; M2 = 2 * M

   ! Calculates squared distance.
   do iatom = 1, natom
   do jatom = 1, natom
      d(iatom,jatom) = (r(iatom,1) - r(jatom,1)) * (r(iatom,1) - r(jatom,1)) + &
                       (r(iatom,2) - r(jatom,2)) * (r(iatom,2) - r(jatom,2)) + &
                       (r(iatom,3) - r(jatom,3)) * (r(iatom,3) - r(jatom,3))
   enddo
   enddo

   ntotal = NX * NY * NZ
   allocate(pote(ntotal))
   pote = 0.0D0

   ! First loop - (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij
            rho_ind = ifunct + (M2 - jfunct) * (jfunct-1) / 2

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal +1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax
               uf    = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                        (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                        (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij
               s0s   = temp * FUNCT(0,uf)

               pote(ntotal) = pote(ntotal) + ccoef * RMM(rho_ind) * s0s
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Second loop - (p|s)
   do ifunct = ns +1, ns + np, 3
   do jfunct = 1, ns
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal +1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij
               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)

               do l2 = 1, 3
                  tna     = (Q(l2) - r(nuc(ifunct),l2)) * s0s - &
                            (Q(l2) - xi(l2)           ) * s1s
                  rho_ind = ifunct + l2 -1 + ((M2 - jfunct) * (jfunct -1)) / 2

                  pote(ntotal) = pote(ntotal) + ccoef * tna * RMM(rho_ind)
               enddo
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Third loop - (p|p)
   do ifunct = ns+1, ns + np, 3
   do jfunct = ns+1, ifunct , 3
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)

         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal +1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s

                  lij = 3
                  if (ifunct .eq. jfunct) lij = l1
                  do l2 = 1, lij
                     t1 = Q(l2) - r(nuc(jfunct),l2)
                     t2 = Q(l2) - xi(l2)

                     tna = (Q(l2) - r(nuc(jfunct),l2)) * p0s - &
                           (Q(l2) - xi(l2))            * p1s
                     if (l1 .eq. l2) tna = tna + (s0s - s1s) / Z2

                     i_ind   = ifunct + l1 -1
                     j_ind   = jfunct + l2 -1
                     rho_ind = i_ind + ((M2 - j_ind) * (j_ind -1)) / 2

                     pote(ntotal) = pote(ntotal) + tna * ccoef * RMM(rho_ind)
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

   ! Fourth loop - (d|s)
   do ifunct = ns + np +1, M, 6
   do jfunct = 1, ns
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal + 1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s

                  do l2 = 1, l1
                     t1  = Q(l2) - r(nuc(ifunct),l2)
                     t2  = Q(l2) - xi(l2)
                     tna = t1 * p0s - t2 * p1s

                     f1 = 1.0D0
                     if (l1 .eq. l2) then
                        tna = tna + (s0s - s1s)   / Z2
                        f1  = SQ3
                     endif

                     i_ind   = ifunct + ll(l1) + l2 -1
                     rho_ind = i_ind + ((M2 - jfunct) * (jfunct -1)) / 2

                     pote(ntotal) = pote(ntotal) + tna * RMM(rho_ind) * ccoef/f1
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

   ! Fifth loop - (d|p)
   do ifunct = ns + np +1, M      , 6
   do jfunct = ns +1     , ns + np, 3
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal + 1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)
               s3s = temp * FUNCT(3,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s
                  p2s = t1 * s2s - t2 * s3s

                  do l2 = 1, l1
                     t1   = Q(l2) - r(nuc(ifunct),l2)
                     t2   = Q(l2) - xi(l2)
                     pj0s = t1 * s0s - t2 * s1s
                     pj1s = t1 * s1s - t2 * s2s
                     d0s  = t1 * p0s - t2 * p1s
                     d1s  = t1 * p1s - t2 * p2s
                     f1   = 1.0D0

                     if (l1 .eq. l2) then
                        f1  = SQ3
                        d0s = d0s + (s0s - s1s) / Z2
                        d1s = d1s + (s1s - s2s) / Z2
                     endif

                     do l3 = 1, 3
                        t1  = Q(l3) - r(nuc(jfunct),l3)
                        t2  = Q(l3) - xi(l3)
                        tna = t1 * d0s - t2 * d1s

                        if (l1 .eq. l3) tna = tna + (pj0s - pj1s) / Z2
                        if (l2 .eq. l3) tna = tna + (p0s  - p1s ) / Z2

                        i_ind   = ifunct + ll(l1) + l2 -1
                        j_ind   = jfunct + l3 -1
                        rho_ind = i_ind + ((M2 - j_ind) * (j_ind -1)) / 2

                        pote(ntotal) = pote(ntotal) + tna * RMM(rho_ind) * &
                                       ccoef / f1
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

   ! Sixth loop - (d|d)
   do ifunct = ns + np +1, M     , 6
   do jfunct = ns + np +1, ifunct, 6
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal  = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal + 1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)
               s3s = temp * FUNCT(3,uf)
               s4s = temp * FUNCT(4,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s
                  p2s = t1 * s2s - t2 * s3s
                  p3s = t1 * s3s - t2 * s4s

                  do l2 = 1, l1
                     t1   = Q(l2) - r(nuc(ifunct),l2)
                     t2   = Q(l2) - xi(l2)
                     pj0s = t1 * s0s - t2 * s1s
                     pj1s = t1 * s1s - t2 * s2s
                     pj2s = t1 * s2s - t2 * s3s
                     d0s  = t1 * p0s - t2 * p1s
                     d1s  = t1 * p1s - t2 * p2s
                     d2s  = t1 * p2s - t2 * p3s
                     f1   = 1.D0

                     if (l1 .eq. l2) then
                        f1  = SQ3
                        d0s = d0s + (s0s - s1s) / Z2
                        d1s = d1s + (s1s - s2s) / Z2
                        d2s = d2s + (s2s - s3s) / Z2
                     endif

                     lij = 3
                     if (ifunct .eq. jfunct) lij = l1

                     do l3 = 1, lij
                        t1   = Q(l3) - r(nuc(jfunct),l3)
                        t2   = Q(l3) - xi(l3)
                        d0p  = t1 * d0s - t2 * d1s
                        d1p  = t1 * d1s - t2 * d2s
                        PI0p = t1 * p0s - t2 * p1s
                        PI1p = t1 * p1s - t2 * p2s
                        pj0p = t1 * pj0s - t2 * pj1s
                        pj1p = t1 * pj1s - t2 * pj2s

                        if (l1 .eq. l3) then
                           d0p  = d0p  + (pj0s - pj1s) / Z2
                           d1p  = d1p  + (pj1s - pj2s) / Z2
                           PI0p = PI0p + (s0s  - s1s ) / Z2
                           PI1p = PI1p + (s1s  - s2s ) / Z2
                        endif
                        if (l2.eq.l3) then
                           d0p  = d0p  + (p0s - p1s) / Z2
                           d1p  = d1p  + (p1s - p2s) / Z2
                           pj0p = pj0p + (s0s - s1s) / Z2
                           pj1p = pj1p + (s1s - s2s) / Z2
                        endif

                        lk = l3
                        if (ifunct .eq. jfunct) lk = min(l3, Ll(l1) - Ll(l3)+l2)
                        do l4 = 1, lk
                           f2 = 1.0D0
                           tna = (Q(l4) - R(nuc(jfunct),l4)) * d0p - &
                                 (Q(l4) - xi(l4)           ) * d1p

                           if (l4 .eq. l1) tna = tna + (pj0p - pj1p) / Z2
                           if (l4 .eq. l2) tna = tna + (PI0p - PI1p) / Z2
                           if (l4 .eq. l3) then
                              f2  = SQ3
                              tna = tna + (d0s - d1s) / Z2
                           endif

                           i_ind   = ifunct + ll(l1) + l2 -1
                           j_ind   = jfunct + ll(l3) + l4 -1
                           rho_ind = i_ind + ((M2 - j_ind) * (j_ind -1)) / 2

                           pote(ntotal) = pote(ntotal) + RMM(rho_ind) * tna * &
                                          ccoef / (f1 * f2)
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

   ! Electron potential is finished, now core potential is calculated.
   ntotal = 0
   do nnx = 1, NX
   do nny = 1, NY
   do nnz = 1, NZ
      ntotal = ntotal + 1

      xi(1) = xMin + (nnx -1) * deltax
      xi(2) = yMin + (nny -1) * deltax
      xi(3) = zMin + (nnz -1) * deltax

      V_nuc = 0.0D0
      do iatom = 1, natom
         V_nuc = V_nuc - Iz(iatom) / &
                 sqrt( (xi(1) - r(iatom,1)) * (xi(1) - r(iatom,1)) + &
                       (xi(2) - r(iatom,2)) * (xi(2) - r(iatom,2)) + &
                       (xi(3) - r(iatom,3)) * (xi(3) - r(iatom,3)) )
      enddo
      pote(ntotal) = pote(ntotal) + V_nuc
   enddo
   enddo
   enddo

   ! Writes a cubefile similar to those created by GAUSSIAN or ORCA, with the
   ! following format (as of GAUSSIAN98); all coordinates are in atomic units:
   !     LINE   FORMAT      CONTENTS
   !   ===============================================================
   !   1     A           TITLE
   !   2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
   !   3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
   !   4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
   !   #ATOMS LINES OF ATOM COORDINATES:
   !   ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
   !   REST: 6E13.5      CUBE DATA (WITH Z INCREMENT MOVING FASTEST, THEN
   !                     Y AND THEN X)
   !
   ! For orbital cube files, #ATOMS will be negative (< 0) with one additional
   ! line after the final atom containing the number of orbitales and their
   ! respective numbers. The orbital number will also be the fastest moving
   ! increment.
   open(unit = 188, file = cube_elec_file)

   write(188,*) 'Elfield'
   write(188,*)
   write(188,678) natom, xMin, yMin, zMin
   write(188,678) NX, deltax, 0.0D0 , 0.0D0
   write(188,678) NY, 0.0D0 , deltax, 0.0D0
   write(188,678) NZ, 0.0D0 , 0.0D0 , deltax
   do iatom = 1, natom
      write(188,677) Iz(iatom), 0.0D0, r(iatom,1), r(iatom,2), r(iatom,3)
   enddo
   write(188,679) pote

   close(188)
   deallocate(pote)

677 format(I5,4(F12.6))
678 format(I5,3(F12.6))
679 format(6(E13.5))
end subroutine elec
!##############################################################################!
