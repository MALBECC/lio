!------------------------------------------------------------------------------!
! Intfld subroutine evaluates the reaction field contribution to 1e matrix     !
! elements. The total dipole moment should be given in input (3 components).   !
! Integrals are evaluated using Obara Saika method.                            !
!                                                                              !
! Input: density basis data, dipole moment components.                         !
! Output: Fock matrix elements.                                                !
!                                                                              !
! Written by D. Estrin: 19/01/1993                                             !
! Refactored by F. Pedron: 08/2018                                             !
!------------------------------------------------------------------------------!
module subm_intfld
contains
subroutine intfld(Fmat, Fmat_b, r, d, Iz, natom, ntatom, open_shell, g, ux, uy,&
                  uz)
   use garcha_mod, only: a, c, nuc, ncont, nshell, NORM, M, Md
   use constants_mod, only: pi32
   implicit none
   integer         , intent(in)    :: natom, ntatom, Iz(natom)
   logical         , intent(in)    :: open_shell
   double precision, intent(in)    :: g, ux, uy, uz, r(ntatom,3), d(natom,natom)
   double precision, intent(inout) :: Fmat(:), Fmat_b(:)

   double precision :: sq3, term, ccoef, f1, f2
   double precision :: aux(3) , aux1(3), aux2(3), aux3(3), aux4(3), aux5(3), &
                       aux6(3), Q(3)

   integer          :: MM, M2, l1, l2, l3, l4, Ll(3), ns, np, nd, ifunct, &
                       jfunct, nci, ncj, i_ind, j_ind, fock_ind

   double precision :: ss, ps, pp, pis, pjs, dp, dijs, sxs, sys, szs, t1, ti, &
                       tj, Z2, Zij

   sq3 = 1.0D0
   if (NORM) sq3 = dsqrt(3.D0)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) / 2
   enddo

   ! M3 is Fock alpha, M5 is Fock beta
   M2 = 2 * M    ;
   ns = nshell(0); np = nshell(1)         ; nd = nshell(2)

   ! (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,nci) / Zij

         Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)

         ss  = exp(-a(ifunct,nci)*a(jfunct,ncj)*d(Nuc(ifunct),Nuc(jfunct))/Zij)&
             * pi32 / (Zij * sqrt(Zij))
         sxs = Q(1) * ss
         sys = Q(2) * ss
         szs = Q(3) * ss

         term     = g * ccoef * (sxs * ux + sys * uy + szs * uz)
         fock_ind = ifunct + ((M2 - jfunct) * (jfunct -1)) / 2
         Fmat(fock_ind) = Fmat(fock_ind) + term
         if (open_shell) Fmat_b(fock_ind) = Fmat_b(fock_ind) + term
      enddo
      enddo
   enddo
   enddo

   ! (p|s)
   do ifunct = ns+1,ns+np,3
   do jfunct = 1,ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,nci) / Zij

         Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)

         ss  = exp(-a(ifunct,nci)*a(jfunct,ncj)*d(Nuc(ifunct),Nuc(jfunct))/Zij)&
             * pi32 / (Zij * sqrt(Zij))
         sxs = Q(1) * ss
         sys = Q(2) * ss
         szs = Q(3) * ss

         ! l2: different p in the p shell ( x,y,z respectively)
         do l1 = 1, 3
            t1 = Q(l1) - r(Nuc(ifunct),l1)
            aux(1)  = t1 * sxs
            aux(2)  = t1 * sys
            aux(3)  = t1 * szs
            aux(l1) = aux(l1) + ss / Z2

            term     = g * ccoef * (aux(1)*ux + aux(2)*uy + aux(3)*uz)
            fock_ind = ifunct + l1-1 + ((M2 - jfunct) * (jfunct -1)) / 2
            Fmat(fock_ind) = Fmat(fock_ind) + term
            if (open_shell) Fmat_b(fock_ind) = Fmat_b(fock_ind) + term
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|p)
   do ifunct = ns+1, ns+np , 3
   do jfunct = ns+1, ifunct, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,nci) / Zij

         Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)

         ss  = exp(-a(ifunct,nci)*a(jfunct,ncj)*d(Nuc(ifunct),Nuc(jfunct))/Zij)&
             * pi32 / (Zij * sqrt(Zij))
         sxs = Q(1) * ss
         sys = Q(2) * ss
         szs = Q(3) * ss

         do l1 = 1, 3
            t1 = Q(l1) - r(Nuc(ifunct),l1)
            ps = ss * t1
            aux(1)  = t1 * sxs
            aux(2)  = t1 * sys
            aux(3)  = t1 * szs
            aux(l1) = aux(l1) + ss / Z2

            do l2 = 1, 3
               t1 = Q(l2) - r(Nuc(jfunct),l2)
               pp = t1 * ps
               aux1(1) = aux(1) * t1
               aux1(2) = aux(2) * t1
               aux1(3) = aux(3) * t1

               if (l1 .eq. l2) then
                  aux1(1) = aux1(1) + sxs / Z2
                  aux1(2) = aux1(2) + sys / Z2
                  aux1(3) = aux1(3) + szs / Z2
                  pp = pp + ss / Z2
               endif
               aux1(l2) = aux1(l2) + ps / Z2

               i_ind = ifunct + l1 -1
               j_ind = jfunct + l2 -1
               if (i_ind .ge. j_ind) then
                  fock_ind = i_ind + ((M2 - j_ind) * (j_ind -1)) / 2
                  term     = g * ccoef * (aux1(1)*ux + aux1(2)*uy + aux1(3)*uz)

                  Fmat(fock_ind) = Fmat(fock_ind) + term
                  if (open_shell) Fmat_b(fock_ind) = Fmat_b(fock_ind) + term
               endif
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|s)
   do ifunct = ns+np+1,M,6
   do jfunct = 1,ns
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,nci) / Zij

         Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)

         ss  = exp(-a(ifunct,nci)*a(jfunct,ncj)*d(Nuc(ifunct),Nuc(jfunct))/Zij)&
            * pi32 / (Zij * sqrt(Zij))
         sxs = Q(1) * ss
         sys = Q(2) * ss
         szs = Q(3) * ss

         do l1 = 1, 3
            t1 = Q(l1) - r(Nuc(ifunct),l1)
            ps = ss * t1
            aux(1)  = t1 * sxs
            aux(2)  = t1 * sys
            aux(3)  = t1 * szs
            aux(l1) = aux(l1) + ss / Z2

            do l2 = 1, l1
               t1 = Q(l2) - r(Nuc(ifunct),l2)
               f1 = 1.0D0
               aux1(1) = aux(1) * t1
               aux1(2) = aux(2) * t1
               aux1(3) = aux(3) * t1

               if (l1 .eq. l2) then
                  aux1(1) = aux1(1) + sxs / Z2
                  aux1(2) = aux1(2) + sys / Z2
                  aux1(3) = aux1(3) + szs / Z2
                  f1 = sq3
               endif
               aux1(l2) = aux1(l2) + ps / Z2
               term     = g * ccoef * (aux1(1)*ux + aux1(2)*uy + aux1(3)*uz) /f1

               fock_ind = ifunct + Ll(l1) + l2 -1 + ((M2-jfunct)*(jfunct -1)) /2
               Fmat(fock_ind) = Fmat(fock_ind) + term
               if (open_shell) Fmat_b(fock_ind) = Fmat_b(fock_ind) + term
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|p)
   do ifunct = ns+np+1,M,6
   do jfunct = ns+1,ns+np,3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,nci) / Zij

         Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)

         ss  = exp(-a(ifunct,nci)*a(jfunct,ncj)*d(Nuc(ifunct),Nuc(jfunct))/Zij)&
             * pi32 / (Zij * sqrt(Zij))
         sxs = Q(1) * ss
         sys = Q(2) * ss
         szs = Q(3) * ss

         do l1 = 1, 3 ! aux: (pi|r|s)
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            pis = ss * t1
            aux(1)  = t1 * sxs
            aux(2)  = t1 * sys
            aux(3)  = t1 * szs
            aux(l1) = aux(l1) + ss / Z2

            do l2 = 1, l1 ! aux1: (pj|r|s), aux2 : (dij|r|s)
               t1   = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = t1 * ss
               dijs = t1 * pis
               f1   = 1.0D0
               aux1(1)  = t1 * sxs
               aux1(2)  = t1 * sys
               aux1(3)  = t1 * szs
               aux1(l2) = aux1(l2) + ss / Z2
               aux2(1)  = t1 * aux(1)
               aux2(2)  = t1 * aux(2)
               aux2(3)  = t1 * aux(3)
               aux2(l2) = aux2(l2) + pis / Z2

               if (l1 .eq. l2) then
                  f1   = sq3
                  dijs = dijs + ss / Z2
                  aux2(1) = aux2(1) + sxs / Z2
                  aux2(2) = aux2(2) + sys / Z2
                  aux2(3) = aux2(3) + szs / Z2
               endif

               do l3 = 1, 3
                  t1 = Q(l3) - r(Nuc(jfunct),l3)
                  aux3(1) = aux2(1) * t1
                  aux3(2) = aux2(2) * t1
                  aux3(3) = aux2(3) * t1
                  if (l1 .eq. l3) then
                     aux3(1) = aux3(1) + aux1(1) / Z2
                     aux3(2) = aux3(2) + aux1(2) / Z2
                     aux3(3) = aux3(3) + aux1(3) / Z2
                  endif
                  if (l2 .eq. l3) then
                     aux3(1) = aux3(1) + aux(1) / Z2
                     aux3(2) = aux3(2) + aux(2) / Z2
                     aux3(3) = aux3(3) + aux(3) / Z2
                  endif
                  aux3(l3) = aux3(l3) + dijs / Z2
                  term     = g*ccoef * (aux3(1)*ux + aux3(2)*uy + aux3(3)*uz)/f1

                  fock_ind = ifunct + Ll(l1) + l2 -1 + &
                             ((M2 - jfunct + l3 -1) * (jfunct + l3 -2)) / 2
                  Fmat(fock_ind) = Fmat(fock_ind) + term
                  if (open_shell) Fmat_b(fock_ind) = Fmat_b(fock_ind) + term
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|d)
   do ifunct = ns+np+1,M,6
   do jfunct = ns+np+1,M,6
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,nci) / Zij

         Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)

         ss  = exp(-a(ifunct,nci)*a(jfunct,ncj)*d(Nuc(ifunct),Nuc(jfunct))/Zij)&
             * pi32 / (Zij * sqrt(Zij))
         sxs = Q(1) * ss
         sys = Q(2) * ss
         szs = Q(3) * ss

         do l1 = 1, 3 ! aux: (pi|r|s)
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            pis = ss * t1
            aux(1)  = t1 * sxs
            aux(2)  = t1 * sys
            aux(3)  = t1 * szs
            aux(l1) = aux(l1) + ss / Z2

            do l2 = 1, l1 ! aux1: (pj|r|s), aux2: (dij|r|s)
               t1   = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = ss * t1
               dijs = t1 * pis
               f1   = 1.D0

               aux1(1)  = t1 * sxs
               aux1(2)  = t1 * sys
               aux1(3)  = t1 * szs
               aux1(l2) = aux1(l2) + ss / Z2
               aux2(1)  = t1 * aux(1)
               aux2(2)  = t1 * aux(2)
               aux2(3)  = t1 * aux(3)
               aux2(l2) = aux2(l2)+pis / Z2

               if (l1 .eq. l2) then
                  f1 = sq3
                  aux2(1) = aux2(1) + sxs / Z2
                  aux2(2) = aux2(2) + sys / Z2
                  aux2(3) = aux2(3) + szs / Z2
                  dijs = dijs + ss / Z2
               endif

               do l3 = 1, 3 ! aux3: (dij|r|pk), aux4: (pi|r|pk), aux5: (pj|r|pk)
                  t1 = Q(l3) - r(Nuc(jfunct),l3)
                  dp = t1 * dijs
                  aux3(1) = aux2(1) * t1
                  aux3(2) = aux2(2) * t1
                  aux3(3) = aux2(3) * t1
                  aux4(1) = aux(1)  * t1
                  aux4(2) = aux(2)  * t1
                  aux4(3) = aux(3)  * t1
                  aux5(1) = aux1(1) * t1
                  aux5(2) = aux1(2) * t1
                  aux5(3) = aux1(3) * t1

                  if (l1 .eq. l3) then
                     dp = dp + pjs / Z2
                     aux3(1) = aux3(1) + aux1(1) / Z2
                     aux3(2) = aux3(2) + aux1(2) / Z2
                     aux3(3) = aux3(3) + aux1(3) / Z2
                     aux4(1) = aux4(1) + sxs / Z2
                     aux4(2) = aux4(2) + sys / Z2
                     aux4(3) = aux4(3) + szs / Z2
                  endif
                  if (l2 .eq. l3) then
                     dp = dp + pis / Z2
                     aux3(1) = aux3(1) + aux(1) / Z2
                     aux3(2) = aux3(2) + aux(2) / Z2
                     aux3(3) = aux3(3) + aux(3) / Z2
                     aux5(1) = aux5(1) + sxs / Z2
                     aux5(2) = aux5(2) + sys / Z2
                     aux5(3) = aux5(3) + szs / Z2
                  endif
                  aux3(l3) = aux3(l3)+dijs / Z2
                  aux4(l3) = aux4(l3)+pis / Z2
                  aux5(l3) = aux5(l3)+pjs / Z2

                  do l4 = 1, l3 ! aux6: (d|r|d)
                     t1 = Q(l4) - r(Nuc(jfunct),l4)
                     f2 = 1.D0
                     aux6(1) = aux3(1) * t1
                     aux6(2) = aux3(2) * t1
                     aux6(3) = aux3(3) * t1

                     if (l1 .eq. l4) then
                        aux6(1) = aux6(1) + aux5(1) / Z2
                        aux6(2) = aux6(2) + aux5(2) / Z2
                        aux6(3) = aux6(3) + aux5(3) / Z2
                     endif
                     if (l2 .eq. l4) then
                        aux6(1) = aux6(1) + aux4(1) / Z2
                        aux6(2) = aux6(2) + aux4(2) / Z2
                        aux6(3) = aux6(3) + aux4(3) / Z2
                     endif
                     if (l3 .eq. l4) then
                        f2 = sq3
                        aux6(1) = aux6(1) + aux2(1) / Z2
                        aux6(2) = aux6(2) + aux2(2) / Z2
                        aux6(3) = aux6(3) + aux2(3) / Z2
                     endif
                     aux6(l4) = aux6(l4) + dp / Z2

                     i_ind = ifunct + Ll(l1) + l2 -1
                     j_ind = jfunct + Ll(l3) + l4 -1
                     if (i_ind .ge. j_ind) then
                        fock_ind = i_ind + ((M2-j_ind)*(j_ind-1))/2
                        term     = (aux6(1) * ux + aux6(2) * uy + aux6(3) * uz)&
                                 * g* ccoef / (f1 * f2)

                        Fmat(fock_ind) = Fmat(fock_ind) + term
                        if (open_shell) Fmat_b(fock_ind) = Fmat_b(fock_ind) &
                                                       + term
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   return
end subroutine intfld
end module subm_intfld
