!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% ML subroutines%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! In this module there are subroutines useful for the construction of a        !
! machine-learning-based kinetic energy functional.                            !
! It is based in the subm_int1 module                                          !
!                                                                              !
! EXTERNAL INPUT: system information.                                          !
!   · natom: number of QM atoms.                                               !
!   · ntatom: total number of atoms (QM+MM)                                    !
!   · r(ntatom,3): atoms' coordinates.                                         !
!   · d(natom,natom): distances between QM atoms.                              !
!   · Iz(natom): nuclear charge for each QM atom.                              !
!                                                                              !
! INTERNAL INPUT: basis set information.                                       !
!   · Md: number of auxiliary basis functions (without contractions)           !
!   · ncontd(Md): number of contractions per auxiliary function.               !
!   · ad(M,nl): auxiliary basis function exponents.                            !
!   · cd(M,nl): auxiliary basis function coefficients.                         !
!   · nshelld(0:3): number of auxiliary basis functions per shell (s,p,d).     !
!   · Nucd(Md): atomic index corresponding to auxiliary function i.            !
!   · NORM: use custom normalization (now default and deprecated option)       !
!                                                                              !
! OUTPUTS:                                                                     !
!   · Smatd(Md,Md): the overlap matrix of the auxiliary basis set              !                                 
!                                                                              !
! Original code:  Dario Estrin, Feb/1992                                       !
! First refactor: Diego Armiño, May/2018                                       !
! Last refactor:  Federico Pedron, Aug/2018                                    !
! And this nasty module was made up by: El Barto                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "datatypes/datatypes.fh"
module ml_sub
   contains

subroutine ovlp_aux(Smatd, d, r, Iz, natom, ntatom )
   use basis_data   , only: Nucd, ad, cd, ncontd, NORM, Md, nshelld
   use liosubs_math , only: FUNCT
   use constants_mod, only: pi, pi32
   use fstsh_data   , only: FSTSH, Sovl_now
   implicit none

   integer,          intent(in) :: natom, ntatom, Iz(natom)
   LIODBLE, intent(in) :: d(natom,natom), r(ntatom,3)

   integer           :: my_natom, igpu, i_ind, j_ind, k_ind, ifunct, jfunct, &
                        iatom, jatom, nci, ncj, l1, l2, l3, l4, l12, l34,    &
                        MMd, nsd, npd, M2d
   LIODBLE  :: ovlap, uf, cc_f, Q(3), temp, sq3, alf, alf2, ccoef,    &
                        t0, t1, t2, f1, f2, tn, tna, Z2, Zij, ss, ps, dd, sks, &
                        p0s, p1s, p2s, p3s, pi0p, pi1p, piks, pikpk, pipk,     &
                        pis, pj0s, pj1s, pj2s, pj0p, pj1p, pjkpk, pjks, pjpk,  &
                        pjs, pks, dijs, dijpk, dijks, dijkpk, d0s, d0p, d1p,   &
                        d1s, d2s
   LIODBLE, allocatable, dimension(:) :: s0s, s1s, s2s, s3s, s4s
   LIODBLE, allocatable :: Smatd(:,:)                  ! Auxiliary basis set overlap matrix

   ! Allocates internal variables.
   allocate(s0s(natom),s1s(natom),s2s(natom), s3s(natom), s4s(natom))
   if (.not.allocated(Smatd)) allocate(Smatd(Md,Md))

   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)

   nsd  = nshelld(0); npd = nshelld(1)
   MMd  = Md*(Md+1)/2
   M2d  = 2*Md

   Smatd = 0.0D0

   ! First loop - (s|s)
   do ifunct = 1, nsd
   do jfunct = 1, ifunct
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij  = ad(ifunct,nci) + ad(jfunct,ncj)
         alf  = ad(ifunct,nci) * ad(jfunct,ncj) / Zij
         alf2 = alf * 2.D0
         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         ovlap = ss

         k_ind = ifunct + ((M2d-jfunct)*(jfunct-1))/2
         Smatd(ifunct,jfunct) = Smatd(ifunct,jfunct) + ccoef * ovlap
         Smatd(jfunct,ifunct) = Smatd(jfunct,ifunct) + ccoef * ovlap

      enddo
      enddo
   enddo
   enddo

   ! Second loop - (p|s)
   do ifunct = nsd+1, nsd+npd, 3
   do jfunct = 1, nsd
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct, nci) + ad(jfunct, ncj)
         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         alf   = ad(ifunct, nci) * ad(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         ! Loops over nuclei, common for all shells
         do iatom = 1, my_natom
            uf = Zij * ((Q(1) - r(iatom,1)) * (Q(1) - r(iatom,1)) + &
                        (Q(2) - r(iatom,2)) * (Q(2) - r(iatom,2)) + &
                        (Q(3) - r(iatom,3)) * (Q(3) - r(iatom,3)))
            temp = 2.D0 * sqrt(Zij/pi) * ss

            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
         enddo

         ! l2: different p in the p shell ( x,y,z respectively)
         do l2 = 1, 3
            t1 = Q(l2) - r(Nucd(ifunct),l2)
            ovlap = t1 * ss

            i_ind = ifunct + l2 -1
            k_ind = i_ind  + ((M2d-jfunct)*(jfunct-1))/2
            Smatd(i_ind ,jfunct) = Smatd(i_ind ,jfunct) + ovlap * ccoef
            Smatd(jfunct, i_ind) = Smatd(jfunct, i_ind) + ovlap * ccoef
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Third loop - (p|p)
   do ifunct = nsd+1, nsd+npd , 3
   do jfunct = nsd+1, ifunct, 3
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct, nci) + ad(jfunct, ncj)
         Z2  = 2.D0 * Zij
         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         alf   = ad(ifunct, nci) * ad(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         ! Loops over nuclei, common for all shells
         do iatom = 1, my_natom
            uf = Zij * ((Q(1) - r(iatom,1)) * (Q(1) - r(iatom,1)) + &
                        (Q(2) - r(iatom,2)) * (Q(2) - r(iatom,2)) + &
                        (Q(3) - r(iatom,3)) * (Q(3) - r(iatom,3)))
            temp = 2.D0 * sqrt(Zij/pi) * ss

            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
         enddo

         do l1 = 1, 3
            t1  = Q(l1) - r(Nucd(ifunct),l1)
            ps  = ss * t1
            pks = sks * t1 + alf2 * ps

            do l2 = 1, 3
               t1 = Q(l2) - r(Nucd(jfunct),l2)

               ovlap = t1 * ps
               tn    = t1 * pks + alf2 * ovlap
               if (l1 .eq. l2) then
                  ovlap = ovlap + ss / Z2
                  tn    = tn + (sks + alf2*ss) / Z2
               endif

               i_ind = ifunct + l1 -1
               j_ind = jfunct + l2 -1
               if (i_ind .ge. j_ind) then
                 k_ind = i_ind + ((M2d-j_ind)*(j_ind-1))/2
                 Smatd(i_ind,j_ind) = Smatd(i_ind,j_ind) + ovlap * ccoef
                 Smatd(j_ind,i_ind) = Smatd(j_ind,i_ind) + ovlap * ccoef
               endif
            enddo
         enddo

         ! Loops over nuclei, specific part
         do iatom = 1, my_natom
            do l1 = 1, 3
               t1  = Q(l1) - r(Nucd(ifunct),l1)
               t2  = Q(l1) - r(iatom,l1)
               p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
               p1s = t1 * s1s(iatom) - t2 * s2s(iatom)

               do l2 = 1, 3
                  t1  = Q(l2) - r(Nucd(jfunct),l2)
                  t2  = Q(l2) - r(iatom,l2)
                  tna = t1 * p0s - t2 * p1s

                  if (l1 .eq. l2) then
                     tna = tna + (s0s(iatom) - s1s(iatom)) / Z2
                  endif

                  i_ind = ifunct + l1 -1
                  j_ind = jfunct + l2 -1
                  if (i_ind .ge. j_ind) then
                     k_ind = i_ind + ((M2d-j_ind)*(j_ind-1))/2
                  endif
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fourth loop - (d|s)
   do ifunct = nsd+npd+1, Md, 6
   do jfunct = 1, nsd
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij  = ad(ifunct, nci) + ad(jfunct, ncj)
         Z2   = 2.D0 * Zij
         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         alf   = ad(ifunct, nci) * ad(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         ! Loops over nuclei, common for all shells
         do iatom = 1, my_natom
            uf = Zij * ((Q(1) - r(iatom,1)) * (Q(1) - r(iatom,1)) + &
                        (Q(2) - r(iatom,2)) * (Q(2) - r(iatom,2)) + &
                        (Q(3) - r(iatom,3)) * (Q(3) - r(iatom,3)))

            temp = 2.D0 * sqrt(Zij/pi) * ss
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
         enddo

         do l1 = 1, 3
            t1  = Q(l1) - r(Nucd(ifunct),l1)
            ps  = ss  * t1
            pks = sks * t1 + alf2 * ps

            do l2 = 1, l1
               t1    = Q(l2) - r(Nucd(ifunct),l2)
               ovlap = t1 * ps
               tn    = t1 * pks

               f1 = 1.D0
               if (l1 .eq. l2) then
                  ovlap = ovlap + ss / Z2
                  tn    = tn    + sks/ Z2 - alf2 * ss / (2.D0*ad(ifunct,nci))
                  f1    = sq3
               endif
               tn = tn + alf2 * ovlap

               ! Ordering of the d shell: xx, yx, yy, zx, zy, zz
               ! (11, 21, 22, 31, 32, 33)
               l12   = l1 * (l1-1)/2 + l2
               i_ind = ifunct + l12 -1
               k_ind = i_ind  + ((M2d-jfunct)*(jfunct-1))/2

               cc_f = ccoef / f1
               Smatd(i_ind,jfunct) = Smatd(i_ind,jfunct)  + ovlap * cc_f
               Smatd(jfunct,i_ind) = Smatd(jfunct,i_ind)  + ovlap * cc_f
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fifth loop - (d|p)
   do ifunct = nsd+npd+1, Md    , 6
   do jfunct = nsd+1   , nsd+npd, 3
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij  = ad(ifunct, nci) + ad(jfunct, ncj)
         Z2   = 2.D0 * Zij
         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         alf   = ad(ifunct, nci) * ad(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         ! Loops over nuclei, common for all shells
         do iatom = 1, my_natom
            uf = Zij * ((Q(1) - r(iatom,1)) * (Q(1) - r(iatom,1)) + &
                        (Q(2) - r(iatom,2)) * (Q(2) - r(iatom,2)) + &
                        (Q(3) - r(iatom,3)) * (Q(3) - r(iatom,3)))

            temp = 2.D0 * sqrt(Zij/pi) * ss
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
            s3s(iatom) = temp * FUNCT(3,uf)
         enddo

         do l1 = 1, 3
            t1   = Q(l1) - r(Nucd(ifunct),l1)
            pis  = ss  * t1
            piks = sks * t1 + alf2 * pis

            do l2 = 1, l1
               t1    = Q(l2) - r(Nucd(ifunct),l2)
               pjs   = ss  * t1
               pjks  = sks * t1 + alf2 * pjs
               dijs  = t1 * pis
               dijks = t1 * piks
               f1    = 1.D0

               if (l1 .eq. l2) then
                  f1    = sq3
                  dijs  = dijs  + ss / Z2
                  dijks = dijks + sks/ Z2 - alf2 * ss / (2.D0*ad(ifunct,nci))
               endif
               dijks = dijks + alf2 * dijs

               do l3 = 1, 3
                  t1    = Q(l3) - r(Nucd(jfunct),l3)
                  ovlap = t1 * dijs
                  tn    = t1 * dijks

                  if (l1 .eq. l3) then
                     ovlap = ovlap + pjs  / Z2
                     tn    = tn    + pjks / Z2
                  endif
                  if (l2 .eq. l3) then
                     ovlap = ovlap + pis  / Z2
                     tn    = tn    + piks / Z2
                  endif
                  tn   = tn + alf2 * ovlap
                  cc_f = ccoef / f1

                  l12 = l1 * (l1-1)/2 + l2
                  i_ind = ifunct + l12 -1
                  j_ind = jfunct + l3  -1
                  k_ind = i_ind + ((M2d-j_ind)*(j_ind-1))/2

                  Smatd(i_ind,j_ind) = Smatd(i_ind,j_ind) + cc_f * ovlap
                  Smatd(j_ind,i_ind) = Smatd(j_ind,i_ind) + cc_f * ovlap
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Sixth and final loop - (d|d)
   do ifunct = nsd+npd+1, Md     , 6
   do jfunct = nsd+npd+1, ifunct, 6
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij  = ad(ifunct, nci) + ad(jfunct, ncj)
         Z2   = 2.D0 * Zij
         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         alf   = ad(ifunct, nci) * ad(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         t0    = ss / Z2
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         ! Loops over nuclei, common for all shells
         do iatom = 1, my_natom
            uf = Zij * ((Q(1) - r(iatom,1)) * (Q(1) - r(iatom,1)) + &
                        (Q(2) - r(iatom,2)) * (Q(2) - r(iatom,2)) + &
                        (Q(3) - r(iatom,3)) * (Q(3) - r(iatom,3)))

            temp = 2.D0 * sqrt(Zij/pi) * ss
            s0s(iatom) = temp * FUNCT(0,uf)
            s1s(iatom) = temp * FUNCT(1,uf)
            s2s(iatom) = temp * FUNCT(2,uf)
            s3s(iatom) = temp * FUNCT(3,uf)
            s4s(iatom) = temp * FUNCT(4,uf)
         enddo

         do l1 = 1, 3
            t1   = Q(l1) - r(Nucd(ifunct),l1)
            pis  = ss  * t1
            piks = sks * t1 + alf2 * pis

            do l2 = 1, l1
               t1   = Q(l2) - r(Nucd(ifunct),l2)
               pjs  = ss  * t1
               pjks = sks * t1 + alf2 * pjs

               f1    = 1.D0
               dijs  = t1 * pis
               dijks = t1 * piks
               if (l1 .eq. l2) then
                  f1    = sq3
                  dijs  = dijs  + t0
                  dijks = dijks + sks / Z2 - alf2 * ss / (2.D0*ad(ifunct,nci))
               endif
               dijks = dijks + alf2 * dijs

               do l3 = 1, 3
                  t2 = Q(l3) - r(Nucd(jfunct),l3)

                  pipk   = t2 * pis
                  pjpk   = t2 * pjs
                  pikpk  = t2 * piks
                  pjkpk  = t2 * pjks
                  dijpk  = t2 * dijs
                  dijkpk = t2 * dijks
                  if (l1 .eq. l3) then
                     pipk   = pipk   + t0
                     dijpk  = dijpk  + pjs  / Z2
                     pikpk  = pikpk  + sks  / Z2
                     dijkpk = dijkpk + pjks / Z2
                  endif
                  if (l2 .eq. l3) then
                     pjpk   = pjpk   + t0
                     dijpk  = dijpk  + pis  / Z2
                     pjkpk  = pjkpk  + sks  / Z2
                     dijkpk = dijkpk + piks / Z2
                  endif
                  pikpk  = pikpk  + alf2 * pipk
                  pjkpk  = pjkpk  + alf2 * pjpk
                  dijkpk = dijkpk + alf2 * dijpk

                  do l4 = 1, l3
                     t1    = Q(l4) - r(Nucd(jfunct),l4)
                     ovlap = t1 * dijpk
                     tn    = t1 * dijkpk
                     f2    = 1.D0

                     if (l1 .eq. l4) then
                        ovlap = ovlap + pjpk  / Z2
                        tn    = tn    + pjkpk / Z2
                     endif
                     if (l2 .eq. l4) then
                        ovlap = ovlap + pipk  / Z2
                        tn    = tn    + pikpk / Z2
                     endif
                     if (l3.eq.l4) then
                        ovlap = ovlap+ dijs  / Z2
                        tn    = tn   + dijks / Z2 - &
                                alf2 * dijs / (2.D0*ad(jfunct,ncj))
                        f2    = sq3
                     endif

                     ! l12 and l34 range from 1 to 6, spanning the d shell in
                     ! the order xx, xy, yy, zx, zy, zz.
                     l12   = l1 * (l1-1)/2 + l2
                     l34   = l3 * (l3-1)/2 + l4
                     i_ind = ifunct + l12 -1
                     j_ind = jfunct + l34 -1
                     if (i_ind .ge. j_ind) then
                        k_ind = i_ind + ((M2d-j_ind)*(j_ind-1))/2
                        tn    = tn + alf2 * ovlap
                        cc_f  = ccoef / (f1 * f2)

                        Smatd(i_ind,j_ind) = Smatd(i_ind,j_ind) + ovlap * cc_f
                        Smatd(j_ind,i_ind) = Smatd(j_ind,i_ind) + ovlap * cc_f
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Avoids double-counting diagonal elements.
   do ifunct = 1, Md
      Smatd(ifunct,ifunct) = Smatd(ifunct,ifunct) / 2.0D0
   enddo

   ! In case of TSH
   if ( FSTSH ) Sovl_now = Smatd

   deallocate(s0s, s1s, s2s, s3s, s4s)
   return;
end subroutine

subroutine print_ml(E1, KinE, Exc, En, Etot, Pmat_vec, af, pro)
   use basis_data, only: Nucd, M, Md, nshelld, MM  
   
   implicit none 
   LIODBLE          :: E1, KinE, Exc, En, Etot 
   LIODBLE          :: Pmat_vec(MM)
   LIODBLE          :: af(Md), pro(Md)
   integer          :: ii, jj, kk, hh 
   
   write(*,*) E1, KinE, Exc, En, Etot, M, Md 
   write(*,*) nshelld(0), nshelld(1), nshelld(2), 0, 0
   do ii = 1, Md
      write(*,'(I6.1, 1x)', advance = "no") Nucd(ii)
   enddo
   write(*,*)  
   do jj = 1, MM 
      write(*,'(ES25.15E3, 1x)', advance = "no") Pmat_vec(jj)
   enddo
   write(*,*)
   do kk = 1, Md 
      write(*,'(ES25.15E3, 1x)', advance = "no") af(kk)
   enddo
   write(*,*)
   do hh = 1, Md 
      write(*,'(ES25.15E3, 1x)', advance = "no") pro(hh)
   enddo
   write(*,*)
end subroutine

end module ml_sub
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
