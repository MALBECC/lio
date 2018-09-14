!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_int1
contains
subroutine int1(En, Fmat, Hmat, Smat, d, r, Iz, natom, ntatom )
!------------------------------------------------------------------------------!
! Calculates 1e integrals using the Obara-Saika recursive method. (See         !
! Helgaker, "Molecular Electronic Structure Theory" (2000). pg 339)            !
!                                                                              !
!      Input : basis function information                                      !
!      Output: Energy, F matrix, and S matrix                                  !
!-------------------------------------------------------------------------------
!      INPUT AND OUTPUT VARIABLES
!-------------------------------------------------------------------------------
!        Smat ............. Overlap matrix
!        Fmat ............. Only Fock/overlap are used here.
!        Hmat ............. 1e matrix elements
!        En ............... Electron-nucleus interaction.
!        natom ............ Number of atoms
!        d ................ Interatomic distances.
!        r ................ Nuclear positions.
!        a ................ Basis exponents.
!        c ................ Basis coefficients.
!        nshell ........... Number of basis functions per shell.
!        Nuc(ifunct) ........... Nucleus corresponding to basis function i.
!        Iz(i) ............ Atomic number Z of nucleus i.
!        ncount(i) ........ Number of contractions of bais function i.
!        M ................ Number of basis functions.
!        NORM ............. Deprecated. Boolean indicating normalization.
!-------------------------------------------------------------------------------
   use garcha_mod   , only: Nuc, a, c, ncont, NORM, M, nshell
   use liotemp      , only: FUNCT
   use constants_mod, only: pi, pi32
   implicit none

   double precision, allocatable, intent(inout) :: Smat(:,:)
   double precision, intent(inout)              :: Fmat(:), Hmat(:), En

   integer,          intent(in) :: natom, ntatom, Iz(natom)
   double precision, intent(in) :: d(natom,natom), r(ntatom,3)

   integer           :: my_natom, igpu, i_ind, j_ind, k_ind, ifunct, jfunct, &
                        iatom, jatom, nci, ncj, l1, l2, l3, l4, l12, l34,    &
                        MM, ns, np, M2, M5, M11
   double precision  :: ovlap, uf, cc_f, Q(3), temp, sq3, alf, alf2, ccoef,    &
                        t0, t1, t2, f1, f2, tn, tna, Z2, Zij, ss, ps, dd, sks, &
                        p0s, p1s, p2s, p3s, pi0p, pi1p, piks, pikpk, pipk,     &
                        pis, pj0s, pj1s, pj2s, pj0p, pj1p, pjkpk, pjks, pjpk,  &
                        pjs, pks, dijs, dijpk, dijks, dijkpk, d0s, d0p, d1p,   &
                        d1s, d2s
   double precision, allocatable, dimension(:) :: s0s, s1s, s2s, s3s, s4s

   ! Allocates internal variables.
   allocate(s0s(natom),s1s(natom),s2s(natom), s3s(natom), s4s(natom))
   if (.not.allocated(Smat)) allocate(Smat(M,M))

   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)

   ns  = nshell(0); np = nshell(1)
   MM  = M*(M+1)/2
   M2  = 2*M

!-------------------------------------------------------------------------------
! Overlap ,Kinetic energy and Nuclear Attraction matrix elements evaluation
!-------------------------------------------------------------------------------
! Overlap matrix will be kept, kinetic energy and nuclear attraction matrix
! elements wont, they're stored in Fock matrix and in the Energy directly
! in order to reduce the memory requirements
!-------------------------------------------------------------------------------

   Smat = 0.0D0
   do ifunct = 1, MM
      Fmat(ifunct) = 0.D0
      Hmat(ifunct) = 0.D0
   enddo

   ! Nuclear repulsion
   En = 0.D0
   do iatom = 1, natom
   do jatom = 1, iatom-1
      En = En + Iz(iatom)*Iz(jatom) / sqrt(d(iatom,jatom))
   enddo
   enddo

   ! Doing nuclear attraction part on GPU - KE part still is done here.
   call aint_query_gpu_level(igpu)
   my_natom = natom
   if (igpu.gt.3) my_natom = 0

   ! First loop - (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         alf  = a(ifunct,nci) * a(jfunct,ncj) / Zij
         alf2 = alf * 2.D0
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         ovlap = ss
         tn    = alf * (3.D0 - alf2*dd) * ovlap

         k_ind = ifunct + ((M2-jfunct)*(jfunct-1))/2
         Fmat(k_ind) = Fmat(k_ind) + ccoef * ovlap
         Smat(ifunct,jfunct) = Smat(ifunct,jfunct) + ccoef * ovlap
         Smat(jfunct,ifunct) = Smat(jfunct,ifunct) + ccoef * ovlap

         ! Loops over nuclei and gets nuclear attraction matrix elements.
         ! tna accumulates nuc. attraction over all nuclei.
         tna = 0.D0
         do iatom = 1, my_natom
            uf = Zij * ((Q(1) - r(iatom,1)) * (Q(1) - r(iatom,1)) + &
                        (Q(2) - r(iatom,2)) * (Q(2) - r(iatom,2)) + &
                        (Q(3) - r(iatom,3)) * (Q(3) - r(iatom,3)))
            s0s(iatom) = Iz(iatom) * 2.0D0 * sqrt(Zij/pi) * ss * FUNCT(0,uf)
            tna = tna - s0s(iatom)
         enddo
         Hmat(k_ind) = Hmat(k_ind) + ccoef * (tn + tna)

      enddo
      enddo
   enddo
   enddo

   ! Second loop - (p|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1, ns
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij = a(ifunct, nci) + a(jfunct, ncj)
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct, nci) * a(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

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
            t1 = Q(l2) - r(Nuc(ifunct),l2)
            ovlap = t1 * ss
            tn    = t1 * sks + alf2 * ovlap

            i_ind = ifunct + l2 -1
            k_ind = i_ind  + ((M2-jfunct)*(jfunct-1))/2
            Fmat(k_ind) = Fmat(k_ind) + ovlap * ccoef
            Smat(i_ind ,jfunct) = Smat(i_ind ,jfunct) + ovlap * ccoef
            Smat(jfunct, i_ind) = Smat(jfunct, i_ind) + ovlap * ccoef

            ! Loops over nuclei, specific part
            tna = 0.D0
            do iatom = 1, my_natom
               tna = tna - &
                     Iz(iatom)*(t1*s0s(iatom) - (Q(l2)-r(iatom,l2))*s1s(iatom))
            enddo
            Hmat(k_ind) = Hmat(k_ind) + ccoef * (tn + tna)
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Third loop - (p|p)
   do ifunct = ns+1, ns+np , 3
   do jfunct = ns+1, ifunct, 3
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij = a(ifunct, nci) + a(jfunct, ncj)
         Z2  = 2.D0 * Zij
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct, nci) * a(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

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
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            ps  = ss * t1
            pks = sks * t1 + alf2 * ps

            do l2 = 1, 3
               t1 = Q(l2) - r(Nuc(jfunct),l2)

               ovlap = t1 * ps
               tn    = t1 * pks + alf2 * ovlap
               if (l1 .eq. l2) then
                  ovlap = ovlap + ss / Z2
                  tn    = tn + (sks + alf2*ss) / Z2
               endif

               i_ind = ifunct + l1 -1
               j_ind = jfunct + l2 -1
               if (i_ind .ge. j_ind) then
                 k_ind = i_ind + ((M2-j_ind)*(j_ind-1))/2
                 Smat(i_ind,j_ind) = Smat(i_ind,j_ind) + ovlap * ccoef
                 Smat(j_ind,i_ind) = Smat(j_ind,i_ind) + ovlap * ccoef
                 Fmat(k_ind) = Fmat(k_ind) + ovlap * ccoef
                 Hmat(k_ind) = Hmat(k_ind) + tn    * ccoef
               endif
            enddo
         enddo

         ! Loops over nuclei, specific part
         do iatom = 1, my_natom
            do l1 = 1, 3
               t1  = Q(l1) - r(Nuc(ifunct),l1)
               t2  = Q(l1) - r(iatom,l1)
               p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
               p1s = t1 * s1s(iatom) - t2 * s2s(iatom)

               do l2 = 1, 3
                  t1  = Q(l2) - r(Nuc(jfunct),l2)
                  t2  = Q(l2) - r(iatom,l2)
                  tna = t1 * p0s - t2 * p1s

                  if (l1 .eq. l2) then
                     tna = tna + (s0s(iatom) - s1s(iatom)) / Z2
                  endif

                  i_ind = ifunct + l1 -1
                  j_ind = jfunct + l2 -1
                  if (i_ind .ge. j_ind) then
                     k_ind = i_ind + ((M2-j_ind)*(j_ind-1))/2
                     Hmat(k_ind) = Hmat(k_ind) - tna * ccoef * Iz(iatom)
                  endif
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fourth loop - (d|s)
   do ifunct = ns+np+1, M, 6
   do jfunct = 1, ns
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct, nci) + a(jfunct, ncj)
         Z2   = 2.D0 * Zij
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct, nci) * a(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

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
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            ps  = ss  * t1
            pks = sks * t1 + alf2 * ps

            do l2 = 1, l1
               t1    = Q(l2) - r(Nuc(ifunct),l2)
               ovlap = t1 * ps
               tn    = t1 * pks

               f1 = 1.D0
               if (l1 .eq. l2) then
                  ovlap = ovlap + ss / Z2
                  tn    = tn    + sks/ Z2 - alf2 * ss / (2.D0*a(ifunct,nci))
                  f1    = sq3
               endif
               tn = tn + alf2 * ovlap

               ! Ordering of the d shell: xx, yx, yy, zx, zy, zz
               ! (11, 21, 22, 31, 32, 33)
               l12   = l1 * (l1-1)/2 + l2
               i_ind = ifunct + l12 -1
               k_ind = i_ind  + ((M2-jfunct)*(jfunct-1))/2

               cc_f = ccoef / f1
               Smat(i_ind,jfunct) = Smat(i_ind,jfunct)  + ovlap * cc_f
               Smat(jfunct,i_ind) = Smat(jfunct,i_ind)  + ovlap * cc_f
               Fmat(k_ind) = Fmat(k_ind) + ovlap * cc_f
               Hmat(k_ind) = Hmat(k_ind) + tn    * cc_f
            enddo
         enddo

         ! Nuclear attraction
         do iatom = 1, my_natom
         do l1 = 1, 3
            t1 = Q(l1) - r(Nuc(ifunct),l1)
            t2 = Q(l1) - r(iatom,l1)
            p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
            p1s = t1 * s1s(iatom) - t2 * s2s(iatom)

            do l2 = 1, l1
               t1  = Q(l2) - r(Nuc(ifunct),l2)
               t2  = Q(l2) - r(iatom,l2)
               tna = t1 * p0s - t2 * p1s
               f1  = 1.D0
               if (l1 .eq. l2) then
                  tna = tna + (s0s(iatom) - s1s(iatom)) / Z2
                  f1  = sq3
               endif

               ! Ordering of the d shell: xx, yx, yy, zx, zy, zz
               ! (11, 21, 22, 31, 32, 33)
               l12   = l1 * (l1-1)/2 + l2
               i_ind = ifunct + l12 -1
               k_ind = i_ind  + ((M2-jfunct)*(jfunct-1))/2
               Hmat(k_ind) = Hmat(k_ind) - &
                                     ccoef * tna * Iz(iatom) / f1
            enddo
         enddo
         enddo
         ! End of nuclear attraction

      enddo
      enddo
   enddo
   enddo

   ! Fifth loop - (d|p)
   do ifunct = ns+np+1, M    , 6
   do jfunct = ns+1   , ns+np, 3
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct, nci) + a(jfunct, ncj)
         Z2   = 2.D0 * Zij
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct, nci) * a(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

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
            t1   = Q(l1) - r(Nuc(ifunct),l1)
            pis  = ss  * t1
            piks = sks * t1 + alf2 * pis

            do l2 = 1, l1
               t1    = Q(l2) - r(Nuc(ifunct),l2)
               pjs   = ss  * t1
               pjks  = sks * t1 + alf2 * pjs
               dijs  = t1 * pis
               dijks = t1 * piks
               f1    = 1.D0

               if (l1 .eq. l2) then
                  f1    = sq3
                  dijs  = dijs  + ss / Z2
                  dijks = dijks + sks/ Z2 - alf2 * ss / (2.D0*a(ifunct,nci))
               endif
               dijks = dijks + alf2 * dijs

               do l3 = 1, 3
                  t1    = Q(l3) - r(Nuc(jfunct),l3)
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
                  k_ind = i_ind + ((M2-j_ind)*(j_ind-1))/2

                  Smat(i_ind,j_ind) = Smat(i_ind,j_ind) + cc_f * ovlap
                  Smat(j_ind,i_ind) = Smat(j_ind,i_ind) + cc_f * ovlap
                  Fmat(k_ind) = Fmat(k_ind) + cc_f * ovlap
                  Hmat(k_ind) = Hmat(k_ind) + cc_f * tn
               enddo
            enddo
         enddo

         ! Nuclear attraction
         do iatom = 1, my_natom
         do l1 = 1, 3
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            t2  = Q(l1) - r(iatom,l1)
            p0s = t1 * s0s(iatom) - t2 * s1s(iatom)
            p1s = t1 * s1s(iatom) - t2 * s2s(iatom)
            p2s = t1 * s2s(iatom) - t2 * s3s(iatom)

            do l2 = 1, l1
               t1   = Q(l2) - r(Nuc(ifunct),l2)
               t2   = Q(l2) - r(iatom,l2)
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
                  t1  = Q(l3) - r(Nuc(jfunct),l3)
                  t2  = Q(l3) - r(iatom,l3)

                  tna = t1 * d0s - t2 * d1s
                  if (l1 .eq. l3) then
                     tna = tna + (pj0s - pj1s) / Z2
                  endif
                  if (l2 .eq. l3) then
                     tna = tna + (p0s - p1s) / Z2
                  endif

                  l12   = l1 * (l1-1)/2 + l2
                  i_ind = ifunct + l12 -1
                  j_ind = jfunct + l3 -1
                  k_ind = i_ind  + ((M2-j_ind)*(j_ind-1))/2
                  Hmat(k_ind) = Hmat(k_ind) - &
                                        ccoef * tna * Iz(iatom) / f1
               enddo
            enddo
         enddo
         enddo
         ! End of nuclear attraction
      enddo
      enddo
   enddo
   enddo

   ! Sixth and final loop - (d|d)
   do ifunct = ns+np+1, M     , 6
   do jfunct = ns+np+1, ifunct, 6
      dd = d(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct, nci) + a(jfunct, ncj)
         Z2   = 2.D0 * Zij
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 a(jfunct,ncj) * r(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct, nci) * a(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         t0    = ss / Z2
         ccoef = c(ifunct,nci) * c(jfunct,ncj)

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
            t1   = Q(l1) - r(Nuc(ifunct),l1)
            pis  = ss  * t1
            piks = sks * t1 + alf2 * pis

            do l2 = 1, l1
               t1   = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = ss  * t1
               pjks = sks * t1 + alf2 * pjs

               f1    = 1.D0
               dijs  = t1 * pis
               dijks = t1 * piks
               if (l1 .eq. l2) then
                  f1    = sq3
                  dijs  = dijs  + t0
                  dijks = dijks + sks / Z2 - alf2 * ss / (2.D0*a(ifunct,nci))
               endif
               dijks = dijks + alf2 * dijs

               do l3 = 1, 3
                  t2 = Q(l3) - r(Nuc(jfunct),l3)

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
                     t1    = Q(l4) - r(Nuc(jfunct),l4)
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
                                alf2 * dijs / (2.D0*a(jfunct,ncj))
                        f2    = sq3
                     endif

                     ! l12 and l34 range from 1 to 6, spanning the d shell in
                     ! the order xx, xy, yy, zx, zy, zz.
                     l12   = l1 * (l1-1)/2 + l2
                     l34   = l3 * (l3-1)/2 + l4
                     i_ind = ifunct + l12 -1
                     j_ind = jfunct + l34 -1
                     if (i_ind .ge. j_ind) then
                        k_ind = i_ind + ((M2-j_ind)*(j_ind-1))/2
                        tn    = tn + alf2 * ovlap
                        cc_f  = ccoef / (f1 * f2)

                        Smat(i_ind,j_ind) = Smat(i_ind,j_ind) + ovlap * cc_f
                        Smat(j_ind,i_ind) = Smat(j_ind,i_ind) + ovlap * cc_f
                        Fmat(k_ind) = Fmat(k_ind) + ovlap * cc_f
                        Hmat(k_ind) = Hmat(k_ind) + tn    * cc_f
                     endif
                  enddo
               enddo
            enddo
         enddo

         ! Nuclear attraction
         do iatom = 1, my_natom
         do l1 = 1, 3
            t1 = Q(l1) - r(Nuc(ifunct),l1)
            t2 = Q(l1) - r(iatom,l1)
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

               do l3 = 1, 3
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

                  do l4 = 1, l3
                     t1  = Q(l4) - R(Nuc(jfunct),l4)
                     t2  = Q(l4) - r(iatom,l4)
                     tna = t1 * d0p - t2 * d1p
                     f2  = 1.D0

                     if (l4 .eq. l1) then
                        tna = tna + (pj0p - pj1p) / Z2
                     endif
                     if (l4 .eq. l2) then
                        tna = tna + (pi0p - pi1p) / Z2
                     endif
                     if (l4 .eq. l3) then
                        f2  = sq3
                        tna = tna + (d0s - d1s) / Z2
                     endif

                     l12   = l1 * (l1-1) / 2 + l2
                     l34   = l3 * (l3-1) / 2 + l4
                     i_ind = ifunct + l12 -1
                     j_ind = jfunct + l34 -1
                     if (i_ind .ge. j_ind) then
                        k_ind = i_ind + ((M2-j_ind)*(j_ind-1))/2
                        Hmat(k_ind) = Hmat(k_ind) - &
                                              ccoef * Iz(iatom) * tna / (f1*f2)
                     endif
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

   ! Avoids double-counting diagonal elements.
   do ifunct = 1, M
      Smat(ifunct,ifunct) = Smat(ifunct,ifunct) / 2.0D0
   enddo

   deallocate(s0s, s1s, s2s, s3s, s4s)
   return;
end subroutine
end module subm_int1
