!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Integrals subroutines - Int2                                                 !
! 2e integrals subroutine, 2 indexes: density fitting functions. All of them   !
! are calculated using the Obara-Saika recursive method, which loops over all  !
! basis functions. Basis functions are supposed to be ordered according to the !
! type, first all s, then all p, then all d... inside each type, they are      !
! ordered in shells (px, py, pz, dxx, dxy, dyy, dzx, dzy, dzz).                !
! Input: Density basis                                                         !
! Output: G matrix, which should be inverted when evaluating Coulomb terms.    !
! Refactored in 2018 by F. Pedron                                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! SACAR XX, solo se usa aca
module subm_int2
contains
subroutine int2()
   use liotemp   , only: FUNCT
   use garcha_mod, only: RMM, ngd, M, Md, ad, nucd, ncontd, r, d, cd, &
                         NORM, pi5, nshelld

   implicit none
   double precision, allocatable :: dgelss_temp(:), inv_work(:), trabajo(:)
   integer, allocatable :: XXX(:), aux(:), XX(:,:)

   double precision :: Q(3)

   ! Ex Implicits
   double precision  :: t0, t1, t2, t3, t4, t5
   double precision  :: ss, s0s, s1s, s2s, s3s, s4s
   double precision  :: ps, pjs, pjp, pj2s, pis, pip, pi2s, pi3s
   double precision  :: alf, cc, ccoef, d1s, d2s, dd, dp, ds
   double precision  :: roz, rcond, f1, f2, sq3
   double precision  :: uf, tmp, tn, tj, ti, t6, Z2, za, zc, Zij

   integer :: igpu, info
   integer :: j_ind, i_ind, k_ind
   integer :: ifunct, jfunct, nci, ncj
   integer :: nsd, npd, ndd
   integer :: lll, l12, l34, l1, l2, l3, l4, lij, lk
   integer :: MM, MMp, MMd, Md2, Md3, Md5
   integer :: M2, M7, M9, M10, M12, M15

   allocate(inv_work(Md), XXX(8*Md), aux(ngd), XX(Md, Md))

   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)

   nsd = nshelld(0); npd = nshelld(1); ndd = nshelld(2)
   Md2 = 2*Md; M2 = 2*M; MM = M*(M+1)/2; MMd = Md*(Md+1)/2

   M7  = 1 + 3*MM ! G
   M9  = M7 + MMd !  Gm
   M15 = M9+MMd+MM+M ! Auxiliars

   do ifunct = 1, MMd
       RMM(M7+ifunct-1) = 0.D0
   enddo

   ! 2 index electron repulsion for density basis set
   ! First loop - (s|s)
   do ifunct = 1, nsd
   do jfunct = 1, ifunct
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         t1  = sqrt(Zij) * t0
         alf = t0 / Zij

         uf    = alf * dd
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
         s0s   = pi5 / t1 * FUNCT(0,uf)

         k_ind = ifunct +((Md2-jfunct)*(jfunct-1))/2
         RMM(M7+k_ind-1) = RMM(M7+k_ind-1) + ccoef * s0s
      enddo
      enddo
   enddo
   enddo

   ! Second loop - (p|s)
   do ifunct = nsd+1, nsd+npd, 3
   do jfunct = 1    , nsd
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         alf = t0 / Zij
         t1  = sqrt(Zij) * t0
         t2  = pi5 / t1

         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         uf    = alf * dd
         s1s   = t2 * FUNCT(1,uf)
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         ! l1: different p in the p shell (x,y,z, respectively)
         do l1 = 1, 3
            t1 = Q(l1) - r(Nucd(ifunct),l1)
            tn = t1 * s1s

            i_ind = ifunct + l1 -1
            k_ind = i_ind  + ((Md2-jfunct)*(jfunct-1))/2
            RMM(M7+k_ind-1) = RMM(M7+k_ind-1) + tn * ccoef
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Third loop - (p|p)
   do ifunct = nsd+1, nsd+npd, 3
   do jfunct = nsd+1, ifunct , 3
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         Z2  = 2.D0 * Zij
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         alf = t0 / Zij
         t1  = sqrt(Zij) * t0
         t2  = pi5 / t1

         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         uf    = alf * dd
         s1s   = t2 * FUNCT(1,uf)
         s2s   = t2 * FUNCT(2,uf)
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         do l1 = 1, 3
            t1 = Q(l1) - r(Nucd(ifunct),l1)
            ps = t1 * s2s

            do l2 = 1, 3
               t1 = Q(l2) - r(Nucd(jfunct),l2)
               tn = t1 * ps
               if (l1 .eq. l2) then
                  tn = tn + s1s / Z2
               endif

               i_ind = ifunct + l1 -1
               j_ind = jfunct + l2 -1
               if (i_ind .ge. j_ind) then
                  k_ind = i_ind + ((Md2-j_ind)*(j_ind-1))/2
                  RMM(M7+k_ind-1) = RMM(M7+k_ind-1) + tn * ccoef
               endif
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fourth loop - (d|s)
   do ifunct = nsd+npd+1, Md, 6
   do jfunct = 1        , nsd
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         alf = t0 / Zij
         roz = ad(jfunct,ncj) / Zij
         t1  = sqrt(Zij) * t0
         t2  = pi5 / t1

         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         uf    = alf * dd
         s0s   = t2 * FUNCT(0,uf)
         s1s   = t2 * FUNCT(1,uf)
         s2s   = t2 * FUNCT(2,uf)
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         do l1 = 1, 3
            t1 = Q(l1) - r(Nucd(ifunct),l1)
            ps = t1 * s2s

            do l2 = 1, l1
               t1 = Q(l2) - r(Nucd(ifunct),l2)
               tn = t1 * ps
               f1 = 1.0D0
               if (l1 .eq. l2) then
                  tn = tn + (s0s - roz*s1s) / (2.D0 * ad(ifunct,nci))
                  f1 = sq3
               endif

               l12   = l1 * (l1-1)/2 + l2
               i_ind = ifunct + l12 -1
               k_ind = i_ind + ((Md2-jfunct)*(jfunct-1))/2
               RMM(M7+k_ind-1) = RMM(M7+k_ind-1) + tn * ccoef / f1
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fifth loop - (d|p)
   do ifunct = nsd+npd+1, Md     , 6
   do jfunct = nsd+1    , nsd+npd, 3
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         Z2  = 2.D0 * Zij
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         alf = t0 / Zij
         t1  = sqrt(Zij) * t0
         t2  = pi5 / t1

         Q(1) = (ad(ifunct,nci) * r(Nucd(ifunct),1) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),1)) / Zij
         Q(2) = (ad(ifunct,nci) * r(Nucd(ifunct),2) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),2)) / Zij
         Q(3) = (ad(ifunct,nci) * r(Nucd(ifunct),3) + &
                 ad(jfunct,ncj) * r(Nucd(jfunct),3)) / Zij

         uf    = alf * dd
         s1s   = t2 * FUNCT(1,uf)
         s2s   = t2 * FUNCT(2,uf)
         s3s   = t2 * FUNCT(3,uf)
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)

         do l1 = 1, 3
            t1   = Q(l1) - r(Nucd(ifunct),l1)
            pis  = t1 * s2s
            pi2s = t1 * s3s

            do l2 = 1, l1
               t1  = Q(l2) - r(Nucd(ifunct),l2)
               pjs = t1 * s2s
               ds  = t1 * pi2s
               f1  = 1.D0

               if (l1 .eq. l2) then
                  f1 = sq3
                  ds = ds + (s1s - alf * s2s/ ad(ifunct,nci)) / &
                            (2.D0 * ad(ifunct,nci))
               endif
               ! Index of p
               do l3 = 1, 3
                  t0 = Q(l3) - r(Nucd(jfunct),l3)
                  tn = t0 * ds

                  if (l1 .eq. l3) then
                     tn = tn + pjs / Z2
                  endif
                  if (l2 .eq. l3) then
                     tn = tn + pis / Z2
                  endif
                  l12   = l1 * (l1-1)/2 + l2
                  i_ind = ifunct + l12 -1
                  j_ind = jfunct + l3  -1

                  k_ind = i_ind + ((Md2-j_ind)*(j_ind-1))/2
                  RMM(M7+k_ind-1) = RMM(M7+k_ind-1) + tn * ccoef / f1
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Sixth and final loop - (d|d)
   do ifunct = nsd+npd+1, Md    , 6
   do jfunct = nsd+npd+1, ifunct, 6
      dd = d(Nucd(ifunct), Nucd(jfunct))

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         Z2  = 2.D0 * Zij
         za  = 2.D0 * ad(ifunct,nci)
         zc  = 2.D0 * ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         alf = t0 / Zij
         t1  = sqrt(Zij) * t0
         t2  = pi5 / t1
         ti  = ad(ifunct,nci) / Zij
         tj  = ad(jfunct,ncj) / Zij

         Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
         Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
         Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

         uf  = alf * dd
         s0s = t2 * FUNCT(0,uf)
         s1s = t2 * FUNCT(1,uf)
         s2s = t2 * FUNCT(2,uf)
         s3s = t2 * FUNCT(3,uf)
         s4s = t2 * FUNCT(4,uf)

         t3    = (s0s - tj * s1s) / za
         t4    = (s1s - tj * s2s) / za
         t5    = (s2s - tj * s3s) / za
         ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
         do l1 = 1, 3
            t1   = Q(l1) - r(Nucd(ifunct),l1)
            pis  = t1 * s2s
            pi2s = t1 * s3s
            pi3s = t1 * s4s

            do l2 = 1, l1
               t1   = Q(l2) - r(Nucd(ifunct),l2)
               pjs  = t1 * s2s
               pj2s = t1 * s3s
               ds   = t1 * pis
               d1s  = t1 * pi2s
               d2s  = t1 * pi3s
               f1   = 1.D0

               if (l1 .eq. l2) then
                  ds  = ds  + t3
                  d1s = d1s + t4
                  d2s = d2s + t5
                  f1  = sq3
               endif
               t6 = (ds - ti * d1s) / zc

               lij = 3
               if (ifunct .eq. jfunct) then
                  lij = l1
               endif
               do l3 = 1, lij
                  t0  = Q(l3) - r(Nucd(jfunct),l3)
                  dp  = t0 * d2s
                  pip = t0 * pi2s
                  pjp = t0 * pj2s

                  if (l1 .eq. l3) then
                     dp  = dp  + pj2s / Z2
                     pip = pip + s2s  / Z2
                  endif
                  if (l2 .eq. l3) then
                     dp  = dp  + pi2s / Z2
                     pjp = pjp + s2s  / Z2
                  endif

                  if (ifunct .eq. jfunct) then
                     lk = l1 * (l1-1)/2 - l3 * (l3-1)/2 + l2
                     lk = min(l3, lk)
                  else
                     lk = l3
                  endif
                  do l4 = 1, lk
                     t0 = Q(l4) - r(Nucd(jfunct),l4)
                     tn = t0 * dp

                     if (l1 .eq. l4) then
                        tn = tn + pjp / Z2
                     endif
                     if (l2 .eq. l4) then
                        tn = tn + pip / Z2
                     endif
                     f2 = 1.D0
                     if (l3 .eq. l4) then
                        tn = tn + t6
                        f2 = sq3
                     endif

                     l12 = l1 * (l1-1)/2 + l2
                     l34 = l3 * (l3-1)/2 + l4
                     i_ind = ifunct + l12 -1
                     j_ind = jfunct + l34 -1
                     k_ind = i_ind  + (Md2-j_ind)*(j_ind-1)/2
                     RMM(M7+k_ind-1) = RMM(M7+k_ind-1) + tn * ccoef / (f1 * f2)
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   do ifunct = 1, Md
   do jfunct = 1, Md
      if (ifunct .ge. jfunct) then
         k_ind = ifunct + (Md*2-jfunct)*(jfunct-1)/2
      else
         k_ind = jfunct + (Md*2-ifunct)*(ifunct-1)/2
      endif

      XX(ifunct ,jfunct) = RMM(M7+k_ind-1)
   enddo
   enddo

   k_ind = 0
   do jfunct = 1     , Md
   do ifunct = jfunct, Md
      k_ind = k_ind +1
      RMM(M7+k_ind-1) = XX(ifunct,jfunct)
   enddo
   enddo

   MMp = Md * (Md+1) / 2
   do k_ind = 1, MMp
      RMM(M9+k_ind-1) = 0.0D0
   enddo

   M10 = M9  + MMd + MM+1
   M12 = M10 + Md
   Md3 = 3 * Md

   aux   = 0.0D0
   Md5   = 5*Md
   rcond = 1.0D-07

   call g2g_timer_sum_start('G condition')
   call dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1, RMM(M10),-1,XXX,info)
   Md5 = RMM(M10)
   allocate(dgelss_temp(Md5))
   call dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1,dgelss_temp,Md5,XXX,info)
   deallocate(dgelss_temp)
   ss = RMM(M9) / RMM(M9+Md-1)
   call g2g_timer_sum_stop('G condition')

   ! Inversion of G matrix, kept in Gm
   call g2g_timer_sum_start('G invert')

   do ifunct = 1, Md
   do jfunct = 1, Md
      if (ifunct .ge. jfunct) then
         k_ind = ifunct + (Md*2-jfunct)*(jfunct-1)/2
      else
         k_ind = jfunct + (Md*2-ifunct)*(ifunct-1)/2
      endif
      XX(ifunct,jfunct) = RMM(M7+k_ind-1)
   enddo
   enddo

   call dsytrf('U',Md,XX,Md,XXX,RMM(M10),-1,info)
   Md5 = RMM(M10)
   allocate(dgelss_temp(Md5))
   call dsytrf('U',Md,XX,Md,XXX,dgelss_temp,Md5,info)
   deallocate(dgelss_temp)
   call dsytri('U',Md,XX,Md,XXX,inv_work,info)

   do ifunct = 1, Md
   do jfunct = 1, ifunct
      k_ind = ifunct + (Md*2-jfunct)*(jfunct-1)/2
      RMM(M9+k_ind-1) = XX(jfunct,ifunct)
   enddo
   enddo

   call g2g_timer_sum_stop('G invert')
   deallocate(inv_work, XXX, aux, xx)
   return
end subroutine
end module subm_int2
