!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% SUBM_INTDIP.F90 %%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine used for the calculation of the dipole moment in NEUTRAL          !
! (non-ionic) systems. Integrals are evaluated using the Obara Saika method.   !
! Inputs the density basis and outputs the dipole moment components.           !
! Original file: 19-1-1993                                                     !
!                                                                              !
! A loop is performed over all basis functions. Basis are supposed to be       !
! ordered according to type: first all s, then all p, then all d, etc.; and    !
! inside each type, they are ordered in shells: px, py, pz, dx2, dxy, dyy, dzx,!
! dzy, dzz, and so on.                                                         !
!                                                                              !
! ns, np, nd are markers for the end of s, p and d sections respectively.      !
! pos(Nuc(ifunct),j) is j component of position of nucleus i, jfunct = 1..3.   !
!                                                                              !
! The output is a Matrix with dipole moment integral in AO basis: uDip         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module subm_intdip
contains

subroutine intdip(uDip, pos, dists, Pmat_v)
   use basis_data   , only: a, c, Nuc, ncont, M, nshell, norm
   use constants_mod, only: pi32

   implicit none
   LIODBLE, intent(in)           :: pos(:,:)
   LIODBLE, intent(in)           :: dists(:,:)
   LIODBLE, intent(in), optional :: Pmat_v(:)
   LIODBLE, intent(out)          :: uDip(:,:)

   LIODBLE, allocatable :: temp_Pmat(:)
   LIODBLE :: aux(3), aux1(3), aux2(3), aux3(3), aux4(3), aux5(3), &
              aux6(3), srs(3), Q(3)
   LIODBLE :: sq3, alf, alf2, cc, cCoef, dd, dp, dijs, f2,  Z2, Zij, &
              ps, pis, pjs, ss, t0, t1
   integer :: M2, ns, np, nd, ii, jj, l1, l2, l3, l4, l12, l34
   integer :: ifunct, jfunct, icont, jcont, k_ind, iCrd
      
   ! Constants
   ns  = nshell(0)   
   np  = nshell(1)    
   nd  = nshell(2)    
   M2  = 2 * M    

   ! Initialization
   uDip = 0.0D0
   sq3  = 1.0D0  
   if (NORM) sq3 = 1.0D0 / sqrt( 3.0D0 )

   ! To get the "clean" dip matrix.
   allocate(temp_Pmat(M * (M +1) / 2))
   temp_Pmat = 1.0D0
   if (present(Pmat_v)) temp_Pmat = Pmat_v

   ! First Loop: <s|s> case.
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      dd = dists(Nuc(ifunct), Nuc(jfunct))

      do icont = 1, ncont(ifunct)
      do jcont = 1, ncont(jfunct)
         Zij   = a(ifunct,icont) + a(jfunct,jcont)
         alf   = a(ifunct,icont) * a(jfunct,jcont) / Zij
         ss    = pi32 * exp(-alf * dd) / (Zij * sqrt(Zij))
         k_ind = ifunct + ((M2-jfunct) * (jfunct-1)) /2
         cc    = c(ifunct,icont) * c(jfunct,jcont) * temp_Pmat(k_ind)

         do iCrd = 1, 3
            Q(iCrd)    = (a(ifunct,icont) * pos(Nuc(ifunct),iCrd)     &
                       +  a(jfunct,jcont) * pos(Nuc(jfunct),iCrd)) / Zij
            srs(iCrd)  = Q(iCrd) * ss
            uDip(iCrd,k_ind) = uDip(iCrd,k_ind) + cc * srs(iCrd)
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Second Loop: <p|s> case.
   do ifunct = ns +1, ns + np, 3
   do jfunct =     1,      ns
      dd = dists(Nuc(ifunct), Nuc(jfunct))
   
      do icont = 1, ncont(ifunct)
      do jcont = 1, ncont(jfunct)
         Zij  = a(ifunct,icont) + a(jfunct,jcont)
         alf  = a(ifunct,icont) * a(jfunct,jcont) / Zij
         ss   = pi32 * exp(-alf * dd) / (Zij * sqrt(Zij))

         do iCrd = 1, 3
            Q(iCrd)   = (a(ifunct,icont) * pos(Nuc(ifunct),iCrd)     &
                      +  a(jfunct,jcont) * pos(Nuc(jfunct),iCrd)) / Zij
            srs(iCrd) = Q(iCrd)*ss
         enddo

         cCoef = c(ifunct,icont) * c(jfunct,jcont)
         Z2    = 2.D0 * Zij
         
         ! Different types of p in the p shell (x, y, z respectively)
         ! The l1-1 term takes into account the different components of 
         ! the shell.
         do l1 = 1, 3        
            k_ind   =  ifunct + ((M2 - jfunct) * (jfunct -1)) /2 + (l1-1)
            aux     = (Q(l1) - pos(Nuc(ifunct),l1)) * srs
            aux(l1) = aux(l1) + ss / Z2

            uDip(:,k_ind) = uDip(:,k_ind) + cCoef * aux * temp_Pmat(k_ind)
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Third Loop: <p|p> case.
   do ifunct = ns +1, ns + np, 3
   do jfunct = ns +1,  ifunct, 3
      dd = dists(Nuc(ifunct), Nuc(jfunct))
      
      do icont = 1, ncont(ifunct)
      do jcont = 1, ncont(jfunct)
         Zij = a(ifunct,icont) + a(jfunct,jcont)
         t0  = a(ifunct,icont) * a(jfunct,jcont)
         alf = t0 / Zij
         ss  = pi32 * exp(-alf * dd) /(Zij * sqrt(Zij))
      
         do iCrd = 1, 3
            Q(iCrd)   = (a(ifunct,icont) * pos(Nuc(ifunct), iCrd)      &
                      +  a(jfunct,jcont) * pos(Nuc(jfunct), iCrd)) / Zij
            srs(iCrd) = Q(iCrd)*ss
         enddo

         Z2    = 2.D0 * Zij
         cCoef = c(ifunct,icont) * c(jfunct,jcont)

         do l1 = 1,3
            t1  = Q(l1) - pos(Nuc(ifunct), l1)
            ps  = ss * t1

            aux     = t1 * srs
            aux(l1) = aux(l1) + ss / Z2
                
            do l2 = 1,3
               aux1 = aux * ( Q(l2) - pos(Nuc(jfunct),l2) )                
               aux1(l2) = aux1(l2) + ps / Z2
               if (l1 == l2) aux1 = aux1 + srs / Z2

               ii = ifunct + l1 -1
               jj = jfunct + l2 -1
                    
               if (jj < ii) then
                  k_ind = ii + ((M2 - jj) * (jj -1))/2

                  uDip(:,k_ind) = uDip(:,k_ind) &
                                + cCoef * aux1 * temp_Pmat(k_ind)
               endif
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fourth Loop: <d|s> case.
   do ifunct = ns + np +1, ns + np + nd, 6
   do jfunct =          1, ns
      dd = dists(Nuc(ifunct), Nuc(jfunct))

      do icont = 1, ncont(ifunct)
      do jcont = 1, ncont(jfunct)
         Zij = a(ifunct,icont) + a(jfunct,jcont)
         alf = a(ifunct,icont) * a(jfunct,jcont) / Zij
         ss  = pi32 * exp(-alf * dd) / (Zij * sqrt(Zij))

         do iCrd = 1, 3
            Q(iCrd)   = (a(ifunct,icont) * pos(Nuc(ifunct), iCrd)     &
                      +  a(jfunct,jcont) * pos(Nuc(jfunct), iCrd)) / Zij
            srs(iCrd) = Q(iCrd) * ss
         enddo

         Z2    = 2.D0 * Zij
         alf2  = 2.D0 * alf
         cCoef = c(ifunct,icont) * c(jfunct,jcont)

         do l1 = 1, 3
            t1 = Q(l1) - pos(Nuc(ifunct),l1)
            ps = ss * t1
            
            aux     = (Q(l1) - pos(Nuc(ifunct),l1)) * srs
            aux(l1) = aux(l1) + ss / Z2

            do l2 = 1, l1
               t1   = Q(l2) - pos(Nuc(ifunct),l2)
               aux1 = aux * t1
               cc   = cCoef
                
               if (l1 == l2) then
                  aux1 = aux1 + srs/ Z2
                  cc   = cc * sq3
               endif

               aux1(l2) = aux1(l2) + ps / Z2

               ! The order of the d-shell should be as follows:
               ! xx, yx, yy, zx, zy, zz (11, 21, 22, 31, 32, 33)
               l12   = l1 * (l1 -1) /2 + l2
               ii    = ifunct + l12 - 1
               k_ind = ii + ((M2 - jfunct) * (jfunct -1)) / 2
         
               uDip(:,k_ind) = uDip(:,k_ind) + cc *aux1 * temp_Pmat(k_ind)
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fifth Loop: <d|p> case.
   do ifunct = ns + np +1, ns + np + nd, 6
   do jfunct =      ns +1,      ns + np, 3
      dd = dists(Nuc(ifunct), Nuc(jfunct))

      do icont = 1,ncont(ifunct)
      do jcont = 1,ncont(jfunct)
         Zij = a(ifunct,icont) + a(jfunct,jcont)
         alf = a(ifunct,icont) * a(jfunct,jcont) / Zij
         ss  = pi32*exp(-alf * dd) / (Zij * sqrt(Zij))

         do iCrd = 1, 3
            Q(iCrd)   = (a(ifunct,icont) * pos(Nuc(ifunct), iCrd)      &
                      +  a(jfunct,jcont) * pos(Nuc(jfunct), iCrd)) / Zij
            srs(iCrd) = Q(iCrd) * ss
         enddo

         Z2    = 2.0D0 * Zij
         cCoef = c(ifunct,icont) * c(jfunct,jcont)

         do l1 = 1, 3
            t1  = Q(l1) - pos(Nuc(ifunct),l1)
            pis = ss * t1

            ! aux : <pi|r|s>
            aux    = t1 * srs
            aux(l1)= aux(l1) + ss / Z2

            do l2 = 1, l1
               t1   = Q(l2) - pos(Nuc(ifunct),l2)
               pjs  = ss * t1 
               dijs = t1 * pis

               ! aux1 : <pj|r|s>
               aux1     = t1 * srs
               aux1(l2) = aux1(l2) + ss / Z2
                
               ! aux2 : <dij|r|s>
               aux2     = t1 *aux
               aux2(l2) = aux2(l2) + pis / Z2
               cc = cCoef

               if (l1 == l2) then
                  cc   = cc * sq3
                  aux2 = aux2 + srs / Z2
                  dijs = dijs + ss  / Z2
               endif

               do l3 = 1, 3
                  t1   = Q(l3) - pos(Nuc(jfunct),l3)
                  aux3 = aux2  * t1

                  if (l1 == l3) aux3 = aux3 + aux1 / Z2
                  if (l2 == l3) aux3 = aux3 + aux  / Z2

                  aux3(l3) = aux3(l3) + dijs / Z2

                  l12   = l1 * (l1 -1) /2 + l2
                  ii    = ifunct + l12 -1
                  jj    = jfunct + l3  -1
                  k_ind = ii + ((M2 - jj) * (jj -1)) / 2

                  uDip(:,k_ind) = uDip(:,k_ind) + cc *aux3 * temp_Pmat(k_ind)
               enddo
            enddo    
         enddo        
      enddo
      enddo
   enddo
   enddo

   ! Sixth Loop: <d|d> case.
   do ifunct = ns + np +1, ns + np + nd, 6
   do jfunct = ns + np +1, ns + np + nd, 6
      dd = dists(Nuc(ifunct), Nuc(jfunct))
 
      do icont = 1, ncont(ifunct)
      do jcont = 1, ncont(jfunct)
         Zij = a(ifunct,icont) + a(jfunct,jcont)
         alf = a(ifunct,icont) * a(jfunct,jcont) / Zij
         ss  = pi32 * exp(-alf * dd) / (Zij * sqrt(Zij))

         do iCrd = 1, 3
            Q(iCrd)   = (a(ifunct,icont) * pos(Nuc(ifunct), iCrd)     &
                      +  a(jfunct,jcont) * pos(Nuc(jfunct), iCrd)) / Zij
            srs(iCrd) = Q(iCrd)*ss
         enddo

         alf2  = 2.0D0 * alf
         Z2    = 2.0D0 * Zij
         cCoef = c(ifunct,icont) * c(jfunct,jcont)

         do l1 = 1, 3
            t1  = Q(l1) - pos(Nuc(ifunct),l1)
            pis = ss * t1
            
            ! aux : <pi|r|s>
            aux     = t1 * srs
            aux(l1) = aux(l1) + ss / Z2

            do l2 = 1,l1
               t1   = Q(l2) - pos(Nuc(ifunct),l2)
               pjs  = ss * t1
               dijs = t1 * pis
                
               ! aux1 : <pj|r|s>
               aux1     = t1 * srs
               aux1(l2) = aux1(l2) + ss / Z2

               ! Aux2 : (dij|r|s)
               aux2     = t1 * aux
               aux2(l2) = aux2(l2) + pis / Z2
               cc = cCoef

               if (l1 == l2) then
                  cc   = cc * sq3
                  aux2 = aux2 + srs / Z2
                  dijs = dijs + ss  / Z2
               endif

               do l3 = 1,3
                  t1 = Q(l3) - pos(Nuc(jfunct),l3)
                  dp = t1 *dijs
                
                  ! aux3 : (dij|r|pk)
                  aux3 = aux2 * t1

                  ! aux4 : (pi|r|pk)
                  aux4 = aux  * t1

                  ! aux5 : (pj|r|pk)
                  aux5 = aux1 * t1

                  if (l1 == l3) then
                     dp   = dp + pjs / Z2
                     aux3 = aux3 + aux1 / Z2
                     aux4 = aux4 + srs  / Z2
                  endif
                  if (l2 == l3) then
                     dp   = dp + pis / Z2
                     aux3 = aux3 + aux / Z2
                     aux5 = aux5 + srs / Z2
                  endif

                  aux3(l3) = aux3(l3) + dijs / Z2
                  aux4(l3) = aux4(l3) + pis  / Z2
                  aux5(l3) = aux5(l3) + pjs  / Z2

                  do l4 = 1, l3
                     t1 = Q(l4) - pos(Nuc(jfunct),l4)

                     ! aux3 : used here for (d|r|d)
                     aux6 = aux3 * t1

                     if (l1 == l4) aux6 = aux6 + aux5 / Z2
                     if (l2 == l4) aux6 = aux6 + aux4 / Z2

                     f2 = 1.D0
                     if (l3 == l4) then
                        f2   = sq3
                        aux6 = aux6 + aux2 / Z2
                     endif

                     aux6(l4) = aux6(l4) + dp / Z2
                     l12 = l1 * (l1 -1) / 2 + l2
                     l34 = l3 * (l3 -1) / 2 + l4
                     ii  = ifunct + l12 -1
                     jj  = jfunct + l34 -1

                     if (jj < ii) then
                        k_ind = ii + ((M2 - jj) * (jj -1)) /2
                        uDip(:,k_ind) = uDip(:,k_ind) &
                                      + cc * f2 * aux6 * temp_Pmat(k_ind)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   deallocate(temp_Pmat)

end subroutine intdip
end module subm_intdip