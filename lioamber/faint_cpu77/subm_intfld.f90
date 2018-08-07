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
subroutine intfld(g,ux,uy,uz)
   use garcha_mod, only: RMM, a, c, d, r, nuc, ncont, nshell, Iz, OPEN, NORM, &
                         M, Md, natom, pi32
   implicit none
   double precision, intent(in) :: g
   double precision, intent(in) :: ux, uy, uz

   double precision  :: aux(3) , aux1(3), aux2(3), aux3(3), aux4(3)
   double precision  :: aux5(3), aux6(3), Q(3)

   integer :: MM
   integer :: M2, M3, M5
   integer :: l1, l2, l3, l4, l12, l34
   integer :: ns, np, nd, ni, nj, i, j, k, ii, jj

   double precision :: sq3, term, alf, ccoef, cc, f1, f2
   double precision :: ss, ps, pp, pis, pjs, dp, dd, dijs
   double precision :: sxs, sys, szs
   double precision :: t0, t1, ti, tj, Z2, Zij

   sq3 = 1.0D0
   if (NORM) sq3 = dsqrt(3.D0)

   MM=M*(M+1)/2
   M2=2*M
   ! first P
   ! now F alpha
   M3=1+MM
   ! now S, F beta also uses the same position after S was used
   M5=M3+MM

   ns = nshell(0); np = nshell(1); nd = nshell(2)

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
         RMM(M5-1 + fock_ind) = RMM(M5-1 + fock_ind) + term
         if (OPEN) RMM(M3-1 + fock_ind) = RMM(M3-1 + fock_ind) + term
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
            RMM(M5-1 + fock_ind) = RMM(M5-1 + fock_ind) + term
            if (OPEN) RMM(M3-1 + fock_ind) = RMM(M3-1 + fock_ind) + term
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|p)
   do ifunct = ns+1,ns+np,3
   do jfunct = ns+1,i,3
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

               if (l1.eq.l2) then
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

                  RMM(M5-1 + fock_ind) = RMM(M5-1 + fock_ind) + term
                  if (OPEN) RMM(M3-1 + fock_ind) = RMM(M3-1 + fock_ind) + term
               endif
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|s)
c
      do 500 ifunct = ns+np+1,M,6
      do 500 jfunct = 1,ns
c
      dd=d(Nuc(ifunct),Nuc(jfunct))
c
      do 500 nci = 1, ncont(ifunct)
      do 500 ncj = 1, ncont(jfunct)
c
      Zij = a(ifunct,nci) + a(jfunct,ncj)
      Z2   = 2.D0 * Zij
      ti = a(ifunct,nci) / Zij
      tj = a(jfunct,nci) / Zij
      Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
      Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
      Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
c
      alf=a(ifunct,nci)*a(jfunct,ncj)/Zij
      ss=pi32*exp(-alf*dd)/(Zij*sqrt(Zij))
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c
      ccoef = c(ifunct,nci) * c(jfunct,ncj)
c
      do 505 l1 = 1, 3
       t1 = Q(l1) - r(Nuc(ifunct),l1)
       ps=ss * t1
        aux(1) = t1 * sxs
        aux(2) = t1 * sys
        aux(3) = t1 * szs
c
        aux(l1) = aux(l1) + ss / Z2
c
      do 505 l2 = 1, l1
c
       t1=Q(l2)-r(Nuc(ifunct),l2)
       aux1(1) = aux(1) * t1
       aux1(2) = aux(2) * t1
       aux1(3) = aux(3) * t1
c
       f1=1.D0
       if (l1.eq.l2) then
        aux1(1) = aux1(1) + sxs / Z2
        aux1(2) = aux1(2) + sys / Z2
        aux1(3) = aux1(3) + szs / Z2
        f1=sq3
       endif
c
       aux1(l2) = aux1(l2)+ps / Z2
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       ii=i+l12-1
c
       k=iifunct + ((M2 - jfunct) * (jfunct -1)) / 2
       cc = ccoef/f1
c
      term=cc*(aux1(1)*ux+aux1(2)*uy+aux1(3)*uz)
      RMM(M5-1 + fock_ind) = RMM(M5-1 + fock_ind) + g*term
      if (OPEN) then
       RMM(M3-1 + fock_ind) = RMM(M3-1 + fock_ind) + g*term
      endif
c
c
 505  continue
 500  continue
c-----------------------------------------------------------------
c
c (d|p) case
c
      do 600 ifunct = ns+np+1,M,6
      do 600 jfunct = ns+1,ns+np,3
c
      dd=d(Nuc(ifunct),Nuc(jfunct))
c
      do 600 nci = 1, ncont(ifunct)
      do 600 ncj = 1, ncont(jfunct)
c
      Zij = a(ifunct,nci) + a(jfunct,ncj)
      Z2   = 2.D0 * Zij
      ti = a(ifunct,nci) / Zij
      tj = a(jfunct,nci) / Zij
      Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
      Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
      Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
c
      alf=a(ifunct,nci)*a(jfunct,ncj)/Zij
      ss=pi32*exp(-alf*dd)/(Zij*sqrt(Zij))
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c
      ccoef = c(ifunct,nci) * c(jfunct,ncj)
c
      do 605 l1 = 1, 3
       t1 = Q(l1) - r(Nuc(ifunct),l1)
       pis=ss * t1
c aux : (pi|r|s)
        aux(1) = t1 * sxs
        aux(2) = t1 * sys
        aux(3) = t1 * szs
        aux(l1) = aux(l1) + ss / Z2
c
      do 605 l2 = 1, l1
       t1=Q(l2)-r(Nuc(ifunct),l2)
       pjs=ss * t1
       dijs=t1*pis
c aux1 : (pj|r|s)
        aux1(1) = t1 * sxs
        aux1(2) = t1 * sys
        aux1(3) = t1 * szs
        aux1(l2) = aux1(l2) + ss / Z2
c
c aux2 : (dij|r|s)
        aux2(1)=t1*aux(1)
        aux2(2)=t1*aux(2)
        aux2(3)=t1*aux(3)
        aux2(l2)=aux2(l2)+pis / Z2
       f1=1.D0
c
       if (l1.eq.l2) then
        f1=sq3
        aux2(1)=aux2(1) + sxs / Z2
        aux2(2)=aux2(2) + sys / Z2
        aux2(3)=aux2(3) + szs / Z2
        dijs=dijs+ss / Z2
       endif
c
      do 605 l3=1,3
c
       t1=Q(l3)-r(Nuc(jfunct),l3)
c
       aux3(1)=aux2(1) * t1
       aux3(2)=aux2(2) * t1
       aux3(3)=aux2(3) * t1
c

       if (l1.eq.l3) then
        aux3(1)=aux3(1)+aux1(1) / Z2
        aux3(2)=aux3(2)+aux1(2) / Z2
        aux3(3)=aux3(3)+aux1(3) / Z2
       endif
c
       if (l2.eq.l3) then
        aux3(1)=aux3(1)+aux(1) / Z2
        aux3(2)=aux3(2)+aux(2) / Z2
        aux3(3)=aux3(3)+aux(3) / Z2
       endif
c
       aux3(l3)=aux3(l3)+dijs / Z2
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
c
c
      term=cc*(aux3(1)*ux+aux3(2)*uy+aux3(3)*uz)
      RMM(M5-1 + fock_ind) = RMM(M5-1 + fock_ind) + g*term
      if (OPEN) then
       RMM(M3-1 + fock_ind) = RMM(M3-1 + fock_ind) + g*term
      endif
c
c
 605  continue
 600  continue
c
c ----------------------------------------------------------------
c (d|d) case
c
      do 700 ifunct = ns+np+1,M,6
      do 700 jfunct = ns+np+1,M,6
c
      dd=d(Nuc(ifunct),Nuc(jfunct))
c
      do 700 nci = 1, ncont(ifunct)
      do 700 ncj = 1, ncont(jfunct)
c
      Zij = a(ifunct,nci) + a(jfunct,ncj)
      Z2   = 2.D0 * Zij
      ti = a(ifunct,nci) / Zij
      tj = a(jfunct,nci) / Zij
      Q(1) = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
      Q(2) = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
      Q(3) = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
c
      alf=a(ifunct,nci)*a(jfunct,ncj)/Zij
      ss=pi32*exp(-alf*dd)/(Zij*sqrt(Zij))
      sxs=Q(1)*ss
      sys=Q(2)*ss
      szs=Q(3)*ss
c
      ccoef = c(ifunct,nci) * c(jfunct,ncj)
c
      do 705 l1 = 1, 3
       t1 = Q(l1) - r(Nuc(ifunct),l1)
       pis=ss * t1
c aux : (pi|r|s)
        aux(1) = t1 * sxs
        aux(2) = t1 * sys
        aux(3) = t1 * szs
        aux(l1) = aux(l1) + ss / Z2
c
      do 705 l2 = 1, l1
       t1=Q(l2)-r(Nuc(ifunct),l2)
       pjs=ss * t1
       dijs=t1*pis
c aux1 : (pj|r|s)
        aux1(1) = t1 * sxs
        aux1(2) = t1 * sys
        aux1(3) = t1 * szs
        aux1(l2) = aux1(l2) + ss / Z2
c
c aux2 : (dij|r|s)
        aux2(1)=t1*aux(1)
        aux2(2)=t1*aux(2)
        aux2(3)=t1*aux(3)
        aux2(l2)=aux2(l2)+pis / Z2
       f1=1.D0
c
       if (l1.eq.l2) then
        f1=sq3
        aux2(1)=aux2(1) + sxs / Z2
        aux2(2)=aux2(2) + sys / Z2
        aux2(3)=aux2(3) + szs / Z2
        dijs=dijs+ss / Z2
       endif
c
      do 705 l3=1,3
c
       t1=Q(l3)-r(Nuc(jfunct),l3)
       dp=t1*dijs
c
c aux3 : (dij|r|pk)
       aux3(1)=aux2(1) * t1
       aux3(2)=aux2(2) * t1
       aux3(3)=aux2(3) * t1
c aux4 : (pi|r|pk)
       aux4(1) = aux(1) * t1
       aux4(2) = aux(2) * t1
       aux4(3) = aux(3) * t1
c aux5 : (pj|r|pk)
       aux5(1) = aux1(1) * t1
       aux5(2) = aux1(2) * t1
       aux5(3) = aux1(3) * t1
c
       if (l1.eq.l3) then
        dp=dp+pjs / Z2
        aux3(1)=aux3(1)+aux1(1) / Z2
        aux3(2)=aux3(2)+aux1(2) / Z2
        aux3(3)=aux3(3)+aux1(3) / Z2
        aux4(1)=aux4(1) + sxs / Z2
        aux4(2)=aux4(2) + sys / Z2
        aux4(3)=aux4(3) + szs / Z2
       endif
c
       if (l2.eq.l3) then
        dp=dp+pis / Z2
        aux3(1)=aux3(1)+aux(1) / Z2
        aux3(2)=aux3(2)+aux(2) / Z2
        aux3(3)=aux3(3)+aux(3) / Z2
        aux5(1)=aux5(1) + sxs / Z2
        aux5(2)=aux5(2) + sys / Z2
        aux5(3)=aux5(3) + szs / Z2
       endif
c
       aux3(l3)=aux3(l3)+dijs / Z2
       aux4(l3)=aux4(l3)+pis / Z2
       aux5(l3)=aux5(l3)+pjs / Z2
c
       do 705 l4=1,l3
       t1=Q(l4)-r(Nuc(jfunct),l4)
c aux3 : used here for (d|r|d)
       aux6(1)=aux3(1) * t1
       aux6(2)=aux3(2) * t1
       aux6(3)=aux3(3) * t1
c
       if (l1.eq.l4) then
        aux6(1)=aux6(1)+aux5(1) / Z2
        aux6(2)=aux6(2)+aux5(2) / Z2
        aux6(3)=aux6(3)+aux5(3) / Z2
       endif
c
       if (l2.eq.l4) then
        aux6(1)=aux6(1)+aux4(1) / Z2
        aux6(2)=aux6(2)+aux4(2) / Z2
        aux6(3)=aux6(3)+aux4(3) / Z2
       endif
c
       f2=1.D0
       if (l3.eq.l4) then
       f2=sq3
        aux6(1)=aux6(1)+aux2(1) / Z2
        aux6(2)=aux6(2)+aux2(2) / Z2
        aux6(3)=aux6(3)+aux2(3) / Z2
       endif
c
       aux6(l4)=aux6(l4)+dp / Z2
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       if (ii.ge.jj) then
       k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/(f1*f2)
c
c
      term=cc*(aux6(1)*ux+aux6(2)*uy+aux6(3)*uz)
      RMM(M5-1 + fock_ind) = RMM(M5-1 + fock_ind) + g*term
      if (OPEN) then
       RMM(M3-1 + fock_ind) = RMM(M3-1 + fock_ind) + g*term
      endif
c
       endif
 705  continue
 700  continue
c--------------------------------------------------------------
c--------------------------------------------------------------
c
      return
      end subroutine
      end module subm_intfld
