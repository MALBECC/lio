!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates the gradients of overlap using the Obara-Saika recursive method.  !
! This is used in forces calculation.                                          !
!                                                                              !
! Input : basis function information, energy weighted density matrix.          !
! Output: gradients from overlap                                               !
!                                                                              !
! Debugged (or supposed to) 29-7-92 by Dario Estrin.                           !
! Refactored by Federico Pedron 08/2018                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_intSG
contains
subroutine intSG(ff)
   use garcha_mod, only: RMM, ll, a, c, d, r, nuc, ncont, nshell, pi32, natom, &
                         M, Md, NORM
   implicit none
   double precision, intent(inout) :: ff(natom,3)
   double precision                :: Q(3)

   integer :: i, j, k, ii, jj, ni, nj
   integer :: l, lk, lij, l1, l2, l3, l4, l5, l12, l34
   integer :: ns, np, nd
   integer :: M15

   double precision  :: ovlap, fsp, sq3, alf, cc, ccoef
   double precision  :: Zij, Z2, fs, fd, f1, f2
   double precision  :: ti, tj, tx, ty, te, t0, t1, t2, t4, t5
   double precision  :: t10, t11, t12, t13, t14, t15, t16, t17
   double precision  :: ss, spi, spj, spk
   double precision  :: ps, pp, pd, pidkl, pipk, pis, pjdkl, pjpk, pjs
   double precision  :: ds, dp, dd, df, dsd, dijpk, dijpl, dijs

   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)

   ns = nshell(0); np = nshell(1); nd = nshell(2)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) / 2
   enddo

   M2 = 2 * M
   ! Energy weighted density matrix
   M15 = 1 + 2 * (M * (M +1)) + Md * (Md +1) + M

   ! (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))

         en_wgt_ind = ifunct + ((M2 - jfunct) * (jfunct -1)) / 2
         te = RMM(M15-1 + en_wgt_ind) * ccoef * 2.0D0 * ss
         t4 = te * a(ifunct,nci)
         t5 = te * a(jfunct,ncj)

         do l1 = 1, 3
            t1 = Q(l1)             - r(Nuc(ifunct),l1)
            tx = r(Nuc(ifunct),l1) - r(Nuc(jfunct),l1)

            ff(Nuc(ifunct),l1) = ff(Nuc(ifunct),l1) + t4 * t1
            ff(Nuc(jfunct),l1) = ff(Nuc(jfunct),l1) + t5 * (t1 + tx)
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1   , ns
      dd=d(Nuc(ifunct),Nuc(jfunct))
      do 300 nci = 1, ncont(ifunct)
      do 300 ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            en_wgt_ind = ifunct + l1 -1 + ((M2 - jfunct) * (jfunct -1)) / 2
            t1 = (Q(l1) - r(Nuc(ifunct),l1))
            te = RMM(M15-1 + en_wgt_ind) * ccoef  * ss
            t4 = 2.D0 * te * a(ifunct,nci)
            t5 = 2.D0 * te * a(jfunct,ncj)

            do l2 = 1, 3
               ds = (Q(l2) - r(Nuc(ifunct),l2)) * t1
               pp = (Q(l2) - r(Nuc(jfunct),l2)) * t1
               if (l1 .eq. l2) then
                  ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) - te
                  ds = ds + 1 / Z2
                  pp = pp + 1 / Z2
               endif
               ff(Nuc(ifunct),l2) = ff(Nuc(ifunct),l2) + t4 * ds
               ff(Nuc(jfunct),l2) = ff(Nuc(jfunct),l2) + t5 * pp
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|p)
   do ifunct = ns+1 , ns+np , 3
   do jfunct = ns+1 , ifunct, 3
      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            pis = ss * (Q(l1) - r(Nuc(ifunct),l1))
            t11 = pis / Z2

            lij = 3
            if (ifunct .eq. jfunct) then
               lij = l1
            endif
            do l2 = 1, lij
               t2  = Q(l2) - r(Nuc(jfunct),l2)
               spj = t2  * ss
               pp  = t2  * pis
               t13 = spj / Z2

               if (l1 .eq. l2) then
                  pp = pp + t10
               endif

               en_wgt_ind = ifunct + l1-1 + ((M2 - (jfunct + l2-1)) * &
                            (jfunct + l2-2)) / 2
               te = RMM(M15-1 + en_wgt_ind) * ccoef
               t5 = te * 2.D0 * a(jfunct,ncj)
               t4 = te * 2.D0 * a(ifunct,nci)

               do l3 = 1, 3
                  dp = (Q(l3) - r(Nuc(ifunct),l3)) * pp
                  if (l1 .eq. l3) then
                     dp = dp + t13
                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * spj
                  endif
                  if (l2 .eq. l3) then
                     dp=dp+t11
                     ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) - te * pis
                  endif
                  pd = dp + (r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)) * pp
                  ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * dp
                  ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * pd
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (d|s)
   do ifunct = ns+np+1, M, 6
   do jfunct = 1      , ns
      dd=d(Nuc(ifunct),Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         ccoef = c(ifunct,nci) * c(jfunct,ncj)
         Zij   = a(ifunct,nci) + a(jfunct,ncj)
         Z2    = 2.0D0 * Zij
         ti    = a(ifunct,nci) / Zij
         tj    = a(jfunct,ncj) / Zij

         Q(1)  = ti * r(Nuc(ifunct),1) + tj * r(Nuc(jfunct),1)
         Q(2)  = ti * r(Nuc(ifunct),2) + tj * r(Nuc(jfunct),2)
         Q(3)  = ti * r(Nuc(ifunct),3) + tj * r(Nuc(jfunct),3)
         rexp  = d(Nuc(ifunct),Nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) /Zij
         ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
         t10   = ss / Z2

         do l1 = 1, 3
            pis = ss * (Q(l1) - r(Nuc(ifunct),l1))
            t12 = pis / Z2

            do l2 = 1, l1
               t1 = Q(l2) - r(Nuc(ifunct),l2)
               pjs  = ss * t1
               dijs = t1 * pis
               t11  = pjs / Z2
               f1   = 1.D0

               if (l1 .eq. l2) then
                  f1 = sq3
                  dijs = dijs + t10
               endif

               l12 = l1 * (l1-1) / 2 + l2
               en_wgt_ind = ifunct + l12-1 + ((M2 - jfunct) * (jfunct-1)) / 2

               te = RMM(M15-1 + en_wgt_ind) * ccoef / f1
               t4 = te * 2.D0 * a(ifunct,nci)
               t5 = te * 2.D0 * a(jfunct,ncj)
               do l3 = 1, 3
                  tx = r(Nuc(ifunct),l3) - r(Nuc(jfunct),l3)
                  dp = (Q(l3) - r(Nuc(jfunct),l3)) * dijs

                  if (l1 .eq. l3) then
                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * pjs
                     dp = dp + t11
                  endif
                  if (l2 .eq. l3) then
                     ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) - te * pis
                     dp = dp + t12
                  endif
                  fs = dp - tx * dijs
                  ff(Nuc(ifunct),l3) = ff(Nuc(ifunct),l3) + t4 * fs
                  ff(Nuc(jfunct),l3) = ff(Nuc(jfunct),l3) + t5 * dp
               enddo
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo
   
c (d|p)  gradients
c
      do 600 i=ns+np+1,M,6
      do 600 j=ns+1,ns+np,3
c
      dd=d(Nuc(ifunct),Nuc(jfunct))
c
      do 600 nci = 1, ncont(ifunct)
      do 600 ncj = 1, ncont(jfunct)
c
      Zij = a(ifunct,nci) + a(jfunct,ncj)
      Z2 = 2.0D0 * Zij
      Q(1)=(a(ifunct,nci)*r(Nuc(ifunct),1)+a(jfunct,ncj)*r(Nuc(jfunct),1)) / Zij
      Q(2)=(a(ifunct,nci)*r(Nuc(ifunct),2)+a(jfunct,ncj)*r(Nuc(jfunct),2)) / Zij
      Q(3)=(a(ifunct,nci)*r(Nuc(ifunct),3)+a(jfunct,ncj)*r(Nuc(jfunct),3)) / Zij
      alf = a(ifunct,nci)*a(jfunct,ncj) / Zij
      ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
c
      ccoef = c(ifunct,nci) * c(jfunct,ncj)
c
      t0=ss / Z2
c
      do 605 l1 = 1, 3
c
       t1=Q(l1)-r(Nuc(ifunct),l1)
       pis=ss*t1
       t11=pis / Z2
      do 605 l2 = 1, l1
       t1 = Q(l2) - r(Nuc(ifunct),l2)
       pjs=ss*t1
       t12=pjs / Z2
c
       f1=1.D0
       dijs=t1*pis
c
       if (l1 .eq. l2) then
        f1=sq3
        dijs=dijs+t0
       endif
c
        t13=dijs / Z2
       do 605 l3 = 1, 3
c
       t2=Q(l3)-r(Nuc(jfunct),l3)
       pipk=t2*pis
       pjpk=t2*pjs
       dijpk=t2*dijs
c
       if (l1 .eq. l3) then
        pipk=pipk+t0
        dijpk=dijpk+t12
       endif
c
       if (l2 .eq. l3) then
        pjpk=pjpk+t0
        dijpk=dijpk+t11
       endif
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
       k=ii+((M2-jj)*(jj-1))/2
c
        cc=ccoef/f1
        te=RMM(M15-1 + en_wgt_ind)*cc
        ty=te*2.D0
        t5=ty*a(jfunct,ncj)
        t4=ty*a(ifunct,nci)

        t14=pipk / Z2
        t15=pjpk / Z2
c gradients
       do 605 l4=1,3
c
       t1=Q(l4)-r(Nuc(jfunct),l4)
       tx=r(Nuc(ifunct),l4)-r(Nuc(jfunct),l4)
       dsd=t1*dijpk
c
       if (l1 .eq. l4) then
        dsd=dsd+t15
        ff(Nuc(ifunct),l4)=ff(Nuc(ifunct),l4)-te*pjpk
       endif
c
       if (l2 .eq. l4) then
        dsd=dsd+t14
       ff(Nuc(ifunct),l4)=ff(Nuc(ifunct),l4)-te*pipk
       endif
c
       if (l3.eq.l4) then
        dsd=dsd+t13
       ff(Nuc(jfunct),l4)=ff(Nuc(jfunct),l4)-te*dijs
       endif
c
       fsp=dsd-tx*dijpk
c
        ff(Nuc(ifunct),l4)=ff(Nuc(ifunct),l4)+t4*fsp
        ff(Nuc(jfunct),l4)=ff(Nuc(jfunct),l4)+t5*dsd
c
 605  continue
c
 600  continue
c-------------------------------------------------------
c (d|d)  gradients
c
      do 700 i=ns+np+1,M,6
      do 700 j=ns+np+1,i,6
c
      dd=d(Nuc(ifunct),Nuc(jfunct))
c
      do 700 nci = 1, ncont(ifunct)
      do 700 ncj = 1, ncont(jfunct)
c
      Zij = a(ifunct,nci) + a(jfunct,ncj)
      Z2 = 2.0D0 * Zij
      Q(1)=(a(ifunct,nci)*r(Nuc(ifunct),1)+a(jfunct,ncj)*r(Nuc(jfunct),1)) / Zij
      Q(2)=(a(ifunct,nci)*r(Nuc(ifunct),2)+a(jfunct,ncj)*r(Nuc(jfunct),2)) / Zij
      Q(3)=(a(ifunct,nci)*r(Nuc(ifunct),3)+a(jfunct,ncj)*r(Nuc(jfunct),3)) / Zij
      alf = a(ifunct,nci)*a(jfunct,ncj) / Zij
      ss    = pi32 * exp(-rexp) / (Zij * sqrt(Zij))
c
      ccoef = c(ifunct,nci) * c(jfunct,ncj)
c
      t0=ss / Z2
c
      do 705 l1 = 1, 3
c
       t1=Q(l1)-r(Nuc(ifunct),l1)
       pis=ss*t1
       t17=pis / Z2
      do 705 l2 = 1, l1
       t1 = Q(l2) - r(Nuc(ifunct),l2)
       pjs=ss*t1
       t16=pjs / Z2
c
       f1=1.D0
       dijs=t1*pis
c
       if (l1 .eq. l2) then
        f1=sq3
        dijs=dijs+t0
       endif
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
       do 705 l3 = 1, lij
c
       t2=Q(l3)-r(Nuc(jfunct),l3)
       spk=t2*ss
       t15=spk / Z2
       pipk=t2*pis
       pjpk=t2*pjs
       dijpk=t2*dijs
c
       if (l1 .eq. l3) then
        pipk=pipk+t0
        dijpk=dijpk+t16
       endif
c
       if (l2 .eq. l3) then
        pjpk=pjpk+t0
        dijpk=dijpk+t17
       endif
c
      t13=dijpk / Z2
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
       do 705 l4=1,lk
c
       f2=1.D0
       t1=Q(l4)-r(Nuc(jfunct),l4)
       ovlap=t1*dijpk
       pjdkl=t1*pjpk
       pidkl=t1*pipk
       dijpl=t1*dijs
c
       if (l1 .eq. l4) then
        pidkl=pidkl+t15
        dijpl=dijpl+t16
        ovlap=ovlap+pjpk / Z2
       endif
c
       if (l2 .eq. l4) then
        pjdkl=pjdkl+t15
        dijpl=dijpl+t17
        ovlap=ovlap+pipk / Z2
       endif
c
       if (l3.eq.l4) then
        pjdkl=pjdkl+t16
        pidkl=pidkl+t17
        ovlap=ovlap+dijs / Z2
        f2=sq3
       endif
c
       t10=pjdkl / Z2
       t11=pidkl / Z2
       t12=dijpl / Z2
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
c
       cc=ccoef/(f1*f2)
        te=RMM(M15-1 + en_wgt_ind)*cc
        ty=te*2.D0
        t5=ty*a(jfunct,ncj)
        t4=ty*a(ifunct,nci)

c gradients
        do 705 l5=1,3
        t1=Q(l5)-r(Nuc(ifunct),l5)
        tx=r(Nuc(ifunct),l5)-r(Nuc(jfunct),l5)
c
        fd=t1*ovlap
c
        if (l1 .eq. l5) then
        fd=fd+t10
        ff(Nuc(ifunct),l5)=ff(Nuc(ifunct),l5)-te*pjdkl
        endif
c
        if (l2 .eq. l5) then
        fd=fd+t11
        ff(Nuc(ifunct),l5)=ff(Nuc(ifunct),l5)-te*pidkl
        endif
c
        if (l3.eq.l5) then
        fd=fd+t12
        ff(Nuc(jfunct),l5)=ff(Nuc(jfunct),l5)-te*dijpl
        endif
c
        if (l4.eq.l5) then
        fd=fd+t13
        ff(Nuc(jfunct),l5)=ff(Nuc(jfunct),l5)-te*dijpk
        endif
c
        df=fd+tx*ovlap
c
        ff(Nuc(ifunct),l5)=ff(Nuc(ifunct),l5)+t4*fd
        ff(Nuc(jfunct),l5)=ff(Nuc(jfunct),l5)+t5*df
c
 705  continue
c
 700  continue
      return
      end subroutine
      end module subm_intSG
