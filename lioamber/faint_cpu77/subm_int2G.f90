! Integrals subroutine -Second part
! Gradients only
! 2 e integrals, 2 index : density fitting functions
! All of them are calculated
! using the Obara-Saika recursive method.
module subm_int2G
contains
subroutine int2G(f)

   use liotemp   , only: FUNCT
   use garcha_mod, only: natom, M, Md, NORM, pi5, r, af, nshelld, ad, &
                         nucd, ncontd, d, cd

   implicit none
   ! aux . things
   double precision  :: Q(3), f(natom,3), ll(3)
   double precision  :: ti, tj, t0, t1, t2, t10, t11, t12, t12b
   double precision  :: t13, t13a, t13b, t14, t14a, t14b, t15, t15a, t15b
   double precision  :: t16, t16b, t17, t17b, t18, t18b, t20, t21, t22
   double precision  :: t23, t24, t25, t26, t27, t27a, t30, t31, t32, t33
   double precision  :: t40, t41, t42, t43

   double precision  :: Z2, Zij, Z2a, zc, zc2
   double precision  :: f1, f2, fs, fp, fd
   double precision  :: cc, cc1, cc2, cci, ccj
   double precision  :: ccoef, factor, roz, roZ2, sq3, alf, u

   double precision  :: sp, spj, sp1j, s1pk, s2pk
   double precision  :: s0s, s1s, s2s, s3s, s4s, s5s
   double precision  :: ps, pp, pd, pis, pip, pid, pjs, pjp, pjd
   double precision  :: pi0s, pi0p, pi0d, pj0s, pj0p, pj0d
   double precision  :: p1p, pi2s, pi2p, pi3s, pi4s, pj2s, pj2p, pj3s
   double precision  :: ds, dp, dpl, dsd
   double precision  :: d0p, d0pl, d1d, d1p, d1s, d2s, d3s, dd, df

   integer :: i, ii, j, jj, ni, nj, nsd, npd, ndd
   integer :: l, lk, lij, l1, l2, l3, l4, l5, l12, l34

   SQ3 = 1.0D0
   if (NORM) SQ3 = sqrt(3.D0)

   do l1 = 1, 3
      Ll(l1) = l1 * (l1 -1) / 2
   enddo
   nsd = nshelld(0); npd = nshelld(1); ndd = nshelld(2)

   ! (s|s)
   do ifunct = 1, nsd
   do jfunct = 1, ifunct
      f1 = 2.0D0
      if (ifunct .eq. jfunct) f1 = 1.0D0

      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         t2  = pi5 / (sqrt(Zij) * t0)

         ti   = ad(ifunct,nci) / Zij
         tj   = ad(jfunct,ncj) / Zij
         Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
         Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
         Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

         uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
         s1s = t2 * FUNCT(1,uf)

         ccoef = - cd(ifunct,nci) * cd(jfunct,ncj) * f1 * af(ifunct) *af(jfunct)
         do l1 = 1, 3
            t1 = Q(l1) - r(Nucd(ifunct),l1)
            t2 = Q(l1) - r(Nucd(jfunct),l1)
            ps = t1 * s1s
            sp = t2 * s1s
            f(Nucd(ifunct),l1) = f(Nucd(ifunct),l1) + ad(ifunct,nci) * ccoef *ps
            f(Nucd(jfunct),l1) = f(Nucd(jfunct),l1) + ad(jfunct,ncj) * ccoef *sp
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|s)
   do ifunct = nsd+1, nsd+npd, 3
   do jfunct = 1    , nsd
      do nci = 1, ncontd(ifunct)
      do ncj = 1, ncontd(jfunct)
         Zij = ad(ifunct,nci) + ad(jfunct,ncj)
         t0  = ad(ifunct,nci) * ad(jfunct,ncj)
         t2  = pi5 / (sqrt(Zij) * t0)

         ti   = ad(ifunct,nci) / Zij
         tj   = ad(jfunct,ncj) / Zij
         Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
         Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
         Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

         uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
         s0s = t2 * FUNCT(0,uf)
         s1s = t2 * FUNCT(1,uf)
         s2s = t2 * FUNCT(2,uf)
         t10 = 0.5D0 * s1s / Zij
         t20 = 0.5D0 *(s0s - tj * s1s) / ad(ifunct,nci)

         ccoef  = cd(ifunct,nci) * cd(jfunct,ncj)
         do l1 = 1, 3
            ps     = (Q(l1) - r(Nucd(ifunct),l1)) * s2s
            ccoef2 = - ccoef * af(jfunct) * af(ifunct + l1-1)

            do l2 = 1, 3
               pp = (Q(l2) - r(Nucd(jfunct),l2)) * ps
               ds = (Q(l2) - r(Nucd(ifunct),l2)) * ps
               if (l1 .eq. l2) then
                  pp = pp  + t10
                  ds = ds  + t20
                  f(Nucd(ifunct),l2) = f(Nucd(ifunct),l2) - ccoef2 * s0s
               endif

               f(Nucd(jfunct),l2) = f(Nucd(jfunct),l2) + 2.0D0 * ccoef2 * pp * &
                                    ad(jfunct,ncj)
               f(Nucd(ifunct),l2) = f(Nucd(ifunct),l2) + 2.0D0 * ccoef2 * ds * &
                                    ad(ifunct,nci)
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! (p|p) case gradients
   do ifunct = nsd+1, nsd+npd, 3
   do jfunct = nsd+1, ifunct , 3
      do nci = 1, ncontd(ifunct)
         do ncj = 1, ncontd(jfunct)
            Zij = ad(ifunct,nci) + ad(jfunct,ncj)
            Z2  = 2.0D0 * Zij
            t0  = ad(ifunct,nci) * ad(jfunct,ncj)
            t2  = pi5 / (sqrt(Zij) * t0)

            ti   = ad(ifunct,nci) / Zij
            tj   = ad(jfunct,ncj) / Zij
            Q(1) = ti * r(Nucd(ifunct),1) + tj * r(Nucd(jfunct),1)
            Q(2) = ti * r(Nucd(ifunct),2) + tj * r(Nucd(jfunct),2)
            Q(3) = ti * r(Nucd(ifunct),3) + tj * r(Nucd(jfunct),3)

            uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
            s1s = t2 * FUNCT(1,uf)
            s2s = t2 * FUNCT(2,uf)
            s3s = t2 * FUNCT(3,uf)
            t25 = s2s  / Z2

            ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
            do l1 = 1, 3
               t1   = Q(l1) - r(Nucd(ifunct),l1)
               pi0s = t1 * s1s
               pis  = t1 * s2s
               pi2s = t1 * s3s
               t15  = 0.5D0 * (pi0s - ti * pis) / ad(jfunct,ncj)
               t17  = pis / Z2
               t20  = pis / Z2

               lij = 3
               if (ifunct .eq. jfunct) lij = l1
               do l2 = 1, lij
                  t2   = Q(l2)-r(Nucd(jfunct),l2)
                  spj  = t2 * s1s
                  sp1j = t2 * s2s
                  p1p  = t2 * pi2s
                  t20  = 0.5D0 * (spj - tj * sp1j) / ad(ifunct,nci)
                  t21  = sp1j / Z2
                  if (l1 .eq. l2) p1p = p1p + t25

                  af_i = ifunct + l1 -1
                  af_j = jfunct + l2 -1
                  f1   = 0.5D0
                  if (af_i .ne. af_j) f1 = 1.0D0
                  ccoef2 = - af(af_i) * af(af_j) * f1 * ccoef

                  do l3 = 1, 3
                     pd = (Q(l3) - r(Nucd(jfunct),l3)) * p1p
                     dp = (Q(l3) - r(Nucd(ifunct),l3)) * p1p

                     if (l1 .eq. l3) then
                        dp = dp + t20
                        pd = pd + t21
                        f(Nucd(ifunct),l3) = f(Nucd(ifunct),l3) - ccoef2 * spj
                     endif
                     if (l2.eq.l3) then
                        dp = dp + t17
                        pd = pd + t15
                        f(Nucd(jfunct),l3) = f(Nucd(jfunct),l3) - ccoef2 * pi0s
                     endif

                     f(Nucd(ifunct),l3) = f(Nucd(ifunct),l3) + 2.D0 * ccoef2 * &
                                          dp * ad(ifunct,nci)
                     f(Nucd(jfunct),l3) = f(Nucd(jfunct),l3) + 2.D0 * ccoef2 * &
                                          pd *ad(jfunct,ncj)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   enddo

   ! (d|s)  gradients
   do ifunct = nsd+npd+1,Md,6
      do jfunct = 1,nsd
         dd=d(Nucd(ifunct),Nucd(jfunct))
         do nci = 1, ncontd(ifunct)
            do ncj = 1, ncontd(jfunct)
               Zij = ad(ifunct,nci) + ad(jfunct,ncj)
               Z2  = 2.0D0 * Zij
               zc=2.D0*ad(ifunct,nci)
               roz=ad(jfunct,ncj) / Zij
               alf=roz*ad(ifunct,nci)
               t0  = ad(ifunct,nci) * ad(jfunct,ncj)

               t2  = pi5 / (sqrt(Zij) * t0)
               Q(1)=(ad(ifunct,nci)*r(Nucd(ifunct),1)+ad(jfunct,ncj)*r(Nucd(jfunct),1)) / Zij
               Q(2)=(ad(ifunct,nci)*r(Nucd(ifunct),2)+ad(jfunct,ncj)*r(Nucd(jfunct),2)) / Zij
               Q(3)=(ad(ifunct,nci)*r(Nucd(ifunct),3)+ad(jfunct,ncj)*r(Nucd(jfunct),3)) / Zij
               uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
               s1s = t2 * FUNCT(1,uf)
               s2s = t2 * FUNCT(2,uf)
               s3s = t2 * FUNCT(3,uf)
               t10=(s1s-alf*s2s/ad(ifunct,nci))/zc
               ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
               do l1 = 1, 3
                  t1 = Q(l1) - r(Nucd(ifunct),l1)
                  pi0s = t1 * s1s
                  pis = t1 * s2s
                  t12=pis / Z2
                  t22=(pi0s-roz*pis)/zc
                  pi2s = t1 * s3s
                  do l2 = 1, l1
                     t1= Q(l2) - r(Nucd(ifunct),l2)
                     pj0s = t1 * s1s
                     pjs = t1 * s2s
                     t11=pjs / Z2
                     t21=(pj0s-roz*pjs)/zc
                     ds = t1 * pi2s
                     f1=1.D0
                     if (l1 .eq. l2) then
                        f1=sq3
                        ds=ds + t10
                     endif
                     l12=l1*(l1-1)/2+l2
                     ii=i+l12-1
                     cc=-af(ii)*af(j)*ccoef/f1
                     cc1=cc*2.D0
                     cci=cc1*ad(ifunct,nci)
                     ccj=cc1*ad(jfunct,ncj)

                     do l3 = 1, 3
                        t0=Q(l3)-r(Nucd(jfunct),l3)
                        t1=Q(l3)-r(Nucd(ifunct),l3)
                        dp=t0*ds
                        fs = t1 * ds
                        if (l1 .eq. l3) then
                           dp = dp + t11
                           fs=fs + t21
                           f(Nucd(ifunct),l3)=f(Nucd(ifunct),l3)-cc*pj0s
                        endif
                        if (l2.eq.l3) then
                           dp = dp + t12
                           fs=fs + t22
                           f(Nucd(ifunct),l3)=f(Nucd(ifunct),l3)-cc*pi0s
                        endif
                        f(Nucd(ifunct),l3)=f(Nucd(ifunct),l3)+cci*fs
                        f(Nucd(jfunct),l3)=f(Nucd(jfunct),l3)+ccj*dp
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   ! (d|p)  gradients
   do ifunct = nsd+npd+1,Md,6
      do jfunct = nsd+1,nsd+npd,3
         dd=d(Nucd(ifunct),Nucd(jfunct))
         do nci = 1, ncontd(ifunct)
            do ncj = 1, ncontd(jfunct)
               Zij = ad(ifunct,nci) + ad(jfunct,ncj)
               Z2  = 2.0D0 * Zij
               t0  = ad(ifunct,nci) * ad(jfunct,ncj)
               alf=t0 / Zij

               t2  = pi5 / (sqrt(Zij) * t0)
               zc=2.D0*ad(ifunct,nci)
               zc2=2.D0*ad(jfunct,ncj)
               roz=ad(jfunct,ncj) / Zij
               roZ2=ad(ifunct,nci) / Zij
               Q(1)=(ad(ifunct,nci)*r(Nucd(ifunct),1)+ad(jfunct,ncj)*r(Nucd(jfunct),1)) / Zij
               Q(2)=(ad(ifunct,nci)*r(Nucd(ifunct),2)+ad(jfunct,ncj)*r(Nucd(jfunct),2)) / Zij
               Q(3)=(ad(ifunct,nci)*r(Nucd(ifunct),3)+ad(jfunct,ncj)*r(Nucd(jfunct),3)) / Zij
               uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
               s0s = t2 * FUNCT(0,uf)
               s1s = t2 * FUNCT(1,uf)
               s2s = t2 * FUNCT(2,uf)
               s3s = t2 * FUNCT(3,uf)
               s4s = t2 * FUNCT(4,uf)
               t13a=s1s / Z2
               t13=s2s / Z2
               t10=(s0s-roz*s1s)/zc
               t11=(s1s-roz*s2s)/zc
               t12=(s2s-roz*s3s)/zc
               ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
               do l1 = 1, 3
                  t1 = Q(l1) - r(Nucd(ifunct),l1)
                  pis = t1 * s2s
                  pi2s = t1 * s3s
                  pi3s = t1 * s4s
                  t14=pi2s / Z2
                  do l2 = 1, l1
                     t1= Q(l2) - r(Nucd(ifunct),l2)
                     pjs = t1 * s2s
                     pj2s = t1 * s3s
                     t15=pj2s / Z2
                     ds = t1 * pis
                     d1s = t1 * pi2s
                     d2s = t1 * pi3s
                     f1=1.D0
                     if (l1 .eq. l2) then
                        ds=ds + t10
                        d1s=d1s + t11
                        d2s=d2s + t12
                        f1=sq3
                     endif
                     t18=(ds-roZ2*d1s)/zc2
                     t23=d1s / Z2
                     do l3 = 1, 3
                        t0=Q(l3)-r(Nucd(jfunct),l3)
                        dp=t0*d2s
                        pip=t0*pi2s
                        pi0p=t0*pis
                        pjp=t0*pj2s
                        pj0p=t0*pjs
                        if (l1 .eq. l3) then
                           dp = dp + t15
                           pip=pip + t13
                           pi0p=pi0p + t13a
                        endif
                        if (l2.eq.l3) then
                           dp = dp + t14
                           pjp=pjp + t13
                           pj0p=pj0p + t13a
                        endif
                        t16=pip / Z2
                        t17=pjp / Z2
                        t21=(pj0p-roz*pjp)/zc
                        t22=(pi0p-roz*pip)/zc
                        l12=l1*(l1-1)/2+l2
                        ii=i+l12-1
                        jj=j+l3-1
                        cc=-af(ii)*af(jj)*ccoef/f1
                        cc1=cc*2.D0
                        cci=cc1*ad(ifunct,nci)
                        ccj=cc1*ad(jfunct,ncj)

                        do l4 = 1, 3
                           t0=Q(l4)-r(Nucd(jfunct),l4)
                           t1=Q(l4)-r(Nucd(ifunct),l4)
                           dsd=t0*dp
                           fp = t1 * dp
                           if (l1 .eq. l4) then
                              dsd=dsd + t17
                              fp=fp + t21
                              f(Nucd(ifunct),l4)=f(Nucd(ifunct),l4)-cc*pj0p
                           endif
                           if (l2.eq.l4) then
                              dsd=dsd + t16
                              fp=fp + t22
                              f(Nucd(ifunct),l4)=f(Nucd(ifunct),l4)-cc*pi0p
                           endif
                           if (l3.eq.l4) then
                              f(Nucd(jfunct),l4)=f(Nucd(jfunct),l4)-cc*ds
                              dsd=dsd + t18
                              fp=fp + t23
                           endif
                           f(Nucd(ifunct),l4)=f(Nucd(ifunct),l4)+cci*fp
                           f(Nucd(jfunct),l4)=f(Nucd(jfunct),l4)+ccj*dsd
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   ! (d|d) gradients
   do ifunct = nsd+npd+1,Md,6
      do jfunct = nsd+npd+1,i,6
         dd=d(Nucd(ifunct),Nucd(jfunct))
         do nci = 1, ncontd(ifunct)
            do ncj = 1, ncontd(jfunct)
               Zij = ad(ifunct,nci) + ad(jfunct,ncj)
               Z2  = 2.0D0 * Zij
               t0  = ad(ifunct,nci) * ad(jfunct,ncj)
               alf=t0 / Zij

               t2  = pi5 / (sqrt(Zij) * t0)
               zc=2.D0*ad(ifunct,nci)
               zc2=2.D0*ad(jfunct,ncj)
               roz=ad(jfunct,ncj) / Zij
               roZ2=ad(ifunct,nci) / Zij
               Q(1)=(ad(ifunct,nci)*r(Nucd(ifunct),1)+ad(jfunct,ncj)*r(Nucd(jfunct),1)) / Zij
               Q(2)=(ad(ifunct,nci)*r(Nucd(ifunct),2)+ad(jfunct,ncj)*r(Nucd(jfunct),2)) / Zij
               Q(3)=(ad(ifunct,nci)*r(Nucd(ifunct),3)+ad(jfunct,ncj)*r(Nucd(jfunct),3)) / Zij
               uf  = d(Nucd(ifunct),Nucd(jfunct)) * t0 / Zij
               s0s = t2 * FUNCT(0,uf)
               s1s = t2 * FUNCT(1,uf)
               s2s = t2 * FUNCT(2,uf)
               s3s = t2 * FUNCT(3,uf)
               s4s = t2 * FUNCT(4,uf)
               s5s = t2 * FUNCT(5,uf)
               t13a=s1s / Z2
               t13=s2s / Z2
               t13b=s3s / Z2
               t10=(s0s-roz*s1s)/zc
               t11=(s1s-roz*s2s)/zc
               t12=(s2s-roz*s3s)/zc
               t12b=(s3s-roz*s4s)/zc
               ccoef = cd(ifunct,nci) * cd(jfunct,ncj)
               do l1 = 1, 3
                  t1 = Q(l1) - r(Nucd(ifunct),l1)
                  pi0s = t1 * s1s
                  pis = t1 * s2s
                  pi2s = t1 * s3s
                  pi3s = t1 * s4s
                  pi4s = t1 * s5s
                  t14a=pis / Z2
                  t14b=pi3s / Z2
                  t14=pi2s / Z2
                  t25=(pi0s-roZ2*pis)/zc2
                  t26=(pis-roZ2*pi2s)/zc2
                  do l2 = 1, l1
                     t1= Q(l2) - r(Nucd(ifunct),l2)
                     pj0s = t1 * s1s
                     pjs = t1 * s2s
                     pj2s = t1 * s3s
                     pj3s = t1 * s4s
                     t15a=pjs / Z2
                     t15b=pj3s / Z2
                     t15=pj2s / Z2
                     t23=(pj0s-roZ2*pjs)/zc2
                     t24=(pjs-roZ2*pj2s)/zc2
                     ds = t1 * pis
                     d1s = t1 * pi2s
                     d2s = t1 * pi3s
                     d3s = t1 * pi4s
                     f1=1.D0
                     if (l1 .eq. l2) then
                        ds=ds + t10
                        d1s=d1s + t11
                        d2s=d2s + t12
                        d3s=d3s + t12b
                        f1=sq3
                     endif
                     t18b=(d1s-roZ2*d2s)/zc2
                     lij=3
                     if (ifunct .eq. jfunct) then
                        lij=l1
                     endif
                     do l3 = 1, lij
                        t0=Q(l3)-r(Nucd(jfunct),l3)
                        s1pk=t0*s2s
                        s2pk=t0*s3s
                        t27a=s1pk / Z2
                        t27=s2pk / Z2
                        d0p=t0*d1s
                        dp=t0*d2s
                        d1p=t0*d3s
                        pi2p=t0*pi3s
                        pip=t0*pi2s
                        pi0p=t0*pis
                        pj2p=t0*pj3s
                        pjp=t0*pj2s
                        pj0p=t0*pjs
                        if (l1 .eq. l3) then
                           d0p = d0p + t15a
                           dp = dp + t15
                           d1p = d1p + t15b
                           pi2p=pi2p + t13b
                           pip=pip + t13
                           pi0p=pi0p + t13a
                        endif
                        if (l2.eq.l3) then
                           d0p = d0p + t14a
                           dp = dp + t14
                           d1p = d1p + t14b
                           pj2p=pj2p + t13b
                           pjp=pjp + t13
                           pj0p=pj0p + t13a
                        endif
                        t16b=pi2p / Z2
                        t17b=pj2p / Z2
                        t33=dp / Z2
                        t43=(d0p-roZ2*dp)/zc2
                        lk=l3
                        if (ifunct .eq. jfunct) then
                           lk=min(l3,Ll(l1)-Ll(l3)+l2)
                        endif

                        do l4 = 1, lk
                           t0=Q(l4)-r(Nucd(jfunct),l4)
                           d1d=t0*d1p
                           d0pl=t0*d1s
                           dpl=t0*d2s
                           pi0d=t0*pip
                           pid=t0*pi2p
                           pj0d=t0*pjp
                           pjd=t0*pj2p
                           if (l1 .eq. l4) then
                              d1d=d1d + t17b
                              d0pl=d0pl + t15a
                              dpl=dpl + t15
                              pi0d = pi0d + t27a
                              pid = pid + t27
                           endif
                           if (l2.eq.l4) then
                              d1d=d1d + t16b
                              d0pl=d0pl + t14a
                              dpl=dpl + t14
                              pj0d = pj0d + t27a
                              pjd = pjd + t27
                           endif
                           f2=1.D0
                           if (l3.eq.l4) then
                              d1d=d1d + t18b
                              pi0d = pi0d + t25
                              pid = pid + t26
                              pj0d = pj0d + t23
                              pjd = pjd + t24
                              f2=sq3
                           endif
                           l12=l1*(l1-1)/2+l2
                           l34=l3*(l3-1)/2+l4
                           ii=i+l12-1
                           jj=j+l34-1
                           factor=1.D0
                           if (ii.ne.jj) then
                              factor=2.D0
                           endif
                           t30=(pj0d-roz*pjd)/zc
                           t31=(pi0d-roz*pid)/zc
                           t32=dpl / Z2
                           t40=pjd / Z2
                           t41=pid / Z2
                           t42=(d0pl-roZ2*dpl)/zc2
                           cc=-0.5D0*factor*af(ii)*af(jj)*ccoef/(f1*f2)
                           cc1=cc*2.D0
                           cci=cc1*ad(ifunct,nci)
                           ccj=cc1*ad(jfunct,ncj)
                           do l5=1,3
                              t1=Q(l5)-r(Nucd(ifunct),l5)
                              t2=Q(l5)-r(Nucd(jfunct),l5)
                              df = t2 * d1d
                              fd = t1 * d1d
                              if (l1 .eq. l5) then
                                 df=df + t40
                                 fd=fd + t30
                                 f(Nucd(ifunct),l5)=f(Nucd(ifunct),l5)-cc*pj0d
                              endif
                              if (l2.eq.l5) then
                                 df=df + t41
                                 fd=fd + t31
                                 f(Nucd(ifunct),l5)=f(Nucd(ifunct),l5)-cc*pi0d
                              endif
                              if (l3.eq.l5) then
                                 df=df + t42
                                 fd=fd + t32
                                 f(Nucd(jfunct),l5)=f(Nucd(jfunct),l5)-cc*d0pl
                              endif
                              if (l4.eq.l5) then
                                 df=df + t43
                                 fd=fd + t33
                                 f(Nucd(jfunct),l5)=f(Nucd(jfunct),l5)-cc*d0p
                              endif
                              f(Nucd(ifunct),l5)=f(Nucd(ifunct),l5)+cci*fd
                              f(Nucd(jfunct),l5)=f(Nucd(jfunct),l5)+ccj*df
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   return
end subroutine
end module subm_int2G
