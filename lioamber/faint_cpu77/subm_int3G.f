c------------------------------------------------------------------
c Integrals subroutine -Third part
c 2 e integrals, 3 index : wavefunction and density fitting functions
c All of them are calculated
c using the Obara-Saika recursive method.
c
c
c loop over all basis functions
c now the basis is supposed to be ordered according to the type,
c all s, then all p, then all d, .....
c inside each type, are ordered in shells
c px,py,pz , dx2,dxy,dyy,dzx,dzy,dzz, .....
c
c ns ... marker for end of s
c np ... marker for end of p
c nd ... marker for end of d
c
c r(Nuc(i),j) j component of position of nucleus i , j=1,3
c Input : G ,F,  standard basis and density basis
c F comes, computed the 1 electron part, and here the
c Coulomb part is added, without storing the integrals
c Output: F updated with Coulomb part, also Coulomb energy
c F also updated with exchange correlation part, also energy
c is updated
c this subroutine calls the fitting for exchange-correlation
c-----------------------------------------------------------------
       module subm_int3G; contains
       subroutine int3G(f,calc_energy)

       use subm_int2G, only: int2G
       use liotemp   , only: FUNCT
       use garcha_mod, only: RMM, ll, a, b, c, d, r, nuc, nucd
     >                     , ad, af, cd, ncont, ncontd, nshell
     >                     , nshelld, natom, ng, ngd, M, Md
     >                     , NORM, pi52, rmax
c

      implicit none
c
      logical, intent(in) :: calc_energy
      real*8  :: f(natom,3)

      real*8  :: Q(3), W(3), Jx(ng), af2(ngd), ftot(3)
      real*8  :: Exc, alf, cc, ccoef, rexp, ro, roz, sq3
      real*8  :: term, u

      real*8  :: d0d, d0p, d0pk, d0pkd, d0pkp, d0pl, d0pld
      real*8  :: d0plp, d0s, d1d, d1p, d1pk, d1pkd, d1pkp
      real*8  :: d1pl, d1pld, d1plp, d1pp, d1s, d1spm
      real*8  :: d2d, d2p, dds, ddp, ddi, ddf, ddd
      real*8  :: dij2plp, dij2pkp, dfs, dfp, dp0p, dp, dijplp
      real*8  :: dfd, dd2p, dijpkp, dp0pm, dp1d, dp1p, dp1pm
      real*8  :: dp1s, dp2p, dpc, dpd, dpf, dpk, dpp, dps, ds
      real*8  :: ds0p, ds1d, ds1p, dsp, dsf, dsd, ds2pl, ds2p
      real*8  :: ds1s, f3, f2, f1, dss, dspl, fpp, fpd, fds
      real*8  :: fdp, fdd, fss, fsp, fsd, fps
      real*8  :: p0pk, p1pk, p1s, p2s, p3s, p4s, p5s, p6s
      real*8  :: pi0dd, pi0d, pds, pdp, pdd, pi0sd, pi0pp
      real*8  :: pi0p, pi0dp, pi0dkl, pi1dp, pi1dkl, pi1dd
      real*8  :: pi0spj, pi1pl, pi1pkpm, pi1pk, pi1d, pi1spl
      real*8  :: pi1sd, pi1pp, pi1plpm, pi1p, pi2pkpl, pi2pk
      real*8  :: pi2p, pi2dklp, pi2dkl, pi2spk, pi2spj, pi2pl
      real*8  :: pi2pkpm, pi3pk, pi3p, pi3dkl, pi2spl, pi2plpm
      real*8  :: pij1s, pidklp, pidkl, pi4pk, pi3pl, pip0d, pip
      real*8  :: pijs, pij3s, pij2s, pis2pk, pis1pk, pipkpl, pipk
      real*8  :: pip1d, pj0dkl, pj0dd, pj0d, pispk, pispj, pj0sd
      real*8  :: pj0s, pj0pp, pj0p, pj0dp, pjs, pjp0d, pj1p, pj1dp
      real*8  :: pj1d, pj1dkl, pj4pk, pjpk, pjp1d, pjp
      real*8  :: pj1dd, pj1plpm, pj1pl, pj1pkpm, pj1pk, pj1spl
      real*8  :: pj1sd, pj1s, pj1pp, pj2pk, pj2p, pj2dklp, pj2pl
      real*8  :: pj2pkpm, pj2pkpl, pj2dkl, pj3dkl, pj2spl, pj2s
      real*8  :: pj2plpm, pj3s, pj3pl, pj3p, pj5s, pj4s, pj3pk
      real*8  :: pjdklp, pjdkl, pjpkpl
      real*8  :: pjs1pk, pp0p, pp0d, pjs2pk, pp1p
      real*8  :: pp1d, pp0pl, pp2p, pp1s, pp1pl, ppp, ppf, ppd
      real*8  :: ps0d, ps, pps, psp, spf, psd, ps1d, dd1s
      real*8  :: dd1pn, s0pk, pss, psf, s2dpm, s2dkl, s1pkpl
      real*8  :: s1pk, s1ds, s1dpm, s2pl, s2pks, s2pkpl, s2pk
      real*8  :: s2pjpk, s4pk, s3pl, s3pks, s3pk, s3dkl, s2ds
      real*8  :: sks, sp0d, sp0js, sp1d, sp1s, sp2js, sp3js
      real*8  :: spd, spjpk, spjs, spk, spp, sps, ss0d, ss0p
      real*8  :: ss0pj, ss1d, ss1p, ss1pj, ss1pk, ss1s, ss2p
      real*8  :: ss2pj, ss2pk, ss2s, ss3s, ss4s, ss5s, ss6s
      real*8  :: ss7s, ssd, ssf, ssp, sspj, sspk, sss
      real*8  :: dd1p, dd1d, dd0pn, dd0p, dd, d5s, d4s, d3s
      real*8  :: d3pl, d3pk, d3p, d4pk, d3d, d2spm, d2s
      real*8  :: d2pl, d2pk

      real*8  :: ta, tb, ti, tj
      real*8  :: te, tee, tw, twe, tx, ty, tye, tz, tze
      real*8  :: t0, t1, t2, t2a, t2b, t3, t3a, t3b, t4, t4b
      real*8  :: t5, t5a, t5b, t5x, t5y, t6, t6a, t6b, t6c, t6d
      real*8  :: t7, t7a, t7b, t7c, t7d, t8, t8a, t8b, t9, t9b
      real*8  :: t10, t10a, t10b, t11, t11a, t11b, t12, t12a, t12b
      real*8  :: t13, t13b, t14, t14b, t15, t15a, t15b, t15p
      real*8  :: t16, t16a, t16b, t17, t17a, t17b, t18, t18a, t18b
      real*8  :: t20, t20b, t21, t21b, t22, t22a, t22b, t22c, t22p
      real*8  :: t23, t23b, t24, t24b, t25, t25b, t26, t26b
      real*8  :: t27, t27b, t28, t28b, t29, t29b, t30, t30a, t30b
      real*8  :: t31, t31b, t32, t32b, t33, t33b, t34, t34b
      real*8  :: t35, t35b, t36, t37, t38, t39, t40, t40a, t40b
      real*8  :: t41, t41b, t50, t50b, t51, t51b, t60, t60b
      real*8  :: t61, t61b, t70, t70b, t80, t80a, t80b

      real*8  :: y1, y1b, y2, y2b, y3, y3b, y4, y4b, y5, y5b
      real*8  :: y6, y6b, y7, y7b, y8, y8b, y9, y9b
      real*8  :: z2, z2a, zc, zij

      real*8  :: y10, y10b, y12, y12b, y13, y13b, y14, y14b
      real*8  :: y15, y15b, y16, y16b, y17, y17b, y18, y18b
      real*8  :: y19, y19b, y20, y21, y22, y23, y24, y25, y26, y27
      real*8  :: y28, y29, y30, y31

      integer :: igpu, MM, MMd, Md2
      integer :: M1, M2, M3, M5, M7, M9, M11, M13, M15, M17, M18
      integer :: i, ii, j, jj, k, kk, ni, nj, nk, kn, k1
      integer :: ns, nsd, nd, ndd, np, npd
      integer :: l, lk, l1, l2, l3, l4, l5, l6, l7
      integer :: lij, l12, l23, l34, l45, l56


c scratch space
c
c auxiliars
c
c------------------------------------------------------------------
c now 16 loops for all combinations, first 2 correspond to 
c wavefunction basis, the third correspond to the density fit
c Rc(k) is constructed adding t(i,j,k)*P(i,j)
c cf(k) , variationally obtained fitting coefficient, is
c obtained by adding R(i)*G-1(i,k)
c if the t(i,j,k) were not stored, then in order to evaluate
c the corresponding part of the Fock matrix, they should be
c calculated again.
c V(i,j) obtained by adding af(k) * t(i,j,k)
c
c------------------------------------------------------------------


      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md
c  pointers
c
c
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
c
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, also F later
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues ), also space used in least-squares
      M13=M11+MM
c aux ( vector for ESSl,also for energy weighted density matrix)
      M15=M13+M
c Least squares
      M17=M15+MM
c MO vectors
      M18=M17+MMd
c
c end ------------------------------------------------
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2
c
      do 2 i=1,M
 2     Jx(i)=(M2-i)*(i-1)/2
c
c      do 5 i=1,natom
c      do 5 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 5    continue
c
c     
c-----------------------------------------------------------------
c
       call int2G(f)
c
c
c------------------------------------------------------------------
c Second Part: initializations
c
c iforce 0, numerical integration for exchange-correlation terms
c iforce 1, approximate analytical expression is used
c
c      if (iforce.eq.0) then
       do k=1,Md
         B(k,2)=0.0D0
       enddo
c      endif
c
c      if (.not.OPEN) then

c        write(*,*) 'exchnum int3G'
      call g2g_timer_start('ExcG')
      call g2g_timer_sum_start('Exchange-correlation gradients')
          if (calc_energy) then
              call g2g_solve_groups(2, Exc, f)
          else
              call g2g_solve_groups(3, Exc, f)
         endif
      call g2g_timer_stop('ExcG')
      call g2g_timer_sum_stop('Exchange-correlation gradients')

      call aint_query_gpu_level(igpu)
      
      call g2g_timer_start('CoulG')
      call g2g_timer_sum_start('Coulomb gradients')
      if (igpu.gt.2) then
        call aint_coulomb_forces(f)
        
      else

c DEBUG DEBUG
c      do k=1,natom
c        write(*,'("fuerza",I,D,D,D)') k,f(k,1),f(k,2),f(k,3)
c      enddo
c
c       else
c         call g2g_timer_sum_start('ExcG')
c          if (calc_energy) then
c              call g2g_solve_groups(2, Exc, f)
c          else
c              call g2g_solve_groups(3, Exc, f)
c         endif
c      call g2g_timer_sum_stop('ExcG')
c
c        NCOa=NCO
c        NCOb=NCO+Nunp
c        call exchnum2op(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
c     >              M18,NCOa,NCOb,Exc,f)
c       endif
c
      do 215 k=1,MM
 215   RMM(M5+k-1)=RMM(M11+k-1)
c
c commented now for debuggings
      do 217 k=1,Md
 217   af2(k)=af(k)+B(k,2)
ct
c
c-------------------------------------------------------------
c (ss|s) and gradients
c
      do 310 i=1,ns
      do 310 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
      kk=i+Jx(j)
c
      do 310 ni=1,ncont(i)
      do 310 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 312
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 311 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 311 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      u=ad(k,nk)*ti*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
c
c construction of Fock matrix part
c
      term=sss*ccoef
      RMM(M5+kk-1)=RMM(M5+kk-1)+af2(k)*term
c
      t1=ccoef*RMM(kk)
      te=t1*af(k)
      ty=2.D0*te
      tw=2.D0*t1*B(k,2)
c
      do 315 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      t3=W(l1)-r(Nucd(k),l1)
      tx=r(Nuc(i),l1)-r(Nuc(j),l1)
c
      pss=t1*sss+t2*ss1s
      sps=pss+tx*sss
      ssp=t3*ss1s
c
c
c forces
       f(Nuc(i),l1)=f(Nuc(i),l1)+(ty+tw)*a(i,ni)*pss
       f(Nuc(j),l1)=f(Nuc(j),l1)+(ty+tw)*a(j,nj)*sps
       f(Nucd(k),l1)=f(Nucd(k),l1)+ty*ad(k,nk)*ssp
c
 315   continue
 311   continue
 312   continue
 310   continue
c
c-------------------------------------------------------------
c
c (ps|s) and all gradients
      do 320 i=ns+1,ns+np,3
      do 320 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 320 ni=1,ncont(i)
      do 320 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 322
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
       sks=pi52*exp(-rexp)/zij
c
      do 321 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 321 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-roz*ss1s)/z2
      tb=ss1s/z2a
c
      do 325 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
c
c Fock matrix part
      ii=i+l1-1
      kk=ii+k1
      term=ccoef*ps
      RMM(M5+kk-1)=RMM(M5+kk-1)+af2(k)*term
c
      p1s=t1*ss1s+t2*ss2s
c
      t1=ccoef*RMM(kk)
      ty=t1*af(k)
      tye=t1*B(k,2)
      tz=2.D0*ty
      tze=2.D0*tye
c
      do 325 l2=1,3
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      t3=W(l2)-r(Nucd(k),l2)
      tx=r(Nuc(i),l2)-r(Nuc(j),l2)
      dss=t1*ps+t2*p1s
      psp=t3*p1s
c
      if (l1.eq.l2) then
       dss=dss+ta
       psp=psp+tb
       f(Nuc(i),l2)=f(Nuc(i),l2)-(ty+tye)*sss
      endif
c
       pps=dss+tx*ps
c
c forces
       f(Nuc(i),l2)=f(Nuc(i),l2)+(tz+tze)*a(i,ni)*dss
       f(Nuc(j),l2)=f(Nuc(j),l2)+(tz+tze)*a(j,nj)*pps
       f(Nucd(k),l2)=f(Nucd(k),l2)+tz*ad(k,nk)*psp
c
 325   continue
 321   continue
 322   continue
 320   continue
c
c-------------------------------------------------------------
c (pp|s) and gradients
      do 330 i=ns+1,ns+np,3
      do 330 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 330 ni=1,ncont(i)
      do 330 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 332
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 331 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 331 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
c
      do 335 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      t8=p1s/z2a
      t5=(ps-roz*p1s)/z2
      p2s=t1*ss2s+t2*ss3s
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
      do 335 l2=1,lij
c
      t2=W(l2)-Q(l2)
      ta=Q(l2)-r(Nuc(j),l2)
c
      sps=ta*sss+t2*ss1s
      sp1s=ta*ss1s+t2*ss2s
      t7=sp1s/z2a
      t6=(sps-roz*sp1s)/z2
      pps=ta*ps+t2*p1s
      pp1s=ta*p1s+t2*p2s
c
      if (l1.eq.l2) then
       pps=pps+t3
       pp1s=pp1s+t4
      endif
c
c Fock matrix
c
      term=pps*ccoef
      ii=i+l1-1
      jj=j+l2-1
      kk=ii+Jx(jj)
      RMM(M5+kk-1)=RMM(M5+kk-1)+af2(k)*term
c
      t1=ccoef*RMM(kk)
      ty=t1*af(k)
      tye=t1*B(k,2)
      tw=ty*2.D0
      twe=tye*2.D0
c gradients
      do 337 l3=1,3
c
      t1=Q(l3)-r(Nuc(i),l3)
      t2=W(l3)-Q(l3)
      tx=r(Nuc(i),l3)-r(Nuc(j),l3)
      tz=W(l3)-r(Nucd(k),l3)
c
      dps=t1*pps+t2*pp1s
      ppp=tz*pp1s
c
      if (l1.eq.l3) then
       dps=dps+t6
       ppp=ppp+t7
       f(Nuc(i),l3)=f(Nuc(i),l3)-(ty+tye)*sps
      endif
c
      if (l2.eq.l3) then
       dps=dps+t5
       ppp=ppp+t8
       f(Nuc(j),l3)=f(Nuc(j),l3)-(ty+tye)*ps
      endif
c
       pds=dps+tx*pps
c
c forces
       f(Nuc(i),l3)=f(Nuc(i),l3)+(tw+twe)*a(i,ni)*dps
       f(Nuc(j),l3)=f(Nuc(j),l3)+(tw+twe)*a(j,nj)*pds
       f(Nucd(k),l3)=f(Nucd(k),l3)+tw*ad(k,nk)*ppp
 337   continue
c
 335   continue
 331   continue
 332   continue
 330   continue
c
c-------------------------------------------------------------
c
c (ds|s) and gradients
      do 340 i=ns+np+1,M,6
      do 340 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 340 ni=1,ncont(i)
      do 340 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 342
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 341 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 341 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ta=(sss-roz*ss1s)/z2
      tb=(ss1s-roz*ss2s)/z2
c
      do 345 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      t12=(ps-roz*p1s)/z2
      t13=p1s/z2a
c
      do 345 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      dss=t1*ps+t2*p1s
      ds1s=t1*p1s+t2*p2s
      pj0s=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      t10=(pj0s-roz*pj1s)/z2
      t11=pj1s/z2a
c
      f1=1.D0
      if (l1.eq.l2) then
       dss=dss+ta
       ds1s=ds1s+tb
       f1=sq3
      endif
c
      cc=ccoef/f1
      term=dss*cc
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      kk=ii+k1
      RMM(M5+kk-1)=RMM(M5+kk-1)+af2(k)*term
c gradients
c
      t1=cc*RMM(kk)
      te=t1*af(k)
      tee=t1*B(k,2)
      ty=2.D0*te
      tye=2.D0*tee
c
       do 345 l3=1,3
        t1=Q(l3)-r(Nuc(j),l3)
        t2=W(l3)-Q(l3)
        t2b=W(l3)-r(Nucd(k),l3)
        tx=r(Nuc(i),l3)-r(Nuc(j),l3)
c
       dps=t1*dss+t2*ds1s
       dsp=t2b*ds1s
c
       if (l1.eq.l3) then
        dps=dps+t10
        dsp=dsp+t11
        f(Nuc(i),l3)=f(Nuc(i),l3)-(te+tee)*pj0s
       endif
c
       if (l2.eq.l3) then
        dps=dps+t12
        dsp=dsp+t13
        f(Nuc(i),l3)=f(Nuc(i),l3)-(te+tee)*ps
       endif
c
        fss=dps-tx*dss
c
        f(Nuc(i),l3)=f(Nuc(i),l3)+a(i,ni)*(ty+tye)*fss
        f(Nuc(j),l3)=f(Nuc(j),l3)+a(j,nj)*(ty+tye)*dps
        f(Nucd(k),l3)=f(Nucd(k),l3)+ad(k,nk)*ty*dsp

 345   continue
 341   continue
 342   continue
 340   continue
c
c-------------------------------------------------------------
c (dp|s)
      do 350 i=ns+np+1,M,6
      do 350 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 350 ni=1,ncont(i)
      do 350 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 352
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 351 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 351 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t4b=(ss2s-roz*ss3s)/z2
c
      do 355 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t5=(ps-roz*p1s)/z2
      t5b=(p1s-roz*p2s)/z2
c
      do 355 l2=1,l1
c
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      t6=(pjs-roz*pj1s)/z2
      t6b=(pj1s-roz*pj2s)/z2
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
c
      f1=1.D0
      if (l1.eq.l2) then
       f1=sq3
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t4b
      endif
c
      t14=(ds-roz*d1s)/z2
      t15=d1s/z2a
c
      do 355 l3=1,3
c
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      dps=t1*ds+t2*d1s
      dp1s=t1*d1s+t2*d2s
      pi0p=t1*ps+t2*p1s
      pi1p=t1*p1s+t2*p2s
      pj0p=t1*pjs+t2*pj1s
      pj1p=t1*pj1s+t2*pj2s
c
      if (l1.eq.l3) then
       dps=dps+t6
       dp1s=dp1s+t6b
       pi0p=pi0p+t3
       pi1p=pi1p+t4
      endif
c
      if (l2.eq.l3) then
       dps=dps+t5
       dp1s=dp1s+t5b
       pj0p=pj0p+t3
       pj1p=pj1p+t4
      endif
c
      t10=(pj0p-roz*pj1p)/z2
      t11=pj1p/z2a
      t12=(pi0p-roz*pi1p)/z2
      t13=pi1p/z2a
c
      l12=Ll(l1)+l2
      ii=i+l12-1
      jj=j+l3-1
c
      cc=ccoef/f1
      term=dps*cc
c
      kk=ii+Jx(jj)
      RMM(M5+kk-1)=RMM(M5+kk-1)+af2(k)*term
c gradients
c
      t1=cc*RMM(kk)
      te=t1*af(k)
      tee=t1*B(k,2)
      ty=2.D0*te
      tye=2.D0*tee
c
       do 355 l4=1,3
        t1=Q(l4)-r(Nuc(j),l4)
        t2=W(l4)-Q(l4)
        t2b=W(l4)-r(Nucd(k),l4)
        tx=r(Nuc(i),l4)-r(Nuc(j),l4)
c
       dds=t1*dps+t2*dp1s
       dpp=t2b*dp1s
c
       if (l1.eq.l4) then
        dds=dds+t10
        dpp=dpp+t11
        f(Nuc(i),l4)=f(Nuc(i),l4)-(te+tee)*pj0p
       endif
c
       if (l2.eq.l4) then
        dds=dds+t12
        dpp=dpp+t13
        f(Nuc(i),l4)=f(Nuc(i),l4)-(te+tee)*pi0p
       endif
c
       if (l3.eq.l4) then
        dds=dds+t14
        dpp=dpp+t15
        f(Nuc(j),l4)=f(Nuc(j),l4)-(te+tee)*ds
       endif
c
        fps=dds-tx*dps
c
        f(Nuc(i),l4)=f(Nuc(i),l4)+a(i,ni)*(ty+tye)*fps
        f(Nuc(j),l4)=f(Nuc(j),l4)+a(j,nj)*(ty+tye)*dds
        f(Nucd(k),l4)=f(Nucd(k),l4)+ad(k,nk)*ty*dpp
c      
 355   continue
 351   continue
 352   continue
 350   continue
c
c-------------------------------------------------------------
c
c (dd|s)
      do 360 i=ns+np+1,M,6
      do 360 j=ns+np+1,i,6
c
      dd=d(Nuc(i),Nuc(j))
c
      do 360 ni=1,ncont(i)
      do 360 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 362
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 361 k=1,nsd
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 361 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t5b=(ss3s-roz*ss4s)/z2
c
      do 365 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      t6=(ps-roz*p1s)/z2
      t7=(p1s-roz*p2s)/z2
      t7b=(p2s-roz*p3s)/z2
c
      do 365 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      t8=(pjs-roz*pj1s)/z2
      t9=(pj1s-roz*pj2s)/z2
      t9b=(pj2s-roz*pj3s)/z2
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t5
       d3s=d3s+t5b
       f1=sq3
      endif
c
      t12=(ds-roz*d1s)/z2
      t12b=(d1s-roz*d2s)/z2
c
c test now
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 365 l3=1,lij
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      s0pk=t1*sss+t2*ss1s
      s1pk=t1*ss1s+t2*ss2s
      s2pk=t1*ss2s+t2*ss3s
      t13=(s0pk-roz*s1pk)/z2
      t13b=(s1pk-roz*s2pk)/z2
c
      pip=t1*ps+t2*p1s
      pi1p=t1*p1s+t2*p2s
      pi2p=t1*p2s+t2*p3s
      pjp=t1*pjs+t2*pj1s
      pj1p=t1*pj1s+t2*pj2s
      pj2p=t1*pj2s+t2*pj3s
      dp=t1*ds+t2*d1s
      d1p=t1*d1s+t2*d2s
      d2p=t1*d2s+t2*d3s
c
      if (l1.eq.l3) then
       pip=pip+t3
       pi1p=pi1p+t4
       pi2p=pi2p+t5
       dp=dp+t8
       d1p=d1p+t9
       d2p=d2p+t9b
      endif
c
      if (l2.eq.l3) then
       pjp=pjp+t3
       pj1p=pj1p+t4
       pj2p=pj2p+t5
       dp=dp+t6
       d1p=d1p+t7
       d2p=d2p+t7b
      endif
c
      t10=(pjp-roz*pj1p)/z2
      t10b=(pj1p-roz*pj2p)/z2
      t11=(pip-roz*pi1p)/z2
      t11b=(pi1p-roz*pi2p)/z2
      t26=(dp-roz*d1p)/z2
      t27=d1p/z2a
      
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      do 365 l4=1,lk
c
      t1=Q(l4)-r(Nuc(j),l4)
      t2=W(l4)-Q(l4)
      dds=t1*dp+t2*d1p
      dd1s=t1*d1p+t2*d2p
      pi0d=t1*pip+t2*pi1p
      pi1d=t1*pi1p+t2*pi2p
      pj0d=t1*pjp+t2*pj1p
      pj1d=t1*pj1p+t2*pj2p
      d0pl=t1*ds+t2*d1s
      d1pl=t1*d1s+t2*d2s
c
      if (l1.eq.l4) then
       dds=dds+t10
       dd1s=dd1s+t10b
       pi0d=pi0d+t13
       pi1d=pi1d+t13b
       d0pl=d0pl+t8
       d1pl=d1pl+t9
      endif
c
      if (l2.eq.l4) then
       dds=dds+t11
       dd1s=dd1s+t11b
       pj0d=pj0d+t13
       pj1d=pj1d+t13b
       d0pl=d0pl+t6
       d1pl=d1pl+t7
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       dds=dds+t12
       dd1s=dd1s+t12b
       pi0d=pi0d+t6
       pi1d=pi1d+t7
       pj0d=pj0d+t8
       pj1d=pj1d+t9
       f2=sq3
      endif
c
      t20=(pj0d-roz*pj1d)/z2
      t21=pj1d/z2a
      t22=(pi0d-roz*pi1d)/z2
      t23=pi1d/z2a
      t24=(d0pl-roz*d1pl)/z2
      t25=d1pl/z2a
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
c
      ii=i+l12-1
      jj=j+l34-1
c
      cc=ccoef/(f1*f2)
      term=dds*cc
c
      kk=ii+Jx(jj)
      RMM(M5+kk-1)=RMM(M5+kk-1)+af2(k)*term
c
c gradients
c
      t1=cc*RMM(kk)
      te=t1*af(k)
      tee=t1*B(k,2)
      ty=2.D0*te
      tye=2.D0*tee
c
      do 365 l5=1,3
      t1=Q(l5)-r(Nuc(i),l5)
      t2=W(l5)-Q(l5)
      t2b=W(l5)-r(Nucd(k),l5)
      tx=r(Nuc(i),l5)-r(Nuc(j),l5)
c
      fds=t1*dds+t2*dd1s
      ddp=t2b*dd1s
c 
      if (l1.eq.l5) then
       fds=fds+t20
       ddp=ddp+t21
       f(Nuc(i),l5)=f(Nuc(i),l5)-(te+tee)*pj0d
      endif
c
      if (l2.eq.l5) then
       fds=fds+t22
       ddp=ddp+t23
       f(Nuc(i),l5)=f(Nuc(i),l5)-(te+tee)*pi0d
      endif
c
      if (l3.eq.l5) then
       fds=fds+t24
       ddp=ddp+t25
       f(Nuc(j),l5)=f(Nuc(j),l5)-(te+tee)*d0pl
      endif
c
      if (l4.eq.l5) then
       fds=fds+t26
       ddp=ddp+t27
       f(Nuc(j),l5)=f(Nuc(j),l5)-(te+tee)*dp
      endif
c
      dfs=fds+tx*dds
      f(Nuc(i),l5)=f(Nuc(i),l5)+(ty+tye)*a(i,ni)*fds
      f(Nuc(j),l5)=f(Nuc(j),l5)+(ty+tye)*a(j,nj)*dfs
      f(Nucd(k),l5)=f(Nucd(k),l5)+ty*ad(k,nk)*ddp
c
 365   continue
 361   continue
 362   continue
 360   continue
c
c-------------------------------------------------------------
c
c (ss|p)  and gradients
      do 370 i=1,ns
      do 370 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
      kn=i+Jx(j)
c
      do 370 ni=1,ncont(i)
      do 370 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 372
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 371 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 371 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=ti
c
      zc=2.D0*ad(k,nk)
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ta=(sss-roz*ss1s)/zc
      tb=ss1s/z2a
c
      do 375 l1=1,3
      t1=W(l1)-r(Nucd(k),l1)
      ssp=t1*ss1s
      ss1p=t1*ss2s
c
      kk=k+l1-1
      term=ssp*ccoef
c
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
      t1=ccoef*RMM(kn)
      te=t1*af(kk)
      ty=2.D0*te
      tw=2.D0*t1*B(kk,2)
c gradients
      do 377 l2=1,3
      t1=W(l2)-r(Nucd(k),l2)
      t2=Q(l2)-r(Nuc(i),l2)
      t3=W(l2)-Q(l2)
      tx=r(Nuc(i),l2)-r(Nuc(j),l2)
c
      ssd=t1*ss1p
      psp=t2*ssp+t3*ss1p
c
      if (l1.eq.l2) then
       ssd=ssd+ta
       psp=psp+tb
       f(Nucd(k),l2)=f(Nucd(k),l2)-te*sss
      endif
c
      spp=psp+tx*ssp
       f(Nuc(i),l2)=f(Nuc(i),l2)+a(i,ni)*(ty+tw)*psp
       f(Nuc(j),l2)=f(Nuc(j),l2)+a(j,nj)*(ty+tw)*spp
       f(Nucd(k),l2)=f(Nucd(k),l2)+ad(k,nk)*ty*ssd
c
 377   continue
 375   continue
 371   continue
 372   continue
 370   continue
c-------------------------------------------------------------
c (ps|p) and gradients
      do 380 i=ns+1,ns+np,3
      do 380 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 380 ni=1,ncont(i)
      do 380 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 382
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 381 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 381 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      zc=2.D0*ad(k,nk)
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      roz=ti
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      t3=ss2s/z2a
      t3a=ss1s/z2a
      ss3s=t2*FUNCT(3,u)
c
      do 385 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      t9=p1s/z2a
      p2s=t1*ss2s+t2*ss3s
      t5=(ps-roz*p1s)/zc
c
      do 385 l2=1,3
      t1=W(l2)-r(Nucd(k),l2)
      ss0pj=t1*ss1s
      sspj=t1*ss2s
      t10=(ss0pj-tj*sspj)/z2
      pispj=t1*p2s
      pi0spj=t1*p1s
c
      if (l1.eq.l2) then
       pi0spj=pi0spj+t3a
       pispj=pispj+t3
      endif
c
c Fock matrix
      ii=i+l1-1
      kk=k+l2-1
      term=pi0spj*ccoef
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
      t1=ccoef*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
      t4=sspj/z2a
c
      do 385 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      t2=Q(l3)-r(Nuc(i),l3)
      t2a=W(l3)-Q(l3)
      tx=r(Nuc(i),l3)-r(Nuc(j),l3)
      psd=t1*pispj
      dsp=t2*pi0spj+t2a*pispj
c
      if (l1.eq.l3) then
       psd=psd+t4
       dsp=dsp+t10
       f(Nuc(i),l3)=f(Nuc(i),l3)-(te+tee)*ss0pj
      endif
c
      if (l2.eq.l3) then
       psd=psd+t5
       dsp=dsp+t9
       f(Nucd(k),l3)=f(Nucd(k),l3)-te*ps
      endif
c
      ppp=dsp+tx*pi0spj
c
      f(Nuc(i),l3)=f(Nuc(i),l3)+a(i,ni)*(ty+tye)*dsp
      f(Nuc(j),l3)=f(Nuc(j),l3)+a(j,nj)*(ty+tye)*ppp
      f(Nucd(k),l3)=f(Nucd(k),l3)+ad(k,nk)*ty*psd
 385  continue
 381  continue
 382  continue
 380  continue
c-------------------------------------------------------------
c (pp|p) and gradients
      do 390 i=ns+1,ns+np,3
      do 390 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 390 ni=1,ncont(i)
      do 390 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 392
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 391 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 391 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      zc=2.D0*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t6a=ss1s/z2a
      t6=ss2s/z2a
c
      do 395 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t8a=p1s/z2a
      t8=p2s/z2a
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 395 l2=1,lij
      t1=Q(l2)-r(Nuc(j),l2)
      t2=W(l2)-Q(l2)
      pijs=t1*ps+t2*p1s
      pij1s=t1*p1s+t2*p2s
      pij2s=t1*p2s+t2*p3s
      spjs=t1*ss1s+t2*ss2s
      sp2js=t1*ss2s+t2*ss3s
      t7=sp2js/z2a
      t7a=spjs/z2a
c
      if (l1.eq.l2) then
       pijs=pijs+t3
       pij1s=pij1s+t4
       pij2s=pij2s+t5
      endif
c
      t11=(pijs-ti*pij1s)/zc
      t12=pij1s/z2a
c
      do 395 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      ppp=t1*pij1s
      spp=t1*spjs
      psp=t1*p1s
      pp1p=t1*pij2s
      spjpk=t1*sp2js
      pispk=t1*p2s
c
      if (l1.eq.l3) then
       ppp=ppp+t7a
       pp1p=pp1p+t7
       pispk=pispk+t6
       psp=psp+t6a
      endif
c
      if (l2.eq.l3) then
       ppp=ppp+t8a
       pp1p=pp1p+t8
       spjpk=spjpk+t6
       spp=spp+t6a
      endif
c
c Fock matrix construction
c
      ii=i+l1-1
      jj=j+l2-1
      kk=k+l3-1
c
      term=ccoef*ppp
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
      t1=ccoef*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
      t9=spjpk/z2a
      t10=pispk/z2a
      t14=(psp-roz*pispk)/z2
      t15=(spp-roz*spjpk)/z2
c
c gradients ----------------
      do 395 l4=1,3
      t1=W(l4)-r(Nucd(k),l4)
      ppd=t1*pp1p
      t20=Q(l4)-r(Nuc(i),l4)
      t30=W(l4)-Q(l4)
      tx=r(Nuc(i),l4)-r(Nuc(j),l4)
      dpp=t20*ppp+t30*pp1p
c
      if (l1.eq.l4) then
       ppd=ppd+t9
       dpp=dpp+t15
       f(Nuc(i),l4)=f(Nuc(i),l4)-(te+tee)*spp
      endif
c
      if (l2.eq.l4) then
       ppd=ppd+t10
       dpp=dpp+t14
       f(Nuc(j),l4)=f(Nuc(j),l4)-(te+tee)*psp
      endif
c
      if (l3.eq.l4) then
       ppd=ppd+t11
       dpp=dpp+t12
       f(Nucd(k),l4)=f(Nucd(k),l4)-te*pijs
      endif
c
       pdp=dpp+tx*ppp
c
       f(Nuc(i),l4)=f(Nuc(i),l4)+a(i,ni)*(ty+tye)*dpp
       f(Nuc(j),l4)=f(Nuc(j),l4)+a(j,nj)*(ty+tye)*pdp
       f(Nucd(k),l4)=f(Nucd(k),l4)+ad(k,nk)*ty*ppd
c
 395   continue
 391   continue
 392   continue
 390   continue
c
c-------------------------------------------------------------
c
c (ds|p)
      do 400 i=ns+np+1,M,6
      do 400 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 400 ni=1,ncont(i)
      do 400 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 402
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 401 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 401 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3a=(sss-roz*ss1s)/z2
      t3=(ss1s-roz*ss2s)/z2
      t3b=(ss2s-roz*ss3s)/z2
      t7=ss1s/z2a
      t7b=ss2s/z2a
c
      do 405 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      t5=p1s/z2a
      p2s=t1*ss2s+t2*ss3s
      t5b=p2s/z2a
      p3s=t1*ss3s+t2*ss4s
c
      do 405 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      t6=pj1s/z2a
      t6b=pj2s/z2a
      d0s=t1*ps+t2*p1s
      ds=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
c
      f1=1.D0
      if (l1.eq.l2) then
       d2s=d2s+t3b
       ds=ds+t3
       d0s=d0s+t3a
       f1=sq3
      endif
c
      t14=ds/z2a
      t15=(d0s-ti*ds)/zc
c
      do 405 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      dsp=t1*ds
      ds1p=t1*d2s
      pi0p=t1*p1s
      pi1p=t1*p2s
      pj0p=t1*pj1s
      pj1p=t1*pj2s
c
      if (l1.eq.l3) then
       dsp=dsp+t6
       ds1p=ds1p+t6b
       pi0p=pi0p+t7
       pi1p=pi1p+t7b
      endif
c
      if (l2.eq.l3) then
       dsp=dsp+t5
       ds1p=ds1p+t5b
       pj0p=pj0p+t7
       pj1p=pj1p+t7b
      endif
c
      t10=(pj0p-roz*pj1p)/z2
      t11=pj1p/z2a
      t12=(pi0p-roz*pi1p)/z2
      t13=pi1p/z2a
c
      l12=Ll(l1)+l2
      ii=i+l12-1
      kk=k+l3-1
c    
      cc=ccoef/f1
      term=dsp*cc
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
c gradients
c
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
      do 405 l4=1,3
        t1=Q(l4)-r(Nuc(j),l4)
        t2=W(l4)-Q(l4)
        t2b=W(l4)-r(Nucd(k),l4)
        tx=r(Nuc(i),l4)-r(Nuc(j),l4)
c
        dpp=t1*dsp+t2*ds1p
        dsd=t2b*ds1p
c
       if (l1.eq.l4) then
         dpp=dpp+t10
         dsd=dsd+t11
         f(Nuc(i),l4)=f(Nuc(i),l4)-(te+tee)*pj0p
       endif
c
       if (l2.eq.l4) then
         dpp=dpp+t12
         dsd=dsd+t13
         f(Nuc(i),l4)=f(Nuc(i),l4)-(te+tee)*pi0p
       endif
c
       if (l3.eq.l4) then
         dpp=dpp+t14
         dsd=dsd+t15
         f(Nucd(k),l4)=f(Nucd(k),l4)-te*d0s
       endif
c
        fsp=dpp-tx*dsp
        f(Nuc(i),l4)=f(Nuc(i),l4)+(ty+tye)*a(i,ni)*fsp
        f(Nuc(j),l4)=f(Nuc(j),l4)+(ty+tye)*a(j,nj)*dpp
        f(Nucd(k),l4)=f(Nucd(k),l4)+ty*ad(k,nk)*dsd
 405   continue
 401   continue
 402   continue
 400   continue
c
c-------------------------------------------------------------
c
c (dp|p) and gradients
      do 410 i=ns+np+1,M,6
      do 410 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 410 ni=1,ncont(i)
      do 410 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 412
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 411 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 411 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      t3a=(sss-roz*ss1s)/z2
      t3=(ss1s-roz*ss2s)/z2
      t4=(ss2s-roz*ss3s)/z2
      t4b=(ss3s-roz*ss4s)/z2
c
      do 415 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      t5a=(ps-roz*p1s)/z2
      t5=(p1s-roz*p2s)/z2
      p3s=t1*ss3s+t2*ss4s
      t5b=(p2s-roz*p3s)/z2
      p4s=t1*ss4s+t2*ss5s
c
      t12=p1s/z2a
      t12b=p2s/z2a
c
      do 415 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      t6a=(pjs-roz*pj1s)/z2
      t6=(pj1s-roz*pj2s)/z2
      t6b=(pj2s-roz*pj3s)/z2
c
      t10=pj1s/z2a
      t10b=pj2s/z2a
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+t3a
       d1s=d1s+t3
       d2s=d2s+t4
       d3s=d3s+t4b
       f1=sq3
      endif
c
      t9=d1s/z2a
      t9b=d2s/z2a
c
      do 415 l3=1,3
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      s1pk=t1*ss1s+t2*ss2s
      s2pk=t1*ss2s+t2*ss3s
      d0p=t1*ds+t2*d1s
      d1p=t1*d1s+t2*d2s
      d2p=t1*d2s+t2*d3s
      pi1p=t1*p1s+t2*p2s
      pi2p=t1*p2s+t2*p3s
      pj1p=t1*pj1s+t2*pj2s
      pj2p=t1*pj2s+t2*pj3s
c   
      t11=s1pk/z2a
      t11b=s2pk/z2a
c
      if (l1.eq.l3) then
       d0p=d0p+t6a
       d1p=d1p+t6
       d2p=d2p+t6b
       pi1p=pi1p+t3
       pi2p=pi2p+t4
      endif
c
      if (l2.eq.l3) then
       d0p=d0p+t5a
       d1p=d1p+t5
       d2p=d2p+t5b
       pj1p=pj1p+t3
       pj2p=pj2p+t4
      endif
c
      t7=pi1p/z2a
      t7b=pi2p/z2a
      t8=pj1p/z2a
      t8b=pj2p/z2a
c
      t26=d1p/z2a
      t27=(d0p-ti*d1p)/zc
c
      do 415 l4=1,3
      t1=W(l4)-r(Nucd(k),l4)
      dpp=t1*d1p
      d1pp=t1*d2p
      pj0pp=t1*pj1p
      pj1pp=t1*pj2p
      pi0pp=t1*pi1p
      pi1pp=t1*pi2p
      d0pl=t1*d1s
      d1pl=t1*d2s
c
      if (l1.eq.l4) then
       dpp=dpp+t8
       d1pp=d1pp+t8b
       d0pl=d0pl+t10
       d1pl=d1pl+t10b
       pi0pp=pi0pp+t11
       pi1pp=pi1pp+t11b
      endif
c
      if (l2.eq.l4) then
       dpp=dpp+t7
       d1pp=d1pp+t7b
       d0pl=d0pl+t12
       d1pl=d1pl+t12b
       pj0pp=pj0pp+t11
       pj1pp=pj1pp+t11b
      endif
c
      if (l3.eq.l4) then
       dpp=dpp+t9
       d1pp=d1pp+t9b
       pi0pp=pi0pp+t12
       pi1pp=pi1pp+t12b
       pj0pp=pj0pp+t10
       pj1pp=pj1pp+t10b
      endif
c     
      t20=(pj0pp-roz*pj1pp)/z2
      t21=pj1pp/z2a
      t22=(pi0pp-roz*pi1pp)/z2
      t23=pi1pp/z2a
      t24=(d0pl-roz*d1pl)/z2
      t25=d1pl/z2a
c
      l12=Ll(l1)+l2
      ii=i+l12-1
      jj=j+l3-1
      kk=k+l4-1
c
      cc=ccoef/f1
      term=dpp*cc
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c gradients
c
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
      do 415 l5=1,3
        t1=Q(l5)-r(Nuc(j),l5)
        t2=W(l5)-Q(l5)
        t2b=W(l5)-r(Nucd(k),l5)
        tx=r(Nuc(i),l5)-r(Nuc(j),l5)
c
        ddp=t1*dpp+t2*d1pp
        dpd=t2b*d1pp
c
       if (l1.eq.l5) then
        ddp=ddp+t20
        dpd=dpd+t21
        f(Nuc(i),l5)=f(Nuc(i),l5)-(te+tee)*pj0pp
       endif
c
       if (l2.eq.l5) then
        ddp=ddp+t22
        dpd=dpd+t23
        f(Nuc(i),l5)=f(Nuc(i),l5)-(te+tee)*pi0pp
       endif
c
       if (l3.eq.l5) then
        ddp=ddp+t24
        dpd=dpd+t25
        f(Nuc(j),l5)=f(Nuc(j),l5)-(te+tee)*d0pl
       endif
c
       if (l4.eq.l5) then
        ddp=ddp+t26
        dpd=dpd+t27
        f(Nucd(k),l5)=f(Nucd(k),l5)-te*d0p
       endif
c
       fpp=ddp-tx*dpp
       f(Nuc(i),l5)=f(Nuc(i),l5)+a(i,ni)*(ty+tye)*fpp
       f(Nuc(j),l5)=f(Nuc(j),l5)+a(j,nj)*(ty+tye)*ddp
       f(Nucd(k),l5)=f(Nucd(k),l5)+ad(k,nk)*ty*dpd
       
 415   continue
 411   continue
 412   continue
 410   continue
c
c-------------------------------------------------------------
c
c (dd|p) and gradients
      do 420 i=ns+np+1,M,6
      do 420 j=ns+np+1,i,6
c
      dd=d(Nuc(i),Nuc(j))
c
      do 420 ni=1,ncont(i)
      do 420 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 422
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 421 k=nsd+1,nsd+npd,3
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 421 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      ss6s=t2*FUNCT(6,u)
      ta=(sss-roz*ss1s)/z2
      t3=(ss1s-roz*ss2s)/z2
      t4=(ss2s-roz*ss3s)/z2
      t5=(ss3s-roz*ss4s)/z2
      t5b=(ss4s-roz*ss5s)/z2
c
      do 425 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      p5s=t1*ss5s+t2*ss6s
      t6a=(ps-roz*p1s)/z2
      t6=(p1s-roz*p2s)/z2
      t7=(p2s-roz*p3s)/z2
      t7b=(p3s-roz*p4s)/z2
c
      do 425 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      pj4s=t1*ss4s+t2*ss5s
      t8a=(pjs-roz*pj1s)/z2
      t8=(pj1s-roz*pj2s)/z2
      t9=(pj2s-roz*pj3s)/z2
      t9b=(pj3s-roz*pj4s)/z2
c
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
      d4s=t1*p4s+t2*p5s
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+ta
       d1s=d1s+t3
       d2s=d2s+t4
       d3s=d3s+t5
       d4s=d4s+t5b
       f1=sq3
      endif
c
      t18a=(ds-roz*d1s)/z2
      t18=(d1s-roz*d2s)/z2
      t18b=(d2s-roz*d3s)/z2
      t35=d1s/z2a
      t35b=d2s/z2a
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 425 l3=1,lij
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      d0pk=t1*ds+t2*d1s
      d1pk=t1*d1s+t2*d2s
      d2pk=t1*d2s+t2*d3s
      d3pk=t1*d3s+t2*d4s
      pjpk=t1*pjs+t2*pj1s
      pj1pk=t1*pj1s+t2*pj2s
      pj2pk=t1*pj2s+t2*pj3s
      pj3pk=t1*pj3s+t2*pj4s
      pipk=t1*ps+t2*p1s
      pi1pk=t1*p1s+t2*p2s
      pi2pk=t1*p2s+t2*p3s
      pi3pk=t1*p3s+t2*p4s
      spk=t1*sss+t2*ss1s
      s1pk=t1*ss1s+t2*ss2s
      s2pk=t1*ss2s+t2*ss3s
      s3pk=t1*ss3s+t2*ss4s
      s4pk=t1*ss4s+t2*ss5s
      t10=(s1pk-roz*s2pk)/z2
      t10b=(s2pk-roz*s3pk)/z2
c
c
      if (l1.eq.l3) then
       d0pk=d0pk+t8a
       d1pk=d1pk+t8
       d2pk=d2pk+t9
       d3pk=d3pk+t9b
       pipk=pipk+ta
       pi1pk=pi1pk+t3
       pi2pk=pi2pk+t4
       pi3pk=pi3pk+t5
      endif
c
      if (l2.eq.l3) then
       d0pk=d0pk+t6a
       d1pk=d1pk+t6
       d2pk=d2pk+t7
       d3pk=d3pk+t7b
       pjpk=pjpk+ta
       pj1pk=pj1pk+t3
       pj2pk=pj2pk+t4
       pj3pk=pj3pk+t5
      endif
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      t16a=(pjpk-roz*pj1pk)/z2
      t16=(pj1pk-roz*pj2pk)/z2
      t16b=(pj2pk-roz*pj3pk)/z2
      t17a=(pipk-roz*pi1pk)/z2
      t17=(pi1pk-roz*pi2pk)/z2
      t17b=(pi2pk-roz*pi3pk)/z2
      t30=pi1pk/z2a
      t30b=pi2pk/z2a
      t31=pj1pk/z2a
      t31b=pj2pk/z2a
c
      do 425 l4=1,lk
      t1=Q(l4)-r(Nuc(j),l4)
      t2=W(l4)-Q(l4)
      d0d=t1*d0pk+t2*d1pk
      d1d=t1*d1pk+t2*d2pk
      d2d=t1*d2pk+t2*d3pk
      pi1pl=t1*p1s+t2*p2s
      pi2pl=t1*p2s+t2*p3s
      pj1pl=t1*pj1s+t2*pj2s
      pj2pl=t1*pj2s+t2*pj3s
c     
      pjdkl=t1*pj1pk+t2*pj2pk
      pj2dkl=t1*pj2pk+t2*pj3pk
      pidkl=t1*pi1pk+t2*pi2pk
      pi2dkl=t1*pi2pk+t2*pi3pk
      d1pl=t1*d1s+t2*d2s
      d2pl=t1*d2s+t2*d3s
c
      s1ds=t1*s1pk+t2*s2pk
      s2ds=t1*s2pk+t2*s3pk
c
      if (l1.eq.l4) then
       d0d=d0d+t16a
       d1d=d1d+t16
       d2d=d2d+t16b
       pidkl=pidkl+t10
       pi2dkl=pi2dkl+t10b
       d1pl=d1pl+t8
       d2pl=d2pl+t9
       pi1pl=pi1pl+t3
       pi2pl=pi2pl+t4
      endif
c
      if (l2.eq.l4) then
       d0d=d0d+t17a
       d1d=d1d+t17
       d2d=d2d+t17b
       pjdkl=pjdkl+t10
       pj2dkl=pj2dkl+t10b
       d1pl=d1pl+t6
       d2pl=d2pl+t7
       pj1pl=pj1pl+t3
       pj2pl=pj2pl+t4
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       d0d=d0d+t18a
       d1d=d1d+t18
       d2d=d2d+t18b
       pjdkl=pjdkl+t8
       pj2dkl=pj2dkl+t9
       pidkl=pidkl+t6
       pi2dkl=pi2dkl+t7
       s1ds=s1ds+t3
       s2ds=s2ds+t4
       f2=sq3
      endif
c
      t11=pjdkl/z2a
      t11b=pj2dkl/z2a
      t12=pidkl/z2a
      t12b=pi2dkl/z2a
      t13=d1pl/z2a
      t13b=d2pl/z2a
      t14=d1pk/z2a
      t14b=d2pk/z2a
      t32=pi1pl/z2a
      t32b=pi2pl/z2a
      t33=pj1pl/z2a
      t33b=pj2pl/z2a
      t34=s1ds/z2a
      t34b=s2ds/z2a
      t28=d1d/z2a
      t29=(d0d-ti*d1d)/zc
c
      do 425 l5=1,3
c
      t1=W(l5)-r(Nucd(k),l5)
      ddp=t1*d1d
      dd1p=t1*d2d
      pj0dp=t1*pjdkl
      pj1dp=t1*pj2dkl
      pi0dp=t1*pidkl
      pi1dp=t1*pi2dkl
      d0plp=t1*d1pl
      d1plp=t1*d2pl
      d0pkp=t1*d1pk
      d1pkp=t1*d2pk
c
      if (l1.eq.l5) then
       ddp=ddp+t11
       dd1p=dd1p+t11b
       pi0dp=pi0dp+t34
       pi1dp=pi1dp+t34b
       d0plp=d0plp+t33
       d1plp=d1plp+t33b
       d0pkp=d0pkp+t31
       d1pkp=d1pkp+t31b
      endif
c
      if (l2.eq.l5) then
       ddp=ddp+t12
       dd1p=dd1p+t12b
       pj0dp=pj0dp+t34
       pj1dp=pj1dp+t34b
       d0plp=d0plp+t32
       d1plp=d1plp+t32b
       d0pkp=d0pkp+t30
       d1pkp=d1pkp+t30b
      endif
c
      if (l3.eq.l5) then
       ddp=ddp+t13
       dd1p=dd1p+t13b
       pj0dp=pj0dp+t33
       pj1dp=pj1dp+t33b
       pi0dp=pi0dp+t32
       pi1dp=pi1dp+t32b
       d0pkp=d0pkp+t35
       d1pkp=d1pkp+t35b
      endif
c
      if (l4.eq.l5) then
       ddp=ddp+t14
       dd1p=dd1p+t14b
       pj0dp=pj0dp+t31
       pj1dp=pj1dp+t31b
       pi0dp=pi0dp+t30
       pi1dp=pi1dp+t30b
       d0plp=d0plp+t35
       d1plp=d1plp+t35b
      endif
c
      t20=(pj0dp-roz*pj1dp)/z2
      t21=pj1dp/z2a
      t22=(pi0dp-roz*pi1dp)/z2
      t23=pi1dp/z2a
      t24=(d0plp-roz*d1plp)/z2
      t25=d1plp/z2a
      t26=(d0pkp-roz*d1pkp)/z2
      t27=d1pkp/z2a
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
      ii=i+l12-1
      jj=j+l34-1
      kk=k+l5-1
c
      cc=ccoef/(f1*f2)
      term=ddp*cc
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
c gradients
c
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
      do 425 l6=1,3
       t1=Q(l6)-r(Nuc(i),l6)
       t2=W(l6)-Q(l6)
       t2b=W(l6)-r(Nucd(k),l6)
       tx=r(Nuc(i),l6)-r(Nuc(j),l6)
c
       fdp=t1*ddp+t2*dd1p
       ddd=t2b*dd1p
c
      if (l1.eq.l6) then
       fdp=fdp+t20
       ddd=ddd+t21
       f(Nuc(i),l6)=f(Nuc(i),l6)-(te+tee)*pj0dp
      endif
c
      if (l2.eq.l6) then
       fdp=fdp+t22
       ddd=ddd+t23
       f(Nuc(i),l6)=f(Nuc(i),l6)-(te+tee)*pi0dp
      endif
c
      if (l3.eq.l6) then
       fdp=fdp+t24
       ddd=ddd+t25
       f(Nuc(j),l6)=f(Nuc(j),l6)-(te+tee)*d0plp
      endif
c
      if (l4.eq.l6) then
       fdp=fdp+t26
       ddd=ddd+t27
       f(Nuc(j),l6)=f(Nuc(j),l6)-(te+tee)*d0pkp
      endif
c
      if (l5.eq.l6) then
       fdp=fdp+t28
       ddd=ddd+t29
       f(Nucd(k),l6)=f(Nucd(k),l6)-te*d0d
      endif
c
      dfp=fdp+tx*ddp
      f(Nuc(i),l6)=f(Nuc(i),l6)+a(i,ni)*(ty+tye)*fdp
      f(Nuc(j),l6)=f(Nuc(j),l6)+a(j,nj)*(ty+tye)*dfp
      f(Nucd(k),l6)=f(Nucd(k),l6)+ad(k,nk)*ty*ddd
c
 425   continue
 421   continue
 422   continue
 420   continue
c
c-------------------------------------------------------------
c
c (ss|d) and gradients
      do 430 i=1,ns
      do 430 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
      kn=i+Jx(j)
c
      do 430 ni=1,ncont(i)
      do 430 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 432
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 431 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 431 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=ti
c
      zc=2.D0*ad(k,nk)
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ta=(sss-roz*ss1s)/zc
      tb=(ss1s-roz*ss2s)/zc
c
      do 435 l1=1,3
      t1=W(l1)-r(Nucd(k),l1)
      ss0p=t1*ss1s
      ss1p=t1*ss2s
      t13=ss1p/z2
      t11=(ss0p-roz*ss1p)/zc
      ss2p=t1*ss3s
c
      do 435 l2=1,l1
      t1=W(l2)-r(Nucd(k),l2)
      ss0pj=t1*ss1s
      ss1pj=t1*ss2s
      t12=ss1pj/z2
      t10=(ss0pj-roz*ss1pj)/zc
      ss0d=t1*ss1p
      ss1d=t1*ss2p
c
      f1=1.D0
      if (l1.eq.l2) then
       ss0d=ss0d+ta
       ss1d=ss1d+tb
       f1=sq3
      endif
c
      cc=ccoef/f1
      term=ss0d*cc
c
      l12=Ll(l1)+l2
      kk=k+l12-1
c
c
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
       do 435 l3=1,3
       t1=Q(l3)-r(Nuc(i),l3)
       t2=W(l3)-Q(l3)
       t3=W(l3)-r(Nucd(k),l3)
       
       tx=r(Nuc(i),l3)-r(Nuc(j),l3)
c
       ssf=t3*ss1d
       psd=t1*ss0d+t2*ss1d
c
       if (l1.eq.l3) then
        ssf=ssf+t10
        psd=psd+t12
        f(Nucd(k),l3)=f(Nucd(k),l3)-te*ss0pj
       endif
c
       if (l2.eq.l3) then
        ssf=ssf+t11
        psd=psd+t13
        f(Nucd(k),l3)=f(Nucd(k),l3)-te*ss0p
       endif
c
       spd=psd+tx*ss0d
c
       f(Nuc(i),l3)=f(Nuc(i),l3)+a(i,ni)*(ty+tye)*psd
       f(Nuc(j),l3)=f(Nuc(j),l3)+a(j,nj)*(ty+tye)*spd
       f(Nucd(k),l3)=f(Nucd(k),l3)+ad(k,nk)*ty*ssf

 435   continue
 431   continue
 432   continue
 430   continue
c
c-------------------------------------------------------------
c
c (ps|d) and gradients
      do 440 i=ns+1,ns+np,3
      do 440 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 440 ni=1,ncont(i)
      do 440 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 442
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 441 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 441 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      zc=2.D0*ad(k,nk)
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
c
      roz=ti
      ro=roz*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      t3a=ss1s/z2a
      t3=ss2s/z2a
      t3b=ss3s/z2a
      t7=(sss-roz*ss1s)/zc
      t7b=(ss1s-roz*ss2s)/zc
c
      do 445 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      t5=(ps-roz*p1s)/zc
      t5b=(p1s-roz*p2s)/zc
c
      do 445 l2=1,3
      t1=W(l2)-r(Nucd(k),l2)
      sspj=t1*ss2s
      ss2pj=t1*ss3s
      pi0spj=t1*p1s
      pispj=t1*p2s
      pi2spj=t1*p3s
c
      if (l1.eq.l2) then
       pi2spj=pi2spj+t3b
       pispj=pispj+t3
       pi0spj=pi0spj+t3a
      endif
c
      t4=sspj/z2a
      t4b=ss2pj/z2a
      t14=pispj/z2a
      t15=(pi0spj-roz*pispj)/zc
c
      do 445 l3=1,l2
      t1=W(l3)-r(Nucd(k),l3)
      ss0d=t1*sspj
      ss1d=t1*ss2pj
      ps0d=t1*pispj
      ps1d=t1*pi2spj
      p0pk=t1*p1s
      p1pk=t1*p2s
c
      if (l1.eq.l3) then
       ps0d=ps0d+t4
       ps1d=ps1d+t4b
       p0pk=p0pk+t3a
       p1pk=p1pk+t3
      endif
c
      f1=1.
      if (l2.eq.l3) then
       ss0d=ss0d+t7
       ss1d=ss1d+t7b
       ps0d=ps0d+t5
       ps1d=ps1d+t5b
       f1=sq3
      endif
c
      t10=(ss0d-tj*ss1d)/z2
      t11=ss1d/z2a
      t12=p1pk/z2a
      t13=(p0pk-roz*p1pk)/zc
c
      l23=l2*(l2-1)/2+l3
      ii=i+l1-1
      kk=k+l23-1
c     
      cc=ccoef/f1
      term=ps0d*cc
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
       do 445 l4=1,3
        t1=Q(l4)-r(Nuc(i),l4)
        t2=W(l4)-Q(l4)
        t2b=W(l4)-r(Nucd(k),l4)
        tx=r(Nuc(i),l4)-r(Nuc(j),l4)
c
        dsd=t1*ps0d+t2*ps1d
        psf=t2b*ps1d
c
       if (l1.eq.l4) then
        dsd=dsd+t10
        psf=psf+t11
        f(Nuc(i),l4)=f(Nuc(i),l4)-(te+tee)*ss0d
       endif
c
       if (l2.eq.l4) then
        dsd=dsd+t12
        psf=psf+t13
        f(Nucd(k),l4)=f(Nucd(k),l4)-te*p0pk
       endif
c
       if (l3.eq.l4) then
        dsd=dsd+t14
        psf=psf+t15
        f(Nucd(k),l4)=f(Nucd(k),l4)-te*pi0spj
       endif
c
       ppd=dsd+tx*ps0d
c
       f(Nuc(i),l4)=f(Nuc(i),l4)+a(i,ni)*(ty+tye)*dsd
       f(Nuc(j),l4)=f(Nuc(j),l4)+a(j,nj)*(ty+tye)*ppd
       f(Nucd(k),l4)=f(Nucd(k),l4)+ad(k,nk)*ty*psf
 445  continue
 441  continue
 442  continue
 440  continue
c
c-------------------------------------------------------------
c
c (pp|d) and gradients
      do 450 i=ns+1,ns+np,3
      do 450 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 450 ni=1,ncont(i)
      do 450 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 452
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 451 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 451 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      zc=2.D0*ad(k,nk)
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t5b=(ss3s-roz*ss4s)/z2
      t6=ss2s/z2a
      t6b=ss3s/z2a
c
      do 455 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      t8a=p1s/z2a
      t8=p2s/z2a
      t8b=p3s/z2a
      t33=(ps-ti*p1s)/zc
      t33b=(p1s-ti*p2s)/zc
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 455 l2=1,lij
      t1=Q(l2)-r(Nuc(j),l2)
      t2=W(l2)-Q(l2)
      pijs=t1*ps+t2*p1s
      pij1s=t1*p1s+t2*p2s
      pij2s=t1*p2s+t2*p3s
      pij3s=t1*p3s+t2*p4s
      sp0js=t1*sss+t2*ss1s
      spjs=t1*ss1s+t2*ss2s
      sp2js=t1*ss2s+t2*ss3s
      sp3js=t1*ss3s+t2*ss4s
      t7a=spjs/z2a
      t7=sp2js/z2a
      t7b=sp3js/z2a
c
      if (l1.eq.l2) then
       pijs=pijs+t3
       pij1s=pij1s+t4
       pij2s=pij2s+t5
       pij3s=pij3s+t5b
      endif
c
      t11=(pijs-ti*pij1s)/zc
      t11b=(pij1s-ti*pij2s)/zc
      t31=(sp0js-ti*spjs)/zc
      t31b=(spjs-ti*sp2js)/zc
c
      do 455 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      pp0p=t1*pij1s
      pp1p=t1*pij2s
      pp2p=t1*pij3s
      spjpk=t1*sp2js
      s2pjpk=t1*sp3js
      pispk=t1*p2s
      pi2spk=t1*p3s
      sspk=t1*ss2s
      ss2pk=t1*ss3s
      t30=sspk/z2a
      t30b=ss2pk/z2a
c
      if (l1.eq.l3) then
       pp0p=pp0p+t7a
       pp1p=pp1p+t7
       pp2p=pp2p+t7b
       pispk=pispk+t6
       pi2spk=pi2spk+t6b
      endif
c
      if (l2.eq.l3) then
       pp0p=pp0p+t8a
       pp1p=pp1p+t8
       pp2p=pp2p+t8b
       spjpk=spjpk+t6
       s2pjpk=s2pjpk+t6b
      endif
c
      t9=spjpk/z2a
      t9b=s2pjpk/z2a
      t10=pispk/z2a
      t10b=pi2spk/z2a
      t26=pp1p/z2a
      t27=(pp0p-ti*pp1p)/zc
c
      do 455 l4=1,l3
      t1=W(l4)-r(Nucd(k),l4)
      pp0d=t1*pp1p
      pp1d=t1*pp2p
      sp0d=t1*spjpk
      sp1d=t1*s2pjpk
      ps0d=t1*pispk
      ps1d=t1*pi2spk
      pp0pl=t1*pij1s
      pp1pl=t1*pij2s
c
      if (l1.eq.l4) then
       pp0d=pp0d+t9
       pp1d=pp1d+t9b
       ps0d=ps0d+t30
       ps1d=ps1d+t30b
       pp0pl=pp0pl+t7a
       pp1pl=pp1pl+t7
      endif
c
      if (l2.eq.l4) then
       pp0d=pp0d+t10
       pp1d=pp1d+t10b
       sp0d=sp0d+t30
       sp1d=sp1d+t30b
       pp0pl=pp0pl+t8a
       pp1pl=pp1pl+t8
      endif
c
      f1=1.D0
      if (l3.eq.l4) then
       pp0d=pp0d+t11
       pp1d=pp1d+t11b
       sp0d=sp0d+t31
       sp1d=sp1d+t31b
       ps0d=ps0d+t33
       ps1d=ps1d+t33b
       f1=sq3
      endif
c
      ii=i+l1-1
      jj=j+l2-1
      l34=l3*(l3-1)/2+l4
      kk=k+l34-1
c
      cc=ccoef/f1
      term=pp0d*cc
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
      t20=(sp0d-sp1d*roz)/z2
      t21=sp1d/z2a
      t22=(ps0d-roz*ps1d)/z2
      t23=ps1d/z2a
      t24=pp1pl/z2a
      t25=(pp0pl-ti*pp1pl)/zc

c gradients
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
       do 455 l5=1,3
        t1=Q(l5)-r(Nuc(i),l5)
        t2=W(l5)-Q(l5)
        t2b=W(l5)-r(Nucd(k),l5)
        tx=r(Nuc(i),l5)-r(Nuc(j),l5)
c
        dpd=t1*pp0d+t2*pp1d
        ppf=t2b*pp1d
c
        if (l1.eq.l5) then
        dpd=dpd+t20
        ppf=ppf+t21
        f(Nuc(i),l5)=f(Nuc(i),l5)-(te+tee)*sp0d
       endif
c
       if (l2.eq.l5) then
        dpd=dpd+t22
        ppf=ppf+t23
        f(Nuc(j),l5)=f(Nuc(j),l5)-(te+tee)*ps0d
       endif
c
       if (l3.eq.l5) then
        dpd=dpd+t24
        ppf=ppf+t25
        f(Nucd(k),l5)=f(Nucd(k),l5)-te*pp0pl
       endif
c
       if (l4.eq.l5) then
        dpd=dpd+t26
        ppf=ppf+t27
        f(Nucd(k),l5)=f(Nucd(k),l5)-te*pp0p
       endif
c
       pdd=dpd+tx*pp0d
c
       f(Nuc(i),l5)=f(Nuc(i),l5)+a(i,ni)*(ty+tye)*dpd
       f(Nuc(j),l5)=f(Nuc(j),l5)+a(j,nj)*(ty+tye)*pdd
       f(Nucd(k),l5)=f(Nucd(k),l5)+ad(k,nk)*ty*ppf
 455   continue
 451   continue
 452   continue
 450   continue
c
c-------------------------------------------------------------
c
c (ds|d) and gradients
      do 460 i=ns+np+1,M,6
      do 460 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
      k1=Jx(j)
c
      do 460 ni=1,ncont(i)
      do 460 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 462
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 461 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 461 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t5b=(ss3s-roz*ss4s)/z2
      t6=ss2s/z2a
      t6b=ss3s/z2a
c
      do 465 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      t7a=p1s/z2a
      t7=p2s/z2a
      t7b=p3s/z2a
      t16=(ps-ti*p1s)/zc
      t16b=(p1s-ti*p2s)/zc
c
      do 465 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
      t8a=pj1s/z2a
      t8=pj2s/z2a
      t8b=pj3s/z2a
      t15=(pjs-ti*pj1s)/zc
      t15b=(pj1s-ti*pj2s)/zc
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t5
       d3s=d3s+t5b
       f1=sq3
      endif
c
      t11=(ds-ti*d1s)/zc
      t11b=(d1s-ti*d2s)/zc
c
      do 465 l3=1,3
      t1=W(l3)-r(Nucd(k),l3)
      ds0p=t1*d1s
      ds1p=t1*d2s
      ds2p=t1*d3s
      pis1pk=t1*p2s
      pis2pk=t1*p3s
      pjs1pk=t1*pj2s
      pjs2pk=t1*pj3s
      ss1pk=t1*ss2s
      ss2pk=t1*ss3s
      t12=ss1pk/z2a
      t12b=ss2pk/z2a
c
      if (l1.eq.l3) then
       ds0p=ds0p+t8a
       ds1p=ds1p+t8
       ds2p=ds2p+t8b
       pis1pk=pis1pk+t6
       pis2pk=pis2pk+t6b
      endif
c
      if (l2.eq.l3) then
       ds0p=ds0p+t7a
       ds1p=ds1p+t7
       ds2p=ds2p+t7b
       pjs1pk=pjs1pk+t6
       pjs2pk=pjs2pk+t6b
      endif
c
      t9=pjs1pk/z2a
      t9b=pjs2pk/z2a
      t10=pis1pk/z2a
      t10b=pis2pk/z2a
      t26=ds1p/z2a
      t27=(ds0p-ti*ds1p)/zc
c
      do 465 l4=1,l3
       t1=W(l4)-r(Nucd(k),l4)
       dsd=t1*ds1p
       ds1d=t1*ds2p
       pj0sd=t1*pjs1pk
       pj1sd=t1*pjs2pk
       pi0sd=t1*pis1pk
       pi1sd=t1*pis2pk
       d0pl=t1*d1s
       d1pl=t1*d2s
c
      if (l1.eq.l4) then
       dsd=dsd+t9
       ds1d=ds1d+t9b
       pi0sd=pi0sd+t12
       pi1sd=pi1sd+t12b
       d0pl=d0pl+t8a
       d1pl=d1pl+t8
      endif
c
      if (l2.eq.l4) then
       dsd=dsd+t10
       ds1d=ds1d+t10b
       pj0sd=pj0sd+t12
       pj1sd=pj1sd+t12b
       d0pl=d0pl+t7a
       d1pl=d1pl+t7
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       dsd=dsd+t11
       ds1d=ds1d+t11b
       pj0sd=pj0sd+t15
       pj1sd=pj1sd+t15b
       pi0sd=pi0sd+t16
       pi1sd=pi1sd+t16b
       f2=sq3
      endif
c
      t20=(pj0sd-roz*pj1sd)/z2
      t21=pj1sd/z2a
      t22=(pi0sd-roz*pi1sd)/z2
      t23=pi1sd/z2a
      t24=d1pl/z2a
      t25=(d0pl-ti*d1pl)/zc
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
c
      ii=i+l12-1
      kk=k+l34-1
c
      cc=ccoef/(f1*f2)
      term=dsd*cc
c  
c
      kn=ii+k1
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c
c gradients
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
       do 465 l5=1,3
        t1=Q(l5)-r(Nuc(j),l5)
        t2=W(l5)-Q(l5)
        t2b=W(l5)-r(Nucd(k),l5)
        tx=r(Nuc(i),l5)-r(Nuc(j),l5)
c
        dpd=t1*dsd+t2*ds1d
        dsf=t2b*ds1d
c
        if (l1.eq.l5) then
        dpd=dpd+t20
        dsf=dsf+t21
        f(Nuc(i),l5)=f(Nuc(i),l5)-(te+tee)*pj0sd
       endif
c
       if (l2.eq.l5) then
        dpd=dpd+t22
        dsf=dsf+t23
        f(Nuc(i),l5)=f(Nuc(i),l5)-(te+tee)*pi0sd
       endif
c
       if (l3.eq.l5) then
        dpd=dpd+t24
        dsf=dsf+t25
        f(Nucd(k),l5)=f(Nucd(k),l5)-te*d0pl
       endif
c
       if (l4.eq.l5) then
        dpd=dpd+t26
        dsf=dsf+t27
        f(Nucd(k),l5)=f(Nucd(k),l5)-te*ds0p
       endif
c
       fsd=dpd-tx*dsd
c
       f(Nuc(i),l5)=f(Nuc(i),l5)+a(i,ni)*(ty+tye)*fsd
       f(Nuc(j),l5)=f(Nuc(j),l5)+a(j,nj)*(ty+tye)*dpd
       f(Nucd(k),l5)=f(Nucd(k),l5)+ad(k,nk)*ty*dsf
c
 465   continue
 461   continue
 462   continue
 460   continue
c
c-------------------------------------------------------------
c
c (dp|d) and gradients
      do 470 i=ns+np+1,M,6
      do 470 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 470 ni=1,ncont(i)
      do 470 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*dd
       if (rexp.gt.rmax) goto 472
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 471 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 471 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      ss6s=t2*FUNCT(6,u)
c
      ta=(sss-roz*ss1s)/z2
      t3=(ss1s-roz*ss2s)/z2
      t4=(ss2s-roz*ss3s)/z2
      t5=(ss3s-roz*ss4s)/z2
      t5b=(ss4s-roz*ss5s)/z2
      t5x=ss2s/z2a
      t5y=ss3s/z2a
c
      do 475 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      p5s=t1*ss5s+t2*ss6s
      t6a=(ps-roz*p1s)/z2
      t6b=(p1s-roz*p2s)/z2
      t6c=(p2s-roz*p3s)/z2
      t6d=(p3s-roz*p4s)/z2
      t9=p2s/z2a
      t9b=p3s/z2a
c
      do 475 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
      d4s=t1*p4s+t2*p5s
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      pj4s=t1*ss4s+t2*ss5s
      t7a=(pjs-roz*pj1s)/z2
      t7b=(pj1s-roz*pj2s)/z2
      t7c=(pj2s-roz*pj3s)/z2
      t7d=(pj3s-roz*pj4s)/z2
c
      t8=pj2s/z2a
      t8b=pj3s/z2a
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+ta
       d1s=d1s+t3
       d2s=d2s+t4
       d3s=d3s+t5
       d4s=d4s+t5b
       f1=sq3
      endif
c
      t10a=d1s/z2a
      t10=d2s/z2a
      t10b=d3s/z2a
      t24=(ds-ti*d1s)/zc
      t24b=(d1s-ti*d2s)/zc
      t28=d1s/z2a
      t28b=d2s/z2a
c
      do 475 l3=1,3
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      s2pks=t1*ss2s+t2*ss3s
      s3pks=t1*ss3s+t2*ss4s
      dp=t1*ds+t2*d1s
      d1p=t1*d1s+t2*d2s
      d2p=t1*d2s+t2*d3s
      d3p=t1*d3s+t2*d4s
      pi0p=t1*ps+t2*p1s
      pi1p=t1*p1s+t2*p2s
      pi2p=t1*p2s+t2*p3s
      pi3p=t1*p3s+t2*p4s
      pj0p=t1*pjs+t2*pj1s
      pj1p=t1*pj1s+t2*pj2s
      pj2p=t1*pj2s+t2*pj3s
      pj3p=t1*pj3s+t2*pj4s
c
      if (l1.eq.l3) then
       dp=dp+t7a
       d1p=d1p+t7b
       d2p=d2p+t7c
       d3p=d3p+t7d
       pi0p=pi0p+ta
       pi1p=pi1p+t3
       pi2p=pi2p+t4
       pi3p=pi3p+t5
      endif
c
      if (l2.eq.l3) then
       dp=dp+t6a
       d1p=d1p+t6b
       d2p=d2p+t6c
       d3p=d3p+t6d
       pj0p=pj0p+ta
       pj1p=pj1p+t3
       pj2p=pj2p+t4
       pj3p=pj3p+t5
      endif
c
      t11a=pi1p/z2a
      t11=pi2p/z2a
      t11b=pi3p/z2a
      t12a=pj1p/z2a
      t12=pj2p/z2a
      t12b=pj3p/z2a
      t13=s2pks/z2a
      t13b=s3pks/z2a
      t22=(pj0p-ti*pj1p)/zc
      t22b=(pj1p-ti*pj2p)/zc
      t26=pj1p/z2a
      t26b=pj2p/z2a
      t25=(pi0p-ti*pi1p)/zc
      t25b=(pi1p-ti*pi2p)/zc
      t27=pi1p/z2a
      t27b=pi2p/z2a
c
      do 475 l4=1,3
      t1=W(l4)-r(Nucd(k),l4)
      dp0p=t1*d1p
      dp1p=t1*d2p
      dp2p=t1*d3p
      pjpkpl=t1*pj2p
      pj2pkpl=t1*pj3p
      pipkpl=t1*pi2p
      pi2pkpl=t1*pi3p
      dspl=t1*d2s
      ds2pl=t1*d3s
      pj1spl=t1*pj2s
      pj2spl=t1*pj3s
      pi1spl=t1*p2s
      pi2spl=t1*p3s
      s1pkpl=t1*s2pks
      s2pkpl=t1*s3pks
c
      if (l1.eq.l4) then
       dp0p=dp0p+t12a
       dp1p=dp1p+t12
       dp2p=dp2p+t12b
       pipkpl=pipkpl+t13
       pi2pkpl=pi2pkpl+t13b
       ds2pl=ds2pl+t8b
       dspl=dspl+t8
       pi1spl=pi1spl+t5x
       pi2spl=pi2spl+t5y
      endif
c
      if (l2.eq.l4) then
       dp0p=dp0p+t11a
       dp1p=dp1p+t11
       dp2p=dp2p+t11b
       pjpkpl=pjpkpl+t13
       pj2pkpl=pj2pkpl+t13b
       dspl=dspl+t9
       ds2pl=ds2pl+t9b
       pj1spl=pj1spl+t5x
       pj2spl=pj2spl+t5y
      endif
c
      if (l3.eq.l4) then
       dp0p=dp0p+t10a
       dp1p=dp1p+t10
       dp2p=dp2p+t10b
       pipkpl=pipkpl+t9
       pi2pkpl=pi2pkpl+t9b
       pjpkpl=pjpkpl+t8
       pj2pkpl=pj2pkpl+t8b
       s1pkpl=s1pkpl+t5x
       s2pkpl=s2pkpl+t5y
      endif
c
      t14=pjpkpl/z2a
      t14b=pj2pkpl/z2a
      t15=pipkpl/z2a
      t15b=pi2pkpl/z2a
      t16=dspl/z2a
      t16b=ds2pl/z2a
      t17=(dp-ti*d1p)/zc
      t17b=(d1p-ti*d2p)/zc
      t20=s1pkpl/z2a
      t20b=s2pkpl/z2a
      t21=pj1spl/z2a
      t21b=pj2spl/z2a
      t23=pi1spl/z2a
      t23b=pi2spl/z2a
      t38=dp1p/z2a
      t39=(dp0p-ti*dp1p)/zc
c
      do 475 l5=1,l4
      t1=W(l5)-r(Nucd(k),l5)
      dpd=t1*dp1p
      dp1d=t1*dp2p
      pip0d=t1*pipkpl
      pip1d=t1*pi2pkpl
      pjp0d=t1*pjpkpl
      pjp1d=t1*pj2pkpl
      d0d=t1*dspl
      d1d=t1*ds2pl
      dp0pm=t1*d1p
      dp1pm=t1*d2p
      
c
      if (l1.eq.l5) then
       dpd=dpd+t14
       dp1d=dp1d+t14b
       pip0d=pip0d+t20
       pip1d=pip1d+t20b
       d0d=d0d+t21
       d1d=d1d+t21b
       dp0pm=dp0pm+t26
       dp1pm=dp1pm+t26b
      endif
c
      if (l2.eq.l5) then
       dpd=dpd+t15
       dp1d=dp1d+t15b
       pjp0d=pjp0d+t20
       pjp1d=pjp1d+t20b
       d0d=d0d+t23
       d1d=d1d+t23b
       dp0pm=dp0pm+t27
       dp1pm=dp1pm+t27b
      endif
c
      if (l3.eq.l5) then
       dpd=dpd+t16
       dp1d=dp1d+t16b
       pjp0d=pjp0d+t21
       pjp1d=pjp1d+t21b
       pip0d=pip0d+t23
       pip1d=pip1d+t23b
       dp0pm=dp0pm+t28
       dp1pm=dp1pm+t28b
      endif
c
      f2=1.D0
      if (l4.eq.l5) then
       dpd=dpd+t17
       dp1d=dp1d+t17b
       pjp0d=pjp0d+t22
       pjp1d=pjp1d+t22b
       pip0d=pip0d+t25
       pip1d=pip1d+t25b
       d0d=d0d+t24
       d1d=d1d+t24b
       f2=sq3
      endif
c
      t30=(pjp0d-roz*pjp1d)/z2
      t31=pjp1d/z2a
      t32=(pip0d-roz*pip1d)/z2
      t33=pip1d/z2a
      t34=(d0d-roz*d1d)/z2
      t35=d1d/z2a
      t36=dp1pm/z2a
      t37=(dp0pm-ti*dp1pm)/zc
c     
      l12=Ll(l1)+l2
      ii=i+l12-1
      jj=j+l3-1
      l45=l4*(l4-1)/2+l5
      kk=k+l45-1
c
      cc=ccoef/(f1*f2)
      term=dpd*cc
c
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c 
c gradients
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
      do 475 l6=1,3
       t1=Q(l6)-r(Nuc(j),l6)
       t2=W(l6)-Q(l6)
       t2b=W(l6)-r(Nucd(k),l6)
       tx=r(Nuc(i),l6)-r(Nuc(j),l6)
c
      ddd=t1*dpd+t2*dp1d
      dpf=t2b*dp1d
c
      if (l1.eq.l6) then
       ddd=ddd+t30
       dpf=dpf+t31
       f(Nuc(i),l6)=f(Nuc(i),l6)-(te+tee)*pjp0d
      endif
c
      if (l2.eq.l6) then
       ddd=ddd+t32
       dpf=dpf+t33
       f(Nuc(i),l6)=f(Nuc(i),l6)-(te+tee)*pip0d
      endif
c
      if (l3.eq.l6) then
       ddd=ddd+t34
       dpf=dpf+t35
       f(Nuc(j),l6)=f(Nuc(j),l6)-(te+tee)*d0d
      endif
c
      if (l4.eq.l6) then
       ddd=ddd+t36
       dpf=dpf+t37
       f(Nucd(k),l6)=f(Nucd(k),l6)-te*dp0pm
      endif
c
      if (l5.eq.l6) then
       ddd=ddd+t38
       dpf=dpf+t39
       f(Nucd(k),l6)=f(Nucd(k),l6)-te*dp0p
      endif
c
      fpd=ddd-tx*dpd
      f(Nuc(i),l6)=f(Nuc(i),l6)+a(i,ni)*(ty+tye)*fpd
      f(Nuc(j),l6)=f(Nuc(j),l6)+a(j,nj)*(ty+tye)*ddd
      f(Nucd(k),l6)=f(Nucd(k),l6)+ty*ad(k,nk)*dpf
c
 475   continue
 471   continue
 472   continue
 470   continue
c
c-------------------------------------------------------------
c
c (dd|d) and gradients
      do 480 i=ns+np+1,M,6
      do 480 j=ns+np+1,i,6
c
      ddi=d(Nuc(i),Nuc(j))
c
      do 480 ni=1,ncont(i)
      do 480 nj=1,ncont(j)
c
       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       ti=a(i,ni)/zij
       tj=a(j,nj)/zij
       alf=a(i,ni)*tj
       rexp=alf*ddi
       if (rexp.gt.rmax) goto 482
c
       Q(1)=ti*r(Nuc(i),1)+tj*r(Nuc(j),1)
       Q(2)=ti*r(Nuc(i),2)+tj*r(Nuc(j),2)
       Q(3)=ti*r(Nuc(i),3)+tj*r(Nuc(j),3)
c
       sks=pi52*exp(-rexp)/zij
c
      do 481 k=nsd+npd+1,Md,6
c
      dpc=(Q(1)-r(Nucd(k),1))**2+(Q(2)-r(Nucd(k),2))**2+
     >    (Q(3)-r(Nucd(k),3))**2
c
      do 481 nk=1,ncontd(k)
c
c
      ccoef=c(i,ni)*c(j,nj)*cd(k,nk)
      t0=ad(k,nk)+zij
      z2a=2.D0*t0
      zc=2.D0*ad(k,nk)
c
      ti=zij/t0
      tj=ad(k,nk)/t0
      W(1)=ti*Q(1)+tj*r(Nucd(k),1)
      W(2)=ti*Q(2)+tj*r(Nucd(k),2)
      W(3)=ti*Q(3)+tj*r(Nucd(k),3)
c
      roz=tj
c
      ro=roz*zij
      u=ro*dpc
      t1=ad(k,nk)*sqrt(t0)
      t2=sks/t1
      sss=t2*FUNCT(0,u)
      ss1s=t2*FUNCT(1,u)
      ss2s=t2*FUNCT(2,u)
      ss3s=t2*FUNCT(3,u)
      ss4s=t2*FUNCT(4,u)
      ss5s=t2*FUNCT(5,u)
      ss6s=t2*FUNCT(6,u)
      ss7s=t2*FUNCT(7,u)
      t3=(sss-roz*ss1s)/z2
      t4=(ss1s-roz*ss2s)/z2
      t5=(ss2s-roz*ss3s)/z2
      t6=(ss3s-roz*ss4s)/z2
      t6b=(ss4s-roz*ss5s)/z2
      t6c=(ss5s-roz*ss6s)/z2
c
      do 485 l1=1,3
      t1=Q(l1)-r(Nuc(i),l1)
      t2=W(l1)-Q(l1)
      ps=t1*sss+t2*ss1s
      p1s=t1*ss1s+t2*ss2s
      p2s=t1*ss2s+t2*ss3s
      p3s=t1*ss3s+t2*ss4s
      p4s=t1*ss4s+t2*ss5s
      p5s=t1*ss5s+t2*ss6s
      p6s=t1*ss6s+t2*ss7s
c
      t7=(ps-roz*p1s)/z2
      t8=(p1s-roz*p2s)/z2
      t9=(p2s-roz*p3s)/z2
      t10=(p3s-roz*p4s)/z2
      t10b=(p4s-roz*p5s)/z2
      y16=p2s/z2a
      y16b=p3s/z2a
c
      do 485 l2=1,l1
      t1=Q(l2)-r(Nuc(i),l2)
      t2=W(l2)-Q(l2)
      pjs=t1*sss+t2*ss1s
      pj1s=t1*ss1s+t2*ss2s
      pj2s=t1*ss2s+t2*ss3s
      pj3s=t1*ss3s+t2*ss4s
      pj4s=t1*ss4s+t2*ss5s
      pj5s=t1*ss5s+t2*ss6s
      ds=t1*ps+t2*p1s
      d1s=t1*p1s+t2*p2s
      d2s=t1*p2s+t2*p3s
      d3s=t1*p3s+t2*p4s
      d4s=t1*p4s+t2*p5s
      d5s=t1*p5s+t2*p6s
c
      t11=(pjs-roz*pj1s)/z2
      t12=(pj1s-roz*pj2s)/z2
      t13=(pj2s-roz*pj3s)/z2
      t14=(pj3s-roz*pj4s)/z2
      t14b=(pj4s-roz*pj5s)/z2
      y19=pj2s/z2a
      y19b=pj3s/z2a
c
      f1=1.D0
      if (l1.eq.l2) then
       ds=ds+t3
       d1s=d1s+t4
       d2s=d2s+t5
       d3s=d3s+t6
       d4s=d4s+t6b
       d5s=d5s+t6c
       f1=sq3
      endif
c
      t16=(ds-roz*d1s)/z2
      t17=(d1s-roz*d2s)/z2
      t18=(d2s-roz*d3s)/z2
      t18b=(d3s-roz*d4s)/z2
      t22a=d2s/z2a
      t22c=d3s/z2a
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
      do 485 l3=1,lij
      t1=Q(l3)-r(Nuc(j),l3)
      t2=W(l3)-Q(l3)
      dpk=t1*ds+t2*d1s
      d1pk=t1*d1s+t2*d2s
      d2pk=t1*d2s+t2*d3s
      d3pk=t1*d3s+t2*d4s
      d4pk=t1*d4s+t2*d5s
c
      pjpk=t1*pjs+t2*pj1s
      pj1pk=t1*pj1s+t2*pj2s
      pj2pk=t1*pj2s+t2*pj3s
      pj3pk=t1*pj3s+t2*pj4s
      pj4pk=t1*pj4s+t2*pj5s
      pipk=t1*ps+t2*p1s
      pi1pk=t1*p1s+t2*p2s
      pi2pk=t1*p2s+t2*p3s
      pi3pk=t1*p3s+t2*p4s
      pi4pk=t1*p4s+t2*p5s
      spk=t1*sss+t2*ss1s
      s1pk=t1*ss1s+t2*ss2s
      s2pk=t1*ss2s+t2*ss3s
      s3pk=t1*ss3s+t2*ss4s
      s4pk=t1*ss4s+t2*ss5s
c
      t15p=(spk-roz*s1pk)/z2
      t15a=(s1pk-roz*s2pk)/z2
      t15=(s2pk-roz*s3pk)/z2
      t15b=(s3pk-roz*s4pk)/z2
      y18=s2pk/z2a
      y18b=s3pk/z2a
c
      if (l1.eq.l3) then
       dpk=dpk+t11
       d1pk=d1pk+t12
       d2pk=d2pk+t13
       d3pk=d3pk+t14
       d4pk=d4pk+t14b
       pipk=pipk+t3
       pi1pk=pi1pk+t4
       pi2pk=pi2pk+t5
       pi3pk=pi3pk+t6
       pi4pk=pi4pk+t6b
      endif
c
      if (l2.eq.l3) then
       dpk=dpk+t7
       d1pk=d1pk+t8
       d2pk=d2pk+t9
       d3pk=d3pk+t10
       d4pk=d4pk+t10b
       pjpk=pjpk+t3
       pj1pk=pj1pk+t4
       pj2pk=pj2pk+t5
       pj3pk=pj3pk+t6
       pj4pk=pj4pk+t6b
      endif
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      t20=pj2pk/z2a
      t20b=pj3pk/z2a
      t21=pi2pk/z2a
      t21b=pi3pk/z2a
      t22p=d1pk/z2a
      t22=d2pk/z2a
      t22b=d3pk/z2a
c
      t24=(pjpk-roz*pj1pk)/z2
      t25=(pj1pk-roz*pj2pk)/z2
      t26=(pj2pk-roz*pj3pk)/z2
      t26b=(pj3pk-roz*pj4pk)/z2
      t27=(pipk-roz*pi1pk)/z2
      t28=(pi1pk-roz*pi2pk)/z2
      t29=(pi2pk-roz*pi3pk)/z2
      t29b=(pi3pk-roz*pi4pk)/z2
c
      y10=t22p
      y10b=t22
      y14=(dpk-ti*d1pk)/zc
      y14b=(d1pk-ti*d2pk)/zc
c
      do 485 l4=1,lk
      t1=Q(l4)-r(Nuc(j),l4)
      t2=W(l4)-Q(l4)
      dd=t1*dpk+t2*d1pk
      d1d=t1*d1pk+t2*d2pk
      d2d=t1*d2pk+t2*d3pk
      d3d=t1*d3pk+t2*d4pk
c
      pj0dkl=t1*pjpk+t2*pj1pk
      pj1dkl=t1*pj1pk+t2*pj2pk
      pj2dkl=t1*pj2pk+t2*pj3pk
      pj3dkl=t1*pj3pk+t2*pj4pk
      pi0dkl=t1*pipk+t2*pi1pk
      pi1dkl=t1*pi1pk+t2*pi2pk
      pi2dkl=t1*pi2pk+t2*pi3pk
      pi3dkl=t1*pi3pk+t2*pi4pk
      d0pl=t1*ds+t2*d1s
      d1pl=t1*d1s+t2*d2s
      d2pl=t1*d2s+t2*d3s
      d3pl=t1*d3s+t2*d4s
c
      s2pl=t1*ss2s+t2*ss3s
      s3pl=t1*ss3s+t2*ss4s
      y17=s2pl/z2a
      y17b=s3pl/z2a
c
      s2dkl=t1*s2pk+t2*s3pk
      s3dkl=t1*s3pk+t2*s4pk
      pj2pl=t1*pj2s+t2*pj3s
      pj3pl=t1*pj3s+t2*pj4s
      pi2pl=t1*p2s+t2*p3s
      pi3pl=t1*p3s+t2*p4s
c
      if (l1.eq.l4) then
       dd=dd+t24
       d1d=d1d+t25
       d2d=d2d+t26
       d3d=d3d+t26b
       pi0dkl=pi0dkl+t15p
       pi1dkl=pi1dkl+t15a
       pi2dkl=pi2dkl+t15
       pi3dkl=pi3dkl+t15b
       d0pl=d0pl+t11
       d1pl=d1pl+t12
       d2pl=d2pl+t13
       d3pl=d3pl+t14
       pi2pl=pi2pl+t5
       pi3pl=pi3pl+t6
      endif
c
      if (l2.eq.l4) then
       dd=dd+t27
       d1d=d1d+t28
       d2d=d2d+t29
       d3d=d3d+t29b
       pj0dkl=pj0dkl+t15p
       pj1dkl=pj1dkl+t15a
       pj2dkl=pj2dkl+t15
       pj3dkl=pj3dkl+t15b
       d0pl=d0pl+t7
       d1pl=d1pl+t8
       d2pl=d2pl+t9
       d3pl=d3pl+t10
       pj2pl=pj2pl+t5
       pj3pl=pj3pl+t6
      endif
c
      f2=1.D0
      if (l3.eq.l4) then
       s2dkl=s2dkl+t5
       s3dkl=s3dkl+t6
       dd=dd+t16
       d1d=d1d+t17
       d2d=d2d+t18
       d3d=d3d+t18b
       pj0dkl=pj0dkl+t11
       pj2dkl=pj2dkl+t13
       pj1dkl=pj1dkl+t12
       pj3dkl=pj3dkl+t14
       pi2dkl=pi2dkl+t9
       pi1dkl=pi1dkl+t8
       pi3dkl=pi3dkl+t10
       pi0dkl=pi0dkl+t7
       f2=sq3
      endif
c
      t30a=pj1dkl/z2a
      t30=pj2dkl/z2a
      t30b=pj3dkl/z2a
      t40a=pi1dkl/z2a
      t40=pi2dkl/z2a
      t40b=pi3dkl/z2a
      t50=s2dkl/z2a
      t50b=s3dkl/z2a
      t60=pj2pl/z2a
      t60b=pj3pl/z2a
      t70=pi2pl/z2a
      t70b=pi3pl/z2a
      t80a=d1pl/z2a
      t80=d2pl/z2a
      t80b=d3pl/z2a
      t23=(dd-ti*d1d)/zc
      t23b=(d1d-ti*d2d)/zc
c
      y1=t30a
      y1b=t30
      y5=t40a
      y5b=t40
      y8=t80a
      y8b=t80
      y12=(pi0dkl-ti*pi1dkl)/zc
      y12b=(pi1dkl-ti*pi2dkl)/zc
      y13=(d0pl-ti*d1pl)/zc
      y13b=(d1pl-ti*d2pl)/zc
      y15=(pj0dkl-ti*pj1dkl)/zc
      y15b=(pj1dkl-ti*pj2dkl)/zc
c
      do 485 l5=1,3
c
      t1=W(l5)-r(Nucd(k),l5)
      dd0p=t1*d1d
      ddp=t1*d2d
      dd2p=t1*d3d
      pjdklp=t1*pj2dkl
      pj2dklp=t1*pj3dkl
      pidklp=t1*pi2dkl
      pi2dklp=t1*pi3dkl
      dijplp=t1*d2pl
      dij2plp=t1*d3pl
      dijpkp=t1*d2pk
      dij2pkp=t1*d3pk
c
      s1dpm=t1*s2dkl
      s2dpm=t1*s3dkl
      pj1plpm=t1*pj2pl
      pj2plpm=t1*pj3pl
      pj1pkpm=t1*pj2pk
      pj2pkpm=t1*pj3pk
      pi1plpm=t1*pi2pl
      pi2plpm=t1*pi3pl
      pi1pkpm=t1*pi2pk
      pi2pkpm=t1*pi3pk
      d1spm=t1*d2s
      d2spm=t1*d3s
c
      if (l1.eq.l5) then
       dd0p=dd0p+t30a
       ddp=ddp+t30
       dd2p=dd2p+t30b
       pidklp=pidklp+t50
       pi2dklp=pi2dklp+t50b
       dijplp=dijplp+t60
       dij2plp=dij2plp+t60b
       dijpkp=dijpkp+t20
       dij2pkp=dij2pkp+t20b
       pi1plpm=pi1plpm+y17
       pi2plpm=pi2plpm+y17b
       pi1pkpm=pi1pkpm+y18
       pi2pkpm=pi2pkpm+y18b
       d1spm=d1spm+y19
       d2spm=d2spm+y19b
      endif
c
      if (l2.eq.l5) then
       dd0p=dd0p+t40a
       ddp=ddp+t40
       dd2p=dd2p+t40b
       pjdklp=pjdklp+t50
       pj2dklp=pj2dklp+t50b
       dijplp=dijplp+t70
       dij2plp=dij2plp+t70b
       dijpkp=dijpkp+t21
       dij2pkp=dij2pkp+t21b
       pj1plpm=pj1plpm+y17
       pj2plpm=pj2plpm+y17b
       pj1pkpm=pj1pkpm+y18
       pj2pkpm=pj2pkpm+y18b
       d1spm=d1spm+y16
       d2spm=d2spm+y16b
      endif
c
      if (l3.eq.l5) then
       dd0p=dd0p+t80a
       ddp=ddp+t80
       dd2p=dd2p+t80b
       pjdklp=pjdklp+t60
       pj2dklp=pj2dklp+t60b
       pidklp=pidklp+t70
       pi2dklp=pi2dklp+t70b
       dijpkp=dijpkp+t22a
       dij2pkp=dij2pkp+t22c
       s1dpm=s1dpm+y17
       s2dpm=s2dpm+y17b
       pj1pkpm=pj1pkpm+y19
       pj2pkpm=pj2pkpm+y19b
       pi1pkpm=pi1pkpm+y16
       pi2pkpm=pi2pkpm+y16b
      endif
c
      if (l4.eq.l5) then
       dd0p=dd0p+t22p
       ddp=ddp+t22
       dd2p=dd2p+t22b
       pjdklp=pjdklp+t20
       pj2dklp=pj2dklp+t20b
       pidklp=pidklp+t21
       pi2dklp=pi2dklp+t21b
       dijplp=dijplp+t22a
       dij2plp=dij2plp+t22c
       s1dpm=s1dpm+y18
       s2dpm=s2dpm+y18b
       pj1plpm=pj1plpm+y19
       pj2plpm=pj2plpm+y19b
       pi1plpm=pi1plpm+y16
       pi2plpm=pi2plpm+y16b
      endif
c
      t31=pjdklp/z2a
      t31b=pj2dklp/z2a
      t41=pidklp/z2a
      t41b=pi2dklp/z2a
      t51=dijplp/z2a
      t51b=dij2plp/z2a
      t61=dijpkp/z2a
      t61b=dij2pkp/z2a
c
      y30=ddp/z2a
      y31=(dd0p-ti*ddp)/zc
c
      y2=s1dpm/z2a
      y2b=s2dpm/z2a
      y3=pj1plpm/z2a
      y3b=pj2plpm/z2a
      y4=pj1pkpm/z2a
      y4b=pj2pkpm/z2a
      y6=pi1plpm/z2a
      y6b=pi2plpm/z2a
      y7=pi1pkpm/z2a
      y7b=pi2pkpm/z2a
      y9=d1spm/z2a
      y9b=d2spm/z2a
c
      do 485 l6=1,l5
c
      
      t1=W(l6)-r(Nucd(k),l6)
      
      ddd=t1*ddp
      dd1d=t1*dd2p
      dd0pn=t1*d1d
      dd1pn=t1*d2d
      pj0dd=t1*pjdklp
      pj1dd=t1*pj2dklp
      pi0dd=t1*pidklp
      pi1dd=t1*pi2dklp
      d0pld=t1*dijplp
      d1pld=t1*dij2plp
      d0pkd=t1*dijpkp
      d1pkd=t1*dij2pkp
c
      if (l1.eq.l6) then
       ddd=ddd+t31
       dd1d=dd1d+t31b
       dd0pn=dd0pn+y1
       dd1pn=dd1pn+y1b
       pi0dd=pi0dd+y2
       pi1dd=pi1dd+y2b
       d0pld=d0pld+y3
       d1pld=d1pld+y3b
       d0pkd=d0pkd+y4
       d1pkd=d1pkd+y4b
      endif
c
      if (l2.eq.l6) then
       ddd=ddd+t41
       dd1d=dd1d+t41b
       dd0pn=dd0pn+y5
       dd1pn=dd1pn+y5b
       d0pld=d0pld+y6
       d1pld=d1pld+y6b
       d0pkd=d0pkd+y7
       d1pkd=d1pkd+y7b
       pj0dd=pj0dd+y2
       pj1dd=pj1dd+y2b
      endif
c
      if (l3.eq.l6) then
       ddd=ddd+t51
       dd1d=dd1d+t51b
       dd0pn=dd0pn+y8
       dd1pn=dd1pn+y8b
       pj0dd=pj0dd+y3
       pj1dd=pj1dd+y3b
       d0pkd=d0pkd+y9
       d1pkd=d1pkd+y9b
       pi0dd=pi0dd+y6
       pi1dd=pi1dd+y6b
      endif
c
      if (l4.eq.l6) then
       ddd=ddd+t61
       dd1d=dd1d+t61b
       pi0dd=pi0dd+y7
       pi1dd=pi1dd+y7b
       d0pld=d0pld+y9
       d1pld=d1pld+y9b
       dd0pn=dd0pn+y10
       dd1pn=dd1pn+y10b
       pj0dd=pj0dd+y4
       pj1dd=pj1dd+y4b
      endif
c
      f3=1.D0
      if (l5.eq.l6) then
       ddd=ddd+t23
       dd1d=dd1d+t23b
       pi0dd=pi0dd+y12
       pi1dd=pi1dd+y12b
       d0pld=d0pld+y13
       d1pld=d1pld+y13b
       d0pkd=d0pkd+y14
       d1pkd=d1pkd+y14b
       pj0dd=pj0dd+y15
       pj1dd=pj1dd+y15b
       f3=sq3
      endif
c
      cc=ccoef/(f1*f2*f3)
      term=ddd*cc
c
      l12=Ll(l1)+l2
      l34=Ll(l3)+l4
      l56=Ll(l5)+l6
      ii=i+l12-1
      jj=j+l34-1
      kk=k+l56-1
c
      y20=(pj0dd-roz*pj1dd)/z2
      y21=pj1dd/z2a
      y22=(pi0dd-roz*pi1dd)/z2
      y23=pi1dd/z2a
      y24=(d0pld-roz*d1pld)/z2
      y25=d1pld/z2a
      y26=(d0pkd-roz*d1pkd)/z2
      y27=d1pkd/z2a
      y28=dd1pn/z2a
      y29=(dd0pn-ti*dd1pn)/zc
c
      kn=ii+Jx(jj)
      RMM(M5+kn-1)=RMM(M5+kn-1)+af2(kk)*term
c gradients
      t1=cc*RMM(kn)
      te=t1*af(kk)
      tee=t1*B(kk,2)
      ty=2.D0*te
      tye=2.D0*tee
c
      do 485 l7=1,3
       t1=Q(l7)-r(Nuc(i),l7)
       t2=W(l7)-Q(l7)
       t2b=W(l7)-r(Nucd(k),l7)
       tx=r(Nuc(i),l7)-r(Nuc(j),l7)
c
       fdd=t1*ddd+t2*dd1d
       ddf=t2b*dd1d
c
       if (l1.eq.l7) then
        fdd=fdd+y20
        ddf=ddf+y21
        f(Nuc(i),l7)=f(Nuc(i),l7)-(te+tee)*pj0dd
       endif
c
       if (l2.eq.l7) then
        fdd=fdd+y22
        ddf=ddf+y23
        f(Nuc(i),l7)=f(Nuc(i),l7)-(te+tee)*pi0dd
       endif
c
       if (l3.eq.l7) then
        fdd=fdd+y24
        ddf=ddf+y25
        f(Nuc(j),l7)=f(Nuc(j),l7)-(te+tee)*d0pld
       endif
c
       if (l4.eq.l7) then
        fdd=fdd+y26
        ddf=ddf+y27
        f(Nuc(j),l7)=f(Nuc(j),l7)-(te+tee)*d0pkd
       endif
c
       if (l5.eq.l7) then
        fdd=fdd+y28
        ddf=ddf+y29
        f(Nucd(k),l7)=f(Nucd(k),l7)-te*dd0pn
       endif
c
       if (l6.eq.l7) then
        fdd=fdd+y30
        ddf=ddf+y31
        f(Nucd(k),l7)=f(Nucd(k),l7)-te*dd0p
       endif
c
       dfd=fdd+tx*ddd
c
       f(Nuc(i),l7)=f(Nuc(i),l7)+a(i,ni)*(ty+tye)*fdd
       f(Nuc(j),l7)=f(Nuc(j),l7)+a(j,nj)*(ty+tye)*dfd
       f(Nucd(k),l7)=f(Nucd(k),l7)+ad(k,nk)*ty*ddf
 485  continue
 481  continue
 482  continue
 480  continue
c
c
c-------------------------------------------------------------
c debuggings and checks
c

c
c     
c     write(*,*) 'Forces for 2 electron terms'
c     write(*,*) Ex,Ea-Eb/2.D0
      
c      do i=1,natom
c       write(22,*) i,f(i,1),f(i,2),f(i,3)
c      enddo
c
c
      endif
      call g2g_timer_stop('CoulG')
      call g2g_timer_sum_stop('Coulomb gradients')
      return
      end subroutine
      end module subm_int3G
