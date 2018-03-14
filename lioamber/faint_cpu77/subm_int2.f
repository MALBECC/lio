c-------------------------------------------------------------------
c Integrals subroutine -Second part
c 2 e integrals, 2 index : density fitting functions
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
c Input :  density basis 
c Output: G matrix  
c G matrix should be inverted, 
c later on, for evaluating  Coulomb terms
c-----------------------------------------------------------------
      module subm_int2; contains
      subroutine int2()
      use liotemp   , only: FUNCT
      use garcha_mod, only: RMM, XX, ngd, M, Md, ad, nucd, ncontd,
     >                      r, d, cd, SVD, NORM, pi5, nshelld

!
!     implicit real*8 (a-h,o-z)
      implicit none
      real*8, dimension(:), allocatable :: dgelss_temp
      real*8, dimension(Md)             :: inv_work
      integer                           :: XXX(8*Md)
!
!     Aux . things
      real*8, dimension(:), allocatable :: trabajo
      real*8                            :: Q(3), aux(ngd), Det(2)
!
!     Ex Implicits
      real*8  :: t0, t1, t2, t3, t4, t5
      real*8  :: ss, s0s, s1s, s2s, s3s, s4s
      real*8  :: ps, pjs, pjp, pj2s, pis, pip, pi2s, pi3s
      real*8  :: alf, cc, ccoef, d1s, d2s, dd, dp, ds
      real*8  :: roz, rcond, f1, f2, sq3
      real*8  :: u, tmp, tn, tj, ti, t6, z2, za, zc, zij

      integer :: igpu, info
      integer :: i, j, ii, jj, iii, jjj, k, kk, kkk
      integer :: nsd, npd, ndd, ni, nj
      integer :: lll, l12, l34, l1, l2, l3, l4, lij, lk
      integer :: MM, MMp, MMd, Md2, Md3, Md5
      integer :: M1, M2, M3, M5, M7, M9, M10, M11, M12, M13, M15

c
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
c      do 50 i=1,natom
c      do 50 j=1,natom
c       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
c     >        (r(i,3)-r(j,3))**2
c 50   continue
c
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md
      M2=2*M
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2

c
c pointers
c
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c now F
      M13=M11+MM
c auxiliar things
      M15=M13+M
c
c
c end ------------------------------------------------
      do 1 k=1,MMd
 1      RMM(M7+k-1)=0.D0
c
c--- 2 index electron repulsion for density basis set
c
c first loop (s|s) case -------------------------------------------
c
      do 200 i=1,nsd
      do 200 j=1,i
c
      dd=d(Nucd(i),Nucd(j))
c
      do 200 ni=1,ncontd(i)
      do 200 nj=1,ncontd(j)
c
c (0|0) calculation
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      
      u=alf*dd
      ccoef=cd(i,ni)*cd(j,nj)
c
      s0s=pi5/t1*FUNCT(0,u)
c
      k=i+((Md2-j)*(j-1))/2
      RMM(M7+k-1)=RMM(M7+k-1)+ccoef*s0s
c
 200  continue
c
c------------------------------------------------------------------
c
c second loop  (p|s) case
c
c
      do 300 i=nsd+1,nsd+npd,3
      do 300 j=1,nsd
c
      dd=d(Nucd(i),Nucd(j))
c
      do 300 ni=1,ncontd(i)
      do 300 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
c
      u=alf*dd
      s1s=t2*FUNCT(1,u)
c
c
      ccoef=cd(i,ni)*cd(j,nj)
c
c l2: different p in the p shell ( x,y,z respectively)
c
      do 305 l1=1,3
        t1=Q(l1)-r(Nucd(i),l1)
        tn=t1*s1s
c
        iii=i+l1-1
c ii index , taking into account different components of the shell
c
        k=iii+((Md2-j)*(j-1))/2
        RMM(M7+k-1)=RMM(M7+k-1)+tn*ccoef
 305   continue
 300   continue
c
c------------------------------------------------------------
c
c (p|p) case
c
      do 400 i=nsd+1,nsd+npd,3
      do 400 j=nsd+1,i,3
c
      dd=d(Nucd(i),Nucd(j))
c
      do 400 ni=1,ncontd(i)
      do 400 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
c
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)

      ccoef=cd(i,ni)*cd(j,nj)
c
      do 405 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       ps=t1*s2s
c
       do 405 l2=1,3
c
       t1=Q(l2)-r(Nucd(j),l2)
       tn=t1*ps
c
       if (l1.eq.l2) then
        tn=tn+s1s/z2
       endif
c
       iii=i+l1-1
       jj=j+l2-1
c
c      this to convert to 1 dimensional array, in diagonal case
c      we calculate more things than necessary . They should be
c      eliminated
       if(iii.ge.jj) then
       k=iii+((Md2-jj)*(jj-1))/2
       RMM(M7+k-1)=RMM(M7+k-1)+tn*ccoef
       endif
 405  continue
c
 400  continue
c-------------------------------------------------------------------
c (d|s) case
      do 500 i=nsd+npd+1,Md,6
      do 500 j=1,nsd
c
      dd=d(Nucd(i),Nucd(j))
c
      do 500 ni=1,ncontd(i)
      do 500 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      roz=ad(j,nj)/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
c
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)

      ccoef=cd(i,ni)*cd(j,nj)
c
      do 505 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       ps=t1*s2s
c
       do 505 l2=1,l1
c
       t1=Q(l2)-r(Nucd(i),l2)
       tn=t1*ps
c
       f1=1.
       if (l1.eq.l2) then
        tn=tn+(s0s-roz*s1s)/(2.D0*ad(i,ni))
        f1=sq3
       endif
c
       l12=l1*(l1-1)/2+l2
       iii=i+l12-1
c
       cc=ccoef/f1
       k=iii+((Md2-j)*(j-1))/2
       RMM(M7+k-1)=RMM(M7+k-1)+tn*cc
 505  continue
c
 500  continue
c-------------------------------------------------------------------
c (d|p) case
      do 600 i=nsd+npd+1,Md,6
      do 600 j=nsd+1,nsd+npd,3
c
      dd=d(Nucd(i),Nucd(j))
c
      do 600 ni=1,ncontd(i)
      do 600 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      Q(1)=(ad(i,ni)*r(Nucd(i),1)+ad(j,nj)*r(Nucd(j),1))/zij
      Q(2)=(ad(i,ni)*r(Nucd(i),2)+ad(j,nj)*r(Nucd(j),2))/zij
      Q(3)=(ad(i,ni)*r(Nucd(i),3)+ad(j,nj)*r(Nucd(j),3))/zij
c
      u=alf*dd
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
c
      ccoef=cd(i,ni)*cd(j,nj)
c
      do 605 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pis=t1*s2s
       pi2s=t1*s3s
c
       do 605 l2=1,l1
c
       t1=Q(l2)-r(Nucd(i),l2)
       pjs=t1*s2s
c
       ds=t1*pi2s
c
       f1=1.D0
       if (l1.eq.l2) then
        f1=sq3
        ds=ds+(s1s-alf*s2s/ad(i,ni))/(2.D0*ad(i,ni))
       endif

c index of p
c
       do 605 l3=1,3
c
       t0=Q(l3)-r(Nucd(j),l3)
       tn=t0*ds
c
       if (l1.eq.l3) then
        tn=tn+pjs/z2
       endif
c
       if (l2.eq.l3) then
        tn=tn+pis/z2
       endif
c
c
       l12=l1*(l1-1)/2+l2
       iii=i+l12-1
       jj=j+l3-1
c
       cc=ccoef/f1
c
       k=iii+((Md2-jj)*(jj-1))/2
       RMM(M7+k-1)=RMM(M7+k-1)+tn*cc
 605  continue
c
 600  continue
c
c-------------------------------------------------------------------
c (d|d) case
      do 700 i=nsd+npd+1,Md,6
      do 700 j=nsd+npd+1,i,6
c
      dd=d(Nucd(i),Nucd(j))
c
      do 700 ni=1,ncontd(i)
      do 700 nj=1,ncontd(j)
c
      zij=ad(i,ni)+ad(j,nj)
      z2=2.D0*zij
      za=2.D0*ad(i,ni)
      zc=2.D0*ad(j,nj)
      t0=ad(i,ni)*ad(j,nj)
      alf=t0/zij
      t1=sqrt(zij)*t0
      t2=pi5/t1
c
      ti=ad(i,ni)/zij
      tj=ad(j,nj)/zij
      Q(1)=ti*r(Nucd(i),1)+tj*r(Nucd(j),1)
      Q(2)=ti*r(Nucd(i),2)+tj*r(Nucd(j),2)
      Q(3)=ti*r(Nucd(i),3)+tj*r(Nucd(j),3)
c
      u=alf*dd
      s0s=t2*FUNCT(0,u)
      s1s=t2*FUNCT(1,u)
      s2s=t2*FUNCT(2,u)
      s3s=t2*FUNCT(3,u)
      s4s=t2*FUNCT(4,u)
c
      t3=(s0s-tj*s1s)/za
      t4=(s1s-tj*s2s)/za
      t5=(s2s-tj*s3s)/za
      ccoef=cd(i,ni)*cd(j,nj)
c
      do 705 l1=1,3
       t1=Q(l1)-r(Nucd(i),l1)
       pis=t1*s2s
       pi2s=t1*s3s
       pi3s=t1*s4s
c
       do 705 l2=1,l1
c
       t1=Q(l2)-r(Nucd(i),l2)
       pjs=t1*s2s
       pj2s=t1*s3s
c
       ds=t1*pis
       d1s=t1*pi2s
       d2s=t1*pi3s
c
       f1=1.D0
       if (l1.eq.l2) then
        ds=ds+t3
        d1s=d1s+t4
        d2s=d2s+t5
        f1=sq3
       endif

c
       t6=(ds-ti*d1s)/zc
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
       do 705 l3=1,lij
c
       t0=Q(l3)-r(Nucd(j),l3)
       dp=t0*d2s
       pip=t0*pi2s
       pjp=t0*pj2s
c
       if (l1.eq.l3) then
        dp=dp+pj2s/z2
        pip=pip+s2s/z2
       endif
c
       if (l2.eq.l3) then
        dp=dp+pi2s/z2
        pjp=pjp+s2s/z2
       endif
c
c
      lk=l3
      if (i.eq.j) then
       lll=l1*(l1-1)/2-l3*(l3-1)/2+l2
       lk=min(l3,lll)
      endif
       do 705 l4=1,lk
c
       t0=Q(l4)-r(Nucd(j),l4)
       tn=t0*dp
c
       if (l1.eq.l4) then
        tn=tn+pjp/z2
       endif
c
       if (l2.eq.l4) then
        tn=tn+pip/z2
       endif
c
       f2=1.D0
       if (l3.eq.l4) then
        tn=tn+t6
        f2=sq3
       endif
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       iii=i+l12-1
       jj=j+l34-1
c
       cc=ccoef/(f1*f2)
c
c      this to convert to 1 dimensional array, in diagonal case
       k=iii+(Md2-jj)*(jj-1)/2
c
       RMM(M7+k-1)=RMM(M7+k-1)+tn*cc
 705  continue
c
 700  continue
c
c
c
       do 216 i=1,Md
       do 216 j=1,Md
c
        if(i.ge.j) then
         k=i+(Md*2-j)*(j-1)/2
         else
         k=j+(Md*2-i)*(i-1)/2
        endif
        XX(i,j)=RMM(M7+k-1)
 216   continue
c
c
       kk=0
      do 112 j=1,Md
      do 112 i=j,Md
       kk=kk+1
       tmp=XX(i,j)
       RMM(M7+kk-1)=tmp
 112  continue
c
      MMp=Md*(Md+1)/2
      do 199 k=1,MMp
 199   RMM(M9+k-1)=0.0D0
c
c     M10=M9+Md
      M10=M9+MMd+MM+1
      M12=M10+Md
      Md3=3*Md
c ESSL OPTION ------------------------------
#ifdef essl
      CALL DGESVF(10,XX,Md,RMM(M9),Md,1,RMM(M10),
     >             Md,Md,RMM(M12),Md3)
c
       ss=RMM(M10)/RMM(M10+Md-1)
c
#endif
c
c LAPACK OPTION ------------------------------
c
#ifdef pack

c
       do i=1,Md
        aux(i)=0.0D0
       enddo
       Md5=5*Md
      rcond=1.0D-07

c
c CH - why call dgelss here? We only want the singular values - couldn't just
c something like dgesvd be called without calculating singular vectors?
c
c      call dgelss(Md,Md,1,XX,Md,aux,Md,RMM(M9),rcond,irank,RMM(M10),
c     >            -1,info)
c      Md5=RMM(M10)
c      allocate(dgelss_temp(Md5))
c      call dgelss(Md,Md,1,XX,Md,aux,Md,RMM(M9),rcond,irank,dgelss_temp,
c     >            Md5,info)
c      deallocate(dgelss_temp)

      call g2g_timer_sum_start('G condition')
#ifdef magma
      call magmaf_dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1,
     >            RMM(M10),-1,XXX,info)
#else
      call dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1,
     >            RMM(M10),-1,XXX,info)
#endif
      Md5=RMM(M10)
      allocate(dgelss_temp(Md5))
#ifdef magma
      call magmaf_dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1,
     >            dgelss_temp,Md5,XXX,info)
#else
      call dgesdd('N',Md,Md,XX,Md,RMM(M9),0,1,0,1,
     >            dgelss_temp,Md5,XXX,info)
#endif
      deallocate(dgelss_temp)

      ss=RMM(M9)/RMM(M9+Md-1)

c
#endif
c	write (*,*) ss, "criterio ajuste base auxiliar, Nick"
       if (ss.gt.1.D14) then
        SVD=.true.
	stop "trata de usar SVD"
       endif

      call g2g_timer_sum_stop('G condition')
c
c------------------------------
c inversion of G matrix , kept in Gm
c
      if (SVD) then
       write(*,900) ss
       call aint_query_gpu_level(igpu)
       if (igpu.eq.5) then
         write(*,*) "G IS ILL-CONDITIONED"
         write(*,*) "THE SVD AUXILIARY DENSITY FIT IS NOT SUPPORTED"
         write(*,*) "IN THE GPU VERSION OF LIO"
         stop
       endif
       
      else
c
c
c
c
c LINPACK OPTION
#ifdef pack
c
      call g2g_timer_sum_start('G invert')

       do i=1,Md
       do j=1,Md

        if(i.ge.j) then
         k=i+(Md*2-j)*(j-1)/2
         else
         k=j+(Md*2-i)*(i-1)/2
        endif
        XX(i,j)=RMM(M7+k-1)
      enddo
      enddo

c      kk=0
c      do 313 j=1,Md
c       do 313 i=1,j
c       kk=kk+1
c       kx=M7+j+(2*Md-i)*(i-1)/2-1
c       RMM(M9+kk-1)=RMM(kx)
c
c 313  continue
c

c      call dppco(RMM(M9),Md,rcond,aux,info)
      call dsytrf('U',Md,XX,Md,XXX,RMM(M10),-1,info)
      Md5=RMM(M10)
      allocate(dgelss_temp(Md5))
      call dsytrf('U',Md,XX,Md,XXX,dgelss_temp,Md5,info)
      deallocate(dgelss_temp)

c      call dpptri('U',Md,RMM(M9), info)
      call dsytri('U',Md,XX,Md,XXX,inv_work,info)

c      call dppdi(RMM(M9),Md,det,1)
c
c      kk=0
c      do 314 j=1,Md
c      do 314 i=1,j
c      kk=kk+1
c      kx=j+(2*Md-i)*(i-1)/2-1
c       RMM(M15+kx)=RMM(M9+kk-1)
c 314  continue
c
c      do 315 kk=1,MMp
c
c 315   RMM(M9+kk-1)=RMM(M15+kk-1)
c
      do i=1,Md
      do j=1,i
        k=i+(Md*2-j)*(j-1)/2
        RMM(M9+k-1) = XX(j,i)
      enddo
      enddo

      call g2g_timer_sum_stop('G invert')
#endif
c
      endif
 900  format('SWITCHING TO SVD rcond=',D10.3)
c
c-------------------------------------------------------------------
      return
      end subroutine
      end module subm_int2

