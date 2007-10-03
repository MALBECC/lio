c function calculating density functionals
c standard, for local density functionals, recalculates
c everything on every iteration
c
      SUBROUTINE DNSOP(densa,densb,F,Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,
     >               M,M18,NCOa,NCOb,RMM)
    
c
      implicit real*8 (a-h,o-z)
      logical NORM
      INCLUDE 'param'
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng)
      dimension r(nt,3),nshell(0:3),Xi(3)
      dimension ds(ntq),F(ng),W(ng),RMM(*)
      dimension indx(ng),FF(ng)
c
      common /Ll/ Ll(3)
      common /Nc/ Ndens
c
c now we should evaluate all same loops as the ones used for
c 1 electron matrix elements, but doing only products
c then, the particular density functional wanted is calculated
c
      fc=1.D0
      if (NORM) then
       fc=1.D0/sqrt(3.D0)
      endif
c
      DENSA=0.0D0
      DENSB=0.0D0
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
      M18b=M18+M*NCOa
c
c basis functions evaluated at r are calculated
c
      do 1 i=1,M
        W(i)=0.D0
 1      F(i)=0.D0
c
c --- s  case -------
      do 10 i=1,ns
c
      di=ds(Nuc(i))
c
      do 15 ni=1,ncont(i)
c
      rexp=a(i,ni)*di
c
      if (rexp.gt.30.D0) go to 16
      t=exp(-rexp)
      F(i)=F(i)+t*c(i,ni)
  16  continue
  15  continue
c
  10  continue
c
c--- p  case -------------
      do 20 i=ns+1,ns+np,3
c
      di=ds(Nuc(i))
c
      do 20 ni=1,ncont(i)

      rexp=a(i,ni)*di
      if (rexp.gt.30.D0) goto 21
      t=exp(-rexp)*c(i,ni)
      do 25 l1=1,3
c
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      ii=i+l1-1
c
      F(ii)=F(ii)+term
  25  continue
c
  21  continue
  20  continue
c
c-- d case  ------------
      do 40 i=ns+np+1,M,6
c
      di=ds(Nuc(i))
c
      do 40 ni=1,ncont(i)
c
      rexp=a(i,ni)*di
c
      if (rexp.gt.30.) goto 41
      t=exp(-rexp)*c(i,ni)
      do 45 l1=1,3
c
      t1=xi(l1)-r(Nuc(i),l1)
      do 45 l2=1,l1
      t2=xi(l2)-r(Nuc(i),l2)
      if (l1.eq.l2) then
       t2=t2*fc
      endif
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
  45  continue
c
  41  continue
  40  continue
c
c construction of sparse F vector   -----------------------
c slower (at least for moderately small molecules)
c
        k=0
      do i=1,M
       if (F(i).ne.0.0D0) then
        k=k+1
        indx(k)=i
       endif
       enddo
c
c ------ DENSA calculation , M18 pointer for Alpha spin coefficients
c NCOa : # of electrons of spin alpha
      kk=M18-1
c
      if (k.lt.M/2) then

      do jj=1,NCOa
       W(jj)=0.0D0     
       do ii=1,k
        ik=indx(ii)
        W(jj)=W(jj)+F(ik)*RMM(kk+ik)
       enddo
       kk=kk+M
      enddo
c
      else
c
      do 153 jj=1,NCOa
       W(jj)=0.0D0
      do 154 ii=1,M
       kk=kk+1
       W(jj)=W(jj)+RMM(kk)*F(ii)
 154    continue
 153    continue
       endif
c
      do 61 ii=1,NCOa
 61    DENSA=DENSA+W(ii)**2
c
c ------ DENSB calculation , Mb pointer for Beta spin coefficients
c NCOb : # of electrons of spin Beta
c
      kk=M18b-1
c
      if (k.lt.M/2) then

      do jj=1,NCOb
       W(jj)=0.0D0
       do ii=1,k
        ik=indx(ii)
        W(jj)=W(jj)+F(ik)*RMM(kk+ik)
       enddo
       kk=kk+M
      enddo
c
      else
c
      do 253 jj=1,NCOb
       W(jj)=0.0D0
      do 254 ii=1,M
       kk=kk+1
       W(jj)=W(jj)+RMM(kk)*F(ii)
 254    continue
 253    continue
       endif
c
      do 71 ii=1,NCOb
 71    DENSB=DENSB+W(ii)**2
c
c------------------------------------
c
c     do i=1,M
c      write(88,*) i,F(i),W(i)
c     enddo
c
c     if (DENSA.ne.0) then
c     write(*,*) DENSA,DENSB
c     endif
      return
c
      end
