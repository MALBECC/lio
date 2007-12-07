c function calculating density functionals
c standard, for local density functionals, recalculates
c everything on every iteration
c
      SUBROUTINE DNS(Dens,F,Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,
     >               M,M18,NCO,RMM)

c
      implicit real*8 (a-h,o-z)
      logical NORM
      INCLUDE 'param'
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng)
      dimension r(nt,3),nshell(0:3),Xi(3)
      dimension ds(ntq),F(M),W(ng),RMM(*)
      dimension indx(ng)
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
      DENS=0.D0
      DENS1=0.0D0
c
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
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
c        write(*,*) 'rexp',i,ni,rexp,c(i,ni)
c        write(*,*) 'func: ', i, ' cont: ', ni, ' a: ',
c     >    a(i,ni),' c: ',c(i,ni)
c
      if (rexp.gt.30.D0) go to 16
      t=exp(-rexp)
      F(i)=F(i)+t*c(i,ni)
  16  continue
  15  continue
c      if (rexp.gt.30.D0) then
c      else
c        write(213,*) 'F(',i-1,') s:',F(i)
c      endif

  10  continue

c DEBUG DEBUG DEBUG
c      do i=1,ns
c        write(*,*) 's:',F(i)
c      enddo
      
c
c--- p  case -------------
*c va de a tres, para luego ir por esos 3
      do 20 i=ns+1,ns+np,3
c
        di=ds(Nuc(i))
c        write(214,*) 'cont:',ncont(i)
c
      do 20 ni=1,ncont(i)

        rexp=a(i,ni)*di
c        write(*,*) 'rexp',i,ni,rexp,c(i,ni)
c        write(*,*) 'func: ', i, ' cont: ', ni, ' a: ',
c     >    a(i,ni), ' c: ', c(i,ni)        
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
c      if (rexp.gt.30.D0) then
c      else      
c        write(214,*) 'F(',i-1,') p:',F(i),F(i+1),F(i+2)
c      endif
  20  continue

c DEBUG DEBUG DEBUG
c      do i=ns+1,ns+np,3
c        t=0.0
c        write(*,*) 'p:',F(i),F(i+1),F(i+2)
c      enddo      
c
c-- d case  ------------
      do 40 i=ns+np+1,M,6
c
      di=ds(Nuc(i))
c
      do 40 ni=1,ncont(i)
c
        rexp=a(i,ni)*di
c        write(*,*) 'rexp',i,ni,rexp,c(i,ni)
c
c      write(*,*) 'func: ', i, ' cont: ', ni, ' a: ',
c     >  a(i,ni), ' c: ', c(i,ni)      
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

c      if (rexp.gt.30.D0) then
c      else
c      write(215,*)'F(',i-1,')d:',F(i),F(i+1),F(i+2),F(i+3),F(i+4),F(i+5)
c      endif
      
  40  continue

c DEBUG DEBUG DEBUG
c      do i=ns+np+1,M,6
c        t=0.0
c        do l1=1,6
c          t=t + F(i + l1 - 1)
c        enddo
c        write(*,*) 'd: F(',i-1,')=',t
c      enddo
c
c now calculation of vector W : density matrix scalar F
c
c

      if (Ndens.eq.1) then
      k=0
      do 50 j=1,M
c
       if (F(j).eq.0.0D0) then
        k=k+M-j+1
        goto 50
       endif
c
      do 51 i=j,M
        k=k+1
c        write(*,*) 'Ndens1 RMM(',k-1,')=',RMM(k)
 51   W(j)=W(j)+RMM(k)*F(i)
 50   continue
c
      do 60 i=1,M
       DENS=DENS+F(i)*W(i)
  60  continue
c
      return
c    
      else
c
c      kk=M18-1
c
c     do 53 jj=1,NCO
c     do 54 ii=1,M
c      kk=kk+1
c      W(jj)=W(jj)+RMM(kk)*F(ii)
c54    continue
c53    continue
c
c
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
      kk=M18-1
c
c DEBUG DEBUG DEBUG Matias
c      if (k.lt.M/2) then
c
c      do jj=1,NCO
c       W(jj)=0.0D0     
c       do ii=1,k
c        ik=indx(ii)
c        W(jj)=W(jj)+F(ik)*RMM(kk+ik)
c       enddo
c       kk=kk+M
c      enddo
c
c      else
c
#ifndef GPU     
      do i=1,M
c        write(*,*) 'F(',i-1,')=',F(i)
      enddo
#endif
      
      do 153 jj=1,NCO
      do 154 ii=1,M
       kk=kk+1
       W(jj)=W(jj)+RMM(kk)*F(ii)
#ifndef GPU       
c       write(*,*) 'W(',jj,') RMM(',kk-M18+1,')=',RMM(kk)
c       write(*,*) F(ii),RMM(kk)*F(ii),ii-1
#endif
 154    continue
 153    continue
c       endif
c ---------------------------------------
        do 61 ii=1,NCO
#ifndef GPU
c       write(*,*) 'W(',ii-1,')=',W(ii)
#endif
 61    DENS=DENS+W(ii)**2
c
      DENS=DENS*2.D0
c
      return
      endif
c

      end
C 
