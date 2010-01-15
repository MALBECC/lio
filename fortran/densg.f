c function calculating density functionals
c to be used only with non local density functionals
c gives , as a result, gradient and second derivatives of the
c density
c 19-1-1993
c
      SUBROUTINE DENSG(Dx,Dy,Dz,Xi,ds,NORM,Nuc,ncont,nshell,
     >                 a,c,r,M,M18,NCO,RMM,natom)
c
      implicit real*8 (a-h,o-z)
      logical NORM
      INCLUDE 'param'
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng)
      dimension r(nt,3),nshell(0:3),Xi(3)
      dimension ds(ntq),F(ng),W(ng)
c gradients
      dimension Fg(ng,3),RMM(*)
      dimension Wx(ntq,ng),Wy(ntq,ng),Wz(ntq,ng)
      dimension Dx(ntq),Dy(ntq),Dz(ntq)
c
      common /Ll/ Ll(3)
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
      do i=1,ntq
       Dx(i)=0.0D0
       Dy(i)=0.0D0
       Dz(i)=0.0D0
      enddo
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
c
c basis functions evaluated at r are calculated
c
      do 2 j=1,natom
       do 1 i=1,M
         Wx(j,i)=0.0D0
         Wy(j,i)=0.0D0
         Wz(j,i)=0.0D0
  1    continue
  2   continue
c
      do i=1,M
        W(i)=0.D0
        Fg(i,1)=0.D0
        Fg(i,2)=0.D0
        Fg(i,3)=0.D0
        F(i)=0.D0
      enddo
c
c --- s  case -------
      do 10 i=1,ns
c
      di=ds(Nuc(i))
c
      tx=(xi(1)-r(Nuc(i),1))*2.D0
      ty=(xi(2)-r(Nuc(i),2))*2.D0
      tz=(xi(3)-r(Nuc(i),3))*2.D0
c
      do 15 ni=1,ncont(i)
c
      rexp=a(i,ni)*di
c
      if (rexp.gt.30.) go to 16
      t=exp(-rexp)
      term=t*c(i,ni)
      F(i)=F(i)+term
      Fg(i,1)=Fg(i,1)+term*a(i,ni)
  16  continue
  15  continue
c
      Fg(i,2)=Fg(i,1)*ty
      Fg(i,3)=Fg(i,1)*tz
      Fg(i,1)=Fg(i,1)*tx
  10  continue
c
c--- p  case -------------
      do 920 i=ns+1,ns+np,3
c
      di=ds(Nuc(i))
c
      tx=(xi(1)-r(Nuc(i),1))*2.D0
      ty=(xi(2)-r(Nuc(i),2))*2.D0
      tz=(xi(3)-r(Nuc(i),3))*2.D0
c
      do 20 ni=1,ncont(i)

      rexp=a(i,ni)*di
      if (rexp.gt.30.) goto 21
      t=exp(-rexp)*c(i,ni)
      l1=1
c
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      ii=i
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t
c
c ---------------
      l1=2
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      ii=i+1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t
c
c ------------------
      l1=3
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      ii=i+2
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t
c
  21  continue
  20  continue
 920  continue
c
c-- d case  ------------
      do 940 i=ns+np+1,M,6
c
      di=ds(Nuc(i))
c
      tx=(xi(1)-r(Nuc(i),1))*2.D0
      ty=(xi(2)-r(Nuc(i),2))*2.D0
      tz=(xi(3)-r(Nuc(i),3))*2.D0
c
      do 40 ni=1,ncont(i)
c
      rexp=a(i,ni)*di
c
      if (rexp.gt.30.) goto 41
      t0=exp(-rexp)*c(i,ni)
c -------------------
      l1=1
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=1
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0*fc
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t*t2
      Fg(ii,l2)=Fg(ii,l2)-t*t1
c
c checked
c-------------
      l1=2
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=1
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t*t2
      Fg(ii,l2)=Fg(ii,l2)-t*t1
c
c-------------
      l1=2
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=2
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0*fc
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t*t2
      Fg(ii,l2)=Fg(ii,l2)-t*t1
c-------------
      l1=3
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=1
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t*t2
      Fg(ii,l2)=Fg(ii,l2)-t*t1
c-------------
      l1=3
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=2
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0
c
      term=t*t1*t2
c
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t*t2
      Fg(ii,l2)=Fg(ii,l2)-t*t1
c-------------
      l1=3
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=3
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0*fc
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)+tx*term
      Fg(ii,2)=Fg(ii,2)+ty*term
      Fg(ii,3)=Fg(ii,3)+tz*term
c
      Fg(ii,l1)=Fg(ii,l1)-t*t2
      Fg(ii,l2)=Fg(ii,l2)-t*t1
c------------
c
  41  continue
  40  continue
 940  continue
c
c
********
        kk=0
c
         do 52 j=1,NCO
         do 52 i=1,M
         kk=kk+1
         tmp=RMM(M18+kk-1)
          W(j)=W(j)+tmp*F(i)
          Wx(Nuc(i),j)=Wx(Nuc(i),j)+tmp*Fg(i,1)
          Wy(Nuc(i),j)=Wy(Nuc(i),j)+tmp*Fg(i,2)
          Wz(Nuc(i),j)=Wz(Nuc(i),j)+tmp*Fg(i,3)
 52     continue
c
       do 61 k=1,natom
         do 60 j=1,NCO
           Dx(k)=Dx(k) + W(j)*Wx(k,j)
           Dy(k)=Dy(k) + W(j)*Wy(k,j)
           Dz(k)=Dz(k) + W(j)*Wz(k,j)
 60      continue
        Dx(k)=4.D0*Dx(k)
        Dy(k)=4.D0*Dy(k)
        Dz(k)=4.D0*Dz(k)
 61   continue
      return
      end


