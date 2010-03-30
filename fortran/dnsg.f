c function calculating density functionals
c to be used only with non local density functionals
c gives , as a result, gradient and second derivatives of the
c density
c 19-1-1993
c
      SUBROUTINE DNSG(Dens,F,Dx,Dy,Dz,Dxx,Dyy,Dzz,Dxy,Dyz,Dxz,
     >           Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,M,M18,NCO,RMM)
c
      implicit real*8 (a-h,o-z)
      logical NORM
      INCLUDE 'param'
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng)
      dimension r(nt,3),nshell(0:3),Xi(3)
      dimension ds(ntq),F(M),W(ng),RMM(*)
c gradients
      dimension Fg(ng,3),NN(ng)
      dimension Wx(ng),Wy(ng),Wz(ng)
      dimension Wxx(ng),Wyy(ng),Wzz(ng),Wxy(ng),Wxz(ng),Wyz(ng)
      dimension Fxx(ng),Fyy(ng),Fzz(ng),Fxy(ng),Fxz(ng),Fyz(ng)
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
      Dx=0.0D0
      Dy=0.0D0
      Dz=0.0D0
      Dxx=0.0D0
      Dyy=0.0D0
      Dzz=0.0D0
      Dxy=0.0D0
      Dyz=0.0D0
      Dxz=0.0D0
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
c
c basis functions evaluated at r are calculated
c
       do 1 i=1,M
        W(i)=0.D0
        Wx(i)=0.0D0
        Wy(i)=0.0D0
        Wz(i)=0.0D0
        Wxx(i)=0.0D0
        Wyy(i)=0.0D0
        Wzz(i)=0.0D0
        Wxy(i)=0.0D0
        Wxz(i)=0.0D0
        Wyz(i)=0.0D0
c
        Fg(i,1)=0.D0
        Fg(i,2)=0.D0
        Fg(i,3)=0.D0
        Fxx(i)=0.0D0
        Fyy(i)=0.0D0
        Fzz(i)=0.0D0
        Fxy(i)=0.0D0
        Fxz(i)=0.0D0
        Fyz(i)=0.0D0
 1      F(i)=0.D0
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
      if (rexp.gt.40.D0) go to 16
c mask to indicate if this function will be used or not
c
c-------------------------------------------------------
      t=exp(-rexp)
      term=t*c(i,ni)
      F(i)=F(i)+term
      Fxx(i)=Fxx(i)+term*a(i,ni)**2
      Fg(i,1)=Fg(i,1)-term*a(i,ni)
  16  continue
  15  continue
c
      Fxy(i)=Fxx(i)*tx*ty
      Fxz(i)=Fxx(i)*tx*tz
      Fyz(i)=Fxx(i)*ty*tz
      t1=2.D0*Fg(i,1)
      Fyy(i)=Fxx(i)*ty**2+t1
      Fzz(i)=Fxx(i)*tz**2+t1
      Fxx(i)=Fxx(i)*tx**2+t1
c
      Fg(i,2)=Fg(i,1)*ty
      Fg(i,3)=Fg(i,1)*tz
      Fg(i,1)=Fg(i,1)*tx
  10  continue
c
c--- p  case -------------
      do 20 i=ns+1,ns+np,3
c
      di=ds(Nuc(i))
c
      tx=(xi(1)-r(Nuc(i),1))*2.D0
      ty=(xi(2)-r(Nuc(i),2))*2.D0
      tz=(xi(3)-r(Nuc(i),3))*2.D0
c
      do 20 ni=1,ncont(i)

      rexp=a(i,ni)*di
      if (rexp.gt.40.D0) go to 21
c mask to indicate if this function will be used or not
c
c----------------------------------------------
      t=exp(-rexp)*c(i,ni)
      l1=1
c
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      term4=t*a(i,ni)
      ii=i
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term2=term*a(i,ni)
      term3=2.D0*term
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t
c
      Fxx(ii)=Fxx(ii)-term3+tx**2*term2-tx*term4
      Fxy(ii)=Fxy(ii)+tx*ty*term2
      Fxz(ii)=Fxz(ii)+tx*tz*term2
      Fyy(ii)=Fyy(ii)-term3+ty**2*term2
      Fyz(ii)=Fyz(ii)+ty*tz*term2
      Fzz(ii)=Fzz(ii)-term3+tz**2*term2
c
      Fxx(ii)=Fxx(ii)-tx*term4
      Fxy(ii)=Fxy(ii)-ty*term4
      Fxz(ii)=Fxz(ii)-tz*term4
c ---------------
      l1=2
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      ii=i+1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term2=term*a(i,ni)
      term3=term*2.D0
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t
c
      Fxx(ii)=Fxx(ii)-term3+tx**2*term2
      Fxy(ii)=Fxy(ii)+tx*ty*term2-term4*tx
      Fxz(ii)=Fxz(ii)+tx*tz*term2
      Fyy(ii)=Fyy(ii)-term3+ty**2*term2-term4*ty
      Fyz(ii)=Fyz(ii)+ty*tz*term2
      Fzz(ii)=Fzz(ii)-term3+tz**2*term2
c
      Fyy(ii)=Fyy(ii)-ty*term4
      Fyz(ii)=Fyz(ii)-tz*term4
c ------------------
      l1=3
      t1=xi(l1)-r(Nuc(i),l1)
      term=t*t1
      ii=i+2
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term2=term*a(i,ni)
      term3=term*2.D0
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t
c
      Fxx(ii)=Fxx(ii)-term3+tx**2*term2
      Fxy(ii)=Fxy(ii)+tx*ty*term2
      Fxz(ii)=Fxz(ii)+tx*tz*term2-term4*tx
      Fyy(ii)=Fyy(ii)-term3+ty**2*term2
      Fyz(ii)=Fyz(ii)+ty*tz*term2-term4*ty
      Fzz(ii)=Fzz(ii)-term3+tz**2*term2-term4*tz
c
      Fzz(ii)=Fzz(ii)-tz*term4
c
  21  continue
  20  continue
c
c-- d case  ------------
      do 40 i=ns+np+1,M,6
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
      if (rexp.gt.40.D0) go to 41
c mask to indicate if this function will be used or not
c
c----------------------------------------------
      t0=exp(-rexp)*c(i,ni)
c -------------------
      l1=1
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=1
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0*fc
       term4=t*a(i,ni)
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term3=2.D0*term
      term2=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t*t2
      Fg(ii,l2)=Fg(ii,l2)+t*t1
c
      Fxx(ii)=Fxx(ii)+ tx**2*term2 +2.D0*t -10.D0*term
      Fxy(ii)=Fxy(ii)+tx*ty*(term2 - term4)
      Fxz(ii)=Fxz(ii)+tx*tz*(term2 -term4)
      Fyy(ii)=Fyy(ii)-term3+ty**2*term2
      Fyz(ii)=Fyz(ii)+ty*tz*term2
      Fzz(ii)=Fzz(ii)-term3+tz**2*term2
c
c checked
c-------------
      l1=2
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=1
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0
       term4=t*a(i,ni)/2.D0
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term3=2.D0*term
      term2=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t*t2
      Fg(ii,l2)=Fg(ii,l2)+t*t1
c
      Fxx(ii)=Fxx(ii)-3.D0*term3+tx**2*term2
      Fxy(ii)=Fxy(ii)+tx*ty*term2-term4*(tx**2+ty**2)+t
      Fxz(ii)=Fxz(ii)+tx*tz*term2-term4*ty*tz
      Fyy(ii)=Fyy(ii)-3.D0*term3+ty**2*term2
      Fyz(ii)=Fyz(ii)+ty*tz*term2-term4*tx*tz
      Fzz(ii)=Fzz(ii)-term3+tz**2*term2
c
c-------------
      l1=2
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=2
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0*fc
      term4=t*a(i,ni)
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term3=2.D0*term
      term2=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t*t2
      Fg(ii,l2)=Fg(ii,l2)+t*t1
c
      Fxx(ii)=Fxx(ii) -term3+tx**2*term2
      Fxy(ii)=Fxy(ii)+tx*ty*(term2 - term4)
      Fxz(ii)=Fxz(ii)+tx*tz*term2 
      Fyy(ii)=Fyy(ii) + ty**2*term2 +2.D0*t -10.D0*term
      Fyz(ii)=Fyz(ii)+ty*tz*(term2 -term4)
      Fzz(ii)=Fzz(ii)-term3+tz**2*term2
c-------------
      l1=3
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=1
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0
      term4=t*a(i,ni)/2.D0
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term3=2.D0*term
      term2=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t*t2
      Fg(ii,l2)=Fg(ii,l2)+t*t1
c
      Fxx(ii)=Fxx(ii) -3.D0*term3+tx**2*term2
      Fxy(ii)=Fxy(ii)+tx*ty*term2  -term4*ty*tz
      Fxz(ii)=Fxz(ii)+tx*tz*term2 -term4*(tx**2 +tz**2) +t
      Fyy(ii)=Fyy(ii) + ty**2*term2 -term3
      Fyz(ii)=Fyz(ii)+ty*tz*term2 -term4*tx*ty
      Fzz(ii)=Fzz(ii)-3.D0*term3+tz**2*term2
c-------------
      l1=3
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=2
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0
      term4=t*a(i,ni)/2.D0
c
      term=t*t1*t2
c
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term3=2.D0*term
      term2=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t*t2
      Fg(ii,l2)=Fg(ii,l2)+t*t1
c
      Fxx(ii)=Fxx(ii) -term3+tx**2*term2
      Fxy(ii)=Fxy(ii)+tx*ty*term2 -term4*tz*tx
      Fxz(ii)=Fxz(ii)+tx*tz*term2 -term4*ty*tx
      Fyy(ii)=Fyy(ii) + ty**2*term2 -3.D0*term3
      Fyz(ii)=Fyz(ii)+ty*tz*term2 -term4*( tz**2 +ty**2) +t
      Fzz(ii)=Fzz(ii) + tz**2*term2 -3.D0*term3
c-------------
      l1=3
c
      t1=xi(l1)-r(Nuc(i),l1)
      l2=3
      t2=xi(l2)-r(Nuc(i),l2)
       t=t0*fc
      term4=t*a(i,ni)
c
      term=t*t1*t2
      l12=Ll(l1)+l2
      ii=i+l12-1
c
      F(ii)=F(ii)+term
      term=term*a(i,ni)
      term3=2.D0*term
      term2=term*a(i,ni)
      Fg(ii,1)=Fg(ii,1)-tx*term
      Fg(ii,2)=Fg(ii,2)-ty*term
      Fg(ii,3)=Fg(ii,3)-tz*term
c
      Fg(ii,l1)=Fg(ii,l1)+t*t2
      Fg(ii,l2)=Fg(ii,l2)+t*t1
c
      Fxx(ii)=Fxx(ii) -term3+tx**2*term2
      Fxy(ii)=Fxy(ii)+tx*ty*term2 
      Fxz(ii)=Fxz(ii)+tx*tz*(term2-term4)
      Fyy(ii)=Fyy(ii) - term3+ty**2*term2
      Fyz(ii)=Fyz(ii)+ty*tz*(term2 -term4)
      Fzz(ii)=Fzz(ii)-10.0*term +tz**2*term2 +2.D0*t
c------------
c
  41  continue
  40  continue
c
c now calculation of vector W : density matrix scalar F
c
      if (Ndens.eq.1) then
      k=0
      do 50 j=1,M
      do 51 i=j,M
       k=k+1
      Wx(j)=Wx(j)+RMM(k)*Fg(i,1)
      Wy(j)=Wy(j)+RMM(k)*Fg(i,2)
      Wz(j)=Wz(j)+RMM(k)*Fg(i,3)
c
      Wxx(j)=Wxx(j)+RMM(k)*Fxx(i)
      Wxy(j)=Wxy(j)+RMM(k)*Fxy(i)
      Wxz(j)=Wxz(j)+RMM(k)*Fxz(i)
      Wyz(j)=Wyz(j)+RMM(k)*Fyz(i)
      Wyy(j)=Wyy(j)+RMM(k)*Fyy(i)
      Wzz(j)=Wzz(j)+RMM(k)*Fzz(i)
c
 51   W(j)=W(j)+RMM(k)*F(i)
 50   continue
c
      do 60 i=1,M
c      write(*,*) W(i)
       DENS=DENS+F(i)*W(i)
c       write(*,*) W(i)
       Dx=Dx+Fg(i,1)*W(i)+F(i)*Wx(i)
       Dy=Dy+Fg(i,2)*W(i)+F(i)*Wy(i)
       Dz=Dz+Fg(i,3)*W(i)+F(i)*Wz(i)
c
       Dxx=Dxx+ 2.D0*Fg(i,1)*Wx(i) + Fxx(i)*W(i) +  F(i)*Wxx(i)
       Dyy=Dyy+ 2.D0*Fg(i,2)*Wy(i) + Fyy(i)*W(i) +  F(i)*Wyy(i)
       Dzz=Dzz+ 2.D0*Fg(i,3)*Wz(i) + Fzz(i)*W(i) +  F(i)*Wzz(i)
       Dxy=Dxy+ Fg(i,1)*Wy(i) +Fg(i,2)*Wx(i) +Fxy(i)*W(i) + F(i)*Wxy(i)
       Dxz=Dxz+ Fg(i,1)*Wz(i) +Fg(i,3)*Wx(i) +Fxz(i)*W(i) + F(i)*Wxz(i)
       Dyz=Dyz+ Fg(i,2)*Wz(i) +Fg(i,3)*Wy(i) +Fyz(i)*W(i) + F(i)*Wyz(i)
c
  60  continue
c
      return
c
      else
c
      kk=0
c
c
c
c
      do 92 j=1,NCO
      do 92 i=1,M
c
      kk=kk+1
      tmp=RMM(M18+kk-1)
      Wx(j)=Wx(j)+tmp*Fg(i,1)
      Wy(j)=Wy(j)+tmp*Fg(i,2)
      Wz(j)=Wz(j)+tmp*Fg(i,3)
c
      Wxx(j)=Wxx(j)+tmp*Fxx(i)
      Wxy(j)=Wxy(j)+tmp*Fxy(i)
      Wxz(j)=Wxz(j)+tmp*Fxz(i)
      Wyz(j)=Wyz(j)+tmp*Fyz(i)
      Wyy(j)=Wyy(j)+tmp*Fyy(i)
      Wzz(j)=Wzz(j)+tmp*Fzz(i)
c
      W(j)=W(j)+tmp*F(i)
c
 92   continue
c
      do 93 i=1,NCO
       DENS=DENS+W(i)**2
       Dx=Dx+W(i)*Wx(i)
       Dy=Dy+W(i)*Wy(i)
       Dz=Dz+W(i)*Wz(i)
c
       Dxx=Dxx+ Wx(i)**2+ W(i)*Wxx(i)
       Dyy=Dyy+ Wy(i)**2+ W(i)*Wyy(i)
       Dzz=Dzz+ Wz(i)**2+ W(i)*Wzz(i)
       Dxy=Dxy+ Wx(i)*Wy(i)+ W(i)*Wxy(i)
       Dxz=Dxz+ Wx(i)*Wz(i)+ W(i)*Wxz(i)
       Dyz=Dyz+ Wy(i)*Wz(i)+ W(i)*Wyz(i)
c
  93  continue
c
      DENS=2.D0*DENS
      Dx=4.D0*Dx
      Dy=4.D0*Dy
      Dz=4.D0*Dz
      Dxx=4.D0*Dxx
      Dyy=4.D0*Dyy
      Dzz=4.D0*Dzz
      Dxy=4.D0*Dxy
      Dxz=4.D0*Dxz
      Dyz=4.D0*Dyz
c
      return
      endif
c
      end
