c-------------------------------------------------------------------
c INTEGRAL CALCULATIONS FOR THE SOLVENT POINT CHARGES ELECTROSTATIC
c INTERACTIONS WITH THE SOLUTE ELECTRONIC DENSITY
c
c Dario Estrin
c Buenos Aires, August. 1994.
c
c 1 e integrals
c using the Obara-Saika recursive method.
c are added into the Fock matrix
c
c It's the same as int.f, but using the solvent atoms partial charges
c
c 
c-------------------------------------------------------------------
      subroutine intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,E1s,FQQ,IT,ITEL,NIN,
     >            IPR1,EAC,NPAS)
c
      implicit real*8 (a-h,o-z)
      logical NORM
      include 'param'
c       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0,a0=0.52918 D0)
      PARAMETER (NAT=6000)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
      dimension r(nt,3),nshell(0:3),pc(nt),Ncaj(3)
      dimension xi(3)
      dimension RMM(*),ttna(natom+nsol*natsol),tterm(natom+nsol*natsol)

      COMMON /TABLE/ STR(880,0:21)
      COMMON/BLOC08/FX(NAT),FY(NAT),FZ(NAT),RCT,RCTSQ,RKAPPA,RKAPPA2
      COMMON/BLOC11/BXLGTH
      INTEGER SPC
      COMMON/tipsol/SPC
c auxiliar quantities

      dimension Q(3),d(ntq,ntq),s0s(nt),s1s(nt),s2s(nt),s3s(nt),
     >   s4s(nt),Ll(3),dw(nt,3),jwatc(nt)
      DIMENSION FQQ(natom+nsol*natsol),EAC(natom+nsol*natsol)

c--------------------------------------------------------------------
c Modificacion para correr SPC o TIP4P desde afuera
c---------------------------------------------------------------------
c     write(*,*)'SPC=',SPC
      IF(SPC.EQ.1)THEN

        alpha=1.00D0
        alpha2=0.00D0

      ELSE

        alpha=0.7439762D0
        alpha2=0.1280119D0

      ENDIF
    
c----------------------------------------------------------------------

c distance between pairs of centers
c
c datos TIP4P ----------------
c corresponde a desplazar carga negativa 0.15A en direccion de los H
c       alpha=0.7439762D0
c       alpha2=0.1280119D0
c caso SPC
c       alpha=1.00D0
c       alpha2=0.00D0
c ----------------------------
c      do i=1, 10000
c       write(45,*) i,RMM(i)
c      enddo
c      goto 999 
      DO I=NATOM+1,NATOM+NSOL*NATSOL
      TTNA(I)=0.D0
      TTERM(I)=0.D0
      FQQ(I)=0.D0
      EAC(I)=0.D0
      ENDDO

      box=BXLGTH/A0 
c      write(*,*) 'box',box,RCTSQ
c      RCTSQ2=(box/2.)**2
      RCTSQ2=RCTSQ/(a0**2)
C      RCTSQ2=1000000000.
      if (NORM) then
      sq3=sqrt(3.D0)
      else
      sq3=1.D0
      endif
c
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2
c
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      M2=2*M
c
c Pointers
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, F also uses the same position after S was used
      M11=M3+MM
c now G
      M7=M11+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c
c      do k=natom+1,natom+nsol*natsol
c       pc(k)=0.D0
c      enddo
c      pause



      do 50 i=1,natom
      do 50 j=1,natom
       d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
 50   continue
c
c---------------------------------------------------------------
c
      E1s=0.0D0
c first loop (s|s) case -------------------------------------------
c
      do 200 i=1,ns
      do 200 j=1,i
c
      dd=d(Nuc(i),Nuc(j))
c
      do 200 ni=1,ncont(i)
      do 200 nj=1,ncont(j)
c
c (0|0) calculation
      zij=a(i,ni)+a(j,nj)
      alf=a(i,ni)*a(j,nj)/zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      ccoef=c(i,ni)*c(j,nj)


c
C-------nano minima distancia y cutoff  
      
      nwatc=0
      do 333 iw=natom+1,natsol*Nsol+natom,3
        do   jw=1,3
      Dw(iw,jw)= (r(iw,jw) - Q(jw))/box
      Ncaj(jw)= ANINT(Dw(iw,jw))
c      Ncaj=0 
      Dw(iw,jw)= (Dw(iw,jw)-Ncaj(jw))*box
       enddo
      DDW= Dw(iw,1)**2+ Dw(iw,2)**2+ Dw(iw,3)**2 

      if(DDW.LT.RCTSQ2) then
      nwatc=nwatc + 1
      jwatc(nwatc)= iw
      endif
      do jw=1,3
      do kw=1,2
      iw2=iw+kw
      Dw(iw2,jw)= (r(iw2,jw) - Q(jw))/box
      Dw(iw2,jw)= (Dw(iw2,jw)- Ncaj(jw))*box
      enddo
      enddo

C---- TIP4P
      do jw=1,3
      Dw(iw,jw)=alpha*Dw(iw,jw)+alpha2*(Dw(iw+1,jw)+dw(iw+2,jw))
      enddo


333   continue
c      write(76,*) nwatc*3
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss
c
      k=i+((M2-j)*(j-1))/2
c 
       tna=0.D0
      do 202 n1=1,nwatc
      do 202 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c       write(76,*) 'O',Dw(j1,1)*A0,Dw(j1,2)*A0,Dw(j1,3)*A0  

        tx=-Dw(j1,1)
        ty=-Dw(j1,2)
        tz=-Dw(j1,3)
c      write(*,*)'dist',j1,tx,ty,tz,PC(j1)
       u=tx**2 + ty**2 + tz**2
       u=u*zij

       s0s(j1)=pc(j1)*temp*FUNCT(0,u)
       tna=tna-s0s(j1)
       FQQ(J1)=FQQ(j1)+temp*FUNCT(0,u)*CCOEF*RMM(K) 
c       EAC(J1)=EAC(J1)-s0s(j1)*ccoef*RMM(k)
       EAC(J1)=EAC(J1)-temp*FUNCT(0,u)*ccoef*RMM(k)*pc(j1)

 202   continue
c       pause

      term=ccoef*tna
      RMM(M11+k-1)=RMM(M11+k-1)+ term 
      E1s=E1s+RMM(k)*term

c      write(*,*)'nwatc y mas..',nwatc,j1,pc(j1),Dw(j1,1),tna,ccoef
    
 200  continue
c      write(*,*)'intsol1 E1s',E1s,rmm(k),term,PC(100)
c      pause
c      write(*,*)'fqq (0|0) ',(fqq(i),i=natom+1,natom+10)
c
c------------------------------------------------------------------
c
c second loop  (p|s) case

c
      do 300 i=ns+1,ns+np,3
      do 300 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
c
      do 300 ni=1,ncont(i)
      do 300 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss
c

C-------nano minima distancia y cutoff
 
      nwatc=0
      do 334 iw=natom+1,natsol*Nsol+natom,3
        do   jw=1,3
      Dw(iw,jw)= (r(iw,jw) - Q(jw))/box
      Ncaj(jw)=ANINT(Dw(iw,jw))
      Dw(iw,jw)= (Dw(iw,jw)- Ncaj(jw))*box
 
       enddo
      DDW= Dw(iw,1)**2+ Dw(iw,2)**2+ Dw(iw,3)**2
 
      if(DDW.LT.RCTSQ2) then
      nwatc=nwatc + 1
      jwatc(nwatc)= iw
      endif
      do jw=1,3
      do kw=1,2
      iw2=iw+kw
      Dw(iw2,jw)= (r(iw2,jw) - Q(jw))/box      
      Dw(iw2,jw)= (Dw(iw2,jw)- Ncaj(jw))*box 
      enddo 
      enddo
      do jw=1,3
      Dw(iw,jw)=alpha*Dw(iw,jw)+alpha2*(Dw(iw+1,jw)+dw(iw+2,jw))
      enddo 

334   continue
c loop over nuclei, part common for all shell
c       write(*,*)'nwatc',nwatc
      do 302 n1=1,nwatc
      do 302 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c       
c 

        tx=-Dw(j1,1)
        ty=-Dw(j1,2)
        tz=-Dw(j1,3)
c
c
c
       u= tx**2 +ty**2 +tz**2
       u=u*zij

       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
c
 302  continue
c
      ccoef=c(i,ni)*c(j,nj)
c
c l2: different p in the p shell ( x,y,z respectively)
c
      do 305 l2=1,3
        ii=i+l2-1
        t1=Q(l2)-r(Nuc(i),l2)
c ii index , taking into account different components of the shell
        k=ii+((M2-j)*(j-1))/2
c
c loop over nuclei, specific part
       tna=0.D0
      do 303 n1=1,nwatc
      do 303 k1=1,natsol
       j1 = jwatc(n1) +k1 -1
c
c


c       
c 
c
       t2=-Dw(j1,l2)
       term=t1*s0s(j1)-t2*s1s(j1)
       FQQ(J1)=FQQ(j1)+TERM*CCOEF*RMM(K)
c       EAC(J1)=EAC(J1)-PC(K1)*term*ccoef*RMM(k)
       EAC(J1)=EAC(J1)-TERM*CCOEF*RMM(k)*pc(j1)
       tna=tna-pc(j1)*term
       if(pc(j1).gt.1.or.pc(j1).lt.-1.5) pause
c       write(*,*)'cacas',j1,term,t1*s0s(j1),t2*s1s(j1)
c
 303  continue

        term=ccoef*tna
        RMM(M11+k-1)=RMM(M11+k-1)+term
        E1s=E1s+RMM(k)*term
 305    continue
c

 300  continue       
c        write(*,*)' 2.E1s ',E1s
c-------------------------------------------------------------------
        
C       goto 999
c 
c (p|p) case
c
      do 400 i=ns+1,ns+np,3
      do 400 j=ns+1,i,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 400 ni=1,ncont(i)
      do 400 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss

C-------nano minima distancia y cutoff
 
      nwatc=0
      do 335 iw=natom+1,natsol*Nsol+natom,3
        do   jw=1,3
      Dw(iw,jw)= (r(iw,jw) - Q(jw))/box
      Ncaj(jw)=ANINT(Dw(iw,jw))
      Dw(iw,jw)= (Dw(iw,jw)- Ncaj(jw))*box
 
       enddo
      DDW= Dw(iw,1)**2+ Dw(iw,2)**2+ Dw(iw,3)**2
 
      if(DDW.LT.RCTSQ2) then
      nwatc=nwatc + 1
      jwatc(nwatc)= iw
      endif
      do jw=1,3
      do kw=1,2
      iw2=iw+kw
      Dw(iw2,jw)= (r(iw2,jw) - Q(jw))/box
      Dw(iw2,jw)= (Dw(iw2,jw)- Ncaj(jw))*box
      enddo
      enddo
      do jw=1,3
      Dw(iw,jw)=alpha*Dw(iw,jw)+alpha2*(Dw(iw+1,jw)+Dw(iw+2,jw))
      enddo


335   continue



c
c loop over nuclei, part common for all shell
      do 402 n1=1,nwatc
      do 402 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
 
cc       
c 

       xi(1)=-Dw(j1,1)
       xi(2)=-Dw(j1,2)
       xi(3)=-Dw(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij

       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
c
 402  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c loop over partial charges ( specific part)
      j1=natom
      do 403 n1=1,nwatc
      do 403 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c TIP4P case
c para el O desplazamiento del sitio
c
c
c
c
      do 406 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
c
       t2=-Dw(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)

c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
       do 406 l2=1,lij
       t1=Q(l2)-r(Nuc(j),l2)
       t2=-Dw(j1,l2)
       tna=t1*p0s-t2*p1s
       TTNA(J1)=TNA
       
c
       if (l1.eq.l2) then
        tna=tna+(s0s(j1)-s1s(j1))/z2
        TTNA(J1)=TTNA(J1)+(s0s(j1)-s1s(j1))/z2
       endif
c
        tna1=tna*pc(j1)
c
       ii=i+l1-1
       jj=j+l2-1
       k=ii+((M2-jj)*(jj-1))/2
       term=-tna1*ccoef
       TTERM(J1)=TTNA(j1)*CCOEF
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term

       FQQ(J1)=FQQ(j1)+TTERM(j1)*RMM(K)
c       EAC(J1)=EAC(J1)+RMM(K)*term
       EAC(J1)=EAC(J1)-TTERM(J1)*RMM(k)*pc(j1)

 406  continue
c

 403   continue
c ---------------
 400  continue
c       write(*,*)'3.E1s ',E1s,rmm(k),term
c
c-------------------------------------------------------------------

       do i=1+natom,natom+nsol*natsol
        ttna(i)=0.d0
        tterm(i)=0.d0
       enddo

c (d|s) case
c
      do 500 i=ns+np+1,M,6   
      do 500 j=1,ns
c
      dd=d(Nuc(i),Nuc(j))
c
      do 500 ni=1,ncont(i)
      do 500 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss

C-------nano minima distancia y cutoff
 
      nwatc=0
      do 336 iw=natom+1,natsol*Nsol+natom,3
        do   jw=1,3
      Dw(iw,jw)= (r(iw,jw) - Q(jw))/box
      Ncaj(jw)=ANINT(Dw(iw,jw))
      Dw(iw,jw)= (Dw(iw,jw)- Ncaj(jw))*box 
       enddo
      DDW= Dw(iw,1)**2+ Dw(iw,2)**2+ Dw(iw,3)**2
 
      if(DDW.LT.RCTSQ2) then
      nwatc=nwatc + 1
      jwatc(nwatc)= iw
      endif
      do jw=1,3
      do kw=1,2
      iw2=iw+kw
      Dw(iw2,jw)= (r(iw2,jw) - Q(jw))/box
      Dw(iw2,jw)= (Dw(iw2,jw)- Ncaj(jw))*box
      enddo
      enddo
      do jw=1,3
      Dw(iw,jw)=alpha*Dw(iw,jw)+alpha2*(Dw(iw+1,jw)+dw(iw+2,jw))
      enddo


336   continue



c
c loop over partial charges, part common for all shell
      do 502 n1=1,nwatc
      do 502 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c TIP4P case
c para el O desplazamiento del sitio
c
c
       xi(1)=-Dw(j1,1)
       xi(2)=-Dw(j1,2)
       xi(3)=-Dw(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
c
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
 502  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
                  
      do 503 n1=1,nwatc
      do 503 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c TIP4P case
c para el O desplazamiento del sitio
c
c
c
      do 506 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=-Dw(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
c
      do 506 l2=1,l1
c
       t1=Q(l2)-r(Nuc(i),l2)
       t2=-Dw(j1,l2)
       tna=t1*p0s-t2*p1s 
       TTNA(J1)=T1*P0S-T2*P1S
c
       f1=1.D0
       if (l1.eq.l2) then
        tna=tna+(s0s(j1)-s1s(j1))/z2
        TTNA(J1)=TTNA(J1)+(S0S(J1)-S1S(J1))/Z2
        f1=sq3
       endif
c
       l12=l1*(l1-1)/2+l2
c ordering of d shell should be:
c xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
c
       ii=i+l12-1
c
       k=ii+((M2-j)*(j-1))/2
       cc=ccoef/f1
       term=-cc*tna*pc(j1)
       TTERM(J1)=CC*TTNA(J1)
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
       FQQ(J1)=FQQ(j1)+TTERM(J1)*RMM(K)
c       EAC(J1)=EAC(J1)+RMM(K)*term
       EAC(J1)=EAC(J1)-TTERM(j1)*RMM(k)*pc(j1)

 506  continue
c
 503  continue
 500  continue
c       write(*,*)'4.E1s ',E1s,rmm(k),term
c-----------------------------------------------------------------
c
c (d|p) case 
c
      do 600 i=ns+np+1,M,6
      do 600 j=ns+1,ns+np,3
c
      dd=d(Nuc(i),Nuc(j))
c
      do 600 ni=1,ncont(i)
      do 600 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       temp=2.D0*sqrt(zij/pi)*ss

C-------nano minima distancia y cutoff
 
      nwatc=0
      do 337 iw=natom+1,natsol*Nsol+natom,3
        do   jw=1,3
      Dw(iw,jw)= (r(iw,jw) - Q(jw))/box
      Ncaj(jw)=ANINT(Dw(iw,jw))
      Dw(iw,jw)= (Dw(iw,jw)- Ncaj(jw))*box 
       enddo
      DDW= Dw(iw,1)**2+ Dw(iw,2)**2+ Dw(iw,3)**2
 
      if(DDW.LT.RCTSQ2) then
      nwatc=nwatc + 1
      jwatc(nwatc)= iw
      endif
      do jw=1,3
      do kw=1,2
      iw2=iw+kw
      Dw(iw2,jw)= (r(iw2,jw) - Q(jw))/box
      Dw(iw2,jw)= (Dw(iw2,jw)- Ncaj(jw))*box
      enddo
      enddo
      do jw=1,3
      Dw(iw,jw)=alpha*Dw(iw,jw)+alpha2*(Dw(iw+1,jw)+dw(iw+2,jw))
      enddo

337   continue



c
c loop over nuclei, part common for all shell
      do 602 n1=1,nwatc
      do 602 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c
       xi(1)=-Dw(j1,1)
       xi(2)=-Dw(j1,2)
       xi(3)=-Dw(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
c
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
       s3s(j1)=temp*FUNCT(3,u)
 602  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
      do 603 n1=1,nwatc
      do 603 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c TIP4P case
c para el O desplazamiento del sitio
c
c
c
      do 606 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=-Dw(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
c
      do 606 l2=1,l1
       t1=Q(l2)-r(Nuc(i),l2)
       t2=-Dw(j1,l2)
       pj0s=t1*s0s(j1)-t2*s1s(j1)
       pj1s=t1*s1s(j1)-t2*s2s(j1)
c
       f1=1.D0
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(j1)-s1s(j1))/z2
        d1s=d1s+(s1s(j1)-s2s(j1))/z2
       endif
c
c
      do 606 l3=1,3
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=-Dw(j1,l3)
       tna=t1*d0s-t2*d1s
c
       if (l1.eq.l3) then
        tna=tna+(pj0s-pj1s)/z2
       endif
c
       if (l2.eq.l3) then
        tna=tna+(p0s-p1s)/z2
       endif
c
c
c
       l12=l1*(l1-1)/2+l2
       ii=i+l12-1
       jj=j+l3-1
c
      k=ii+((M2-jj)*(jj-1))/2
       cc=ccoef/f1
       term=-cc*tna*pc(j1)
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
       FQQ(J1)=FQQ(j1)-TERM*RMM(K)/PC(J1)
       EAC(J1)=EAC(J1)+RMM(K)*term

 606  continue
c
 603  continue
c
 600  continue
c       write(*,*)'6.E1s ',E1s,rmm(k),term
c
c-------------------------------------------------------------------
c
c (d|d) case
c
      do 700 i=ns+np+1,M,6
      do 700 j=ns+np+1,i,6
c
      dd=d(Nuc(i),Nuc(j))
c
      do 700 ni=1,ncont(i)
      do 700 nj=1,ncont(j)
c
      zij=a(i,ni)+a(j,nj)
      z2=2.D0*zij
      Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
      Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
      Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
      alf=a(i,ni)*a(j,nj)/zij
      ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
      temp=2.D0*sqrt(zij/pi)*ss

 
C-------nano minima distancia y cutoff
 
      nwatc=0
      do 338 iw=natom+1,natsol*Nsol+natom,3
        do   jw=1,3
      Dw(iw,jw)= (r(iw,jw) - Q(jw))/box
      Ncaj(jw)=ANINT(Dw(iw,jw))
      Dw(iw,jw)= (Dw(iw,jw)- Ncaj(jw))*box 
       enddo
      DDW= Dw(iw,1)**2+ Dw(iw,2)**2+ Dw(iw,3)**2
 
      if(DDW.LT.RCTSQ2) then
      nwatc=nwatc + 1
      jwatc(nwatc)= iw
      endif
      do jw=1,3
      do kw=1,2
      iw2=iw+kw
      Dw(iw2,jw)= (r(iw2,jw) - Q(jw))/box
      Dw(iw2,jw)= (Dw(iw2,jw)- Ncaj(jw))*box
      enddo
      enddo
      do jw=1,3
      Dw(iw,jw)=alpha*Dw(iw,jw)+alpha2*(Dw(iw+1,jw)+dw(iw+2,jw))
      enddo

338   continue




c loop over nuclei, part common for all shell
                  
      do 702 n1=1,nwatc
      do 702 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c
       xi(1)=-Dw(j1,1)
       xi(2)=-Dw(j1,2)
       xi(3)=-Dw(j1,3)
c
       u=xi(1)**2+xi(2)**2+xi(3)**2
       u=u*zij
c
       s0s(j1)=temp*FUNCT(0,u)
       s1s(j1)=temp*FUNCT(1,u)
       s2s(j1)=temp*FUNCT(2,u)
       s3s(j1)=temp*FUNCT(3,u)
       s4s(j1)=temp*FUNCT(4,u)
 702  continue
c
c
      ccoef=c(i,ni)*c(j,nj)
c
c Loop over partial charges
                
      do 703 n1=1,nwatc
      do 703 k1=1,natsol
       j1=jwatc(n1)+k1-1
c
c
      do 706 l1=1,3
       t1=Q(l1)-r(Nuc(i),l1)
       t2=-Dw(j1,l1)
       p0s=t1*s0s(j1)-t2*s1s(j1)
       p1s=t1*s1s(j1)-t2*s2s(j1)
       p2s=t1*s2s(j1)-t2*s3s(j1)
       p3s=t1*s3s(j1)-t2*s4s(j1)
c
      do 706 l2=1,l1
       f1=1.D0
       t1=Q(l2)-r(Nuc(i),l2)
       t2=-Dw(j1,l2)
       pj0s=t1*s0s(j1)-t2*s1s(j1)
       pj1s=t1*s1s(j1)-t2*s2s(j1)
       pj2s=t1*s2s(j1)-t2*s3s(j1)
c
       d0s=t1*p0s-t2*p1s
       d1s=t1*p1s-t2*p2s
       d2s=t1*p2s-t2*p3s
c
c
       if (l1.eq.l2) then
        f1=sq3
        d0s=d0s+(s0s(j1)-s1s(j1))/z2
        d1s=d1s+(s1s(j1)-s2s(j1))/z2
        d2s=d2s+(s2s(j1)-s3s(j1))/z2
       endif
c
c
      lij=3
      if (i.eq.j) then
       lij=l1
      endif
c
       do 706 l3=1,lij
c
       t1=Q(l3)-r(Nuc(j),l3)
       t2=-Dw(j1,l3)
c
       d0p=t1*d0s-t2*d1s
       d1p=t1*d1s-t2*d2s
c
       pi0p=t1*p0s-t2*p1s
       pi1p=t1*p1s-t2*p2s
       pj0p=t1*pj0s-t2*pj1s
       pj1p=t1*pj1s-t2*pj2s
c
       if (l1.eq.l3) then
        d0p=d0p+(pj0s-pj1s)/z2
        d1p=d1p+(pj1s-pj2s)/z2
        pi0p=pi0p+(s0s(j1)-s1s(j1))/z2
        pi1p=pi1p+(s1s(j1)-s2s(j1))/z2
       endif
c
       if (l2.eq.l3) then
        d0p=d0p+(p0s-p1s)/z2
        d1p=d1p+(p1s-p2s)/z2
        pj0p=pj0p+(s0s(j1)-s1s(j1))/z2
        pj1p=pj1p+(s1s(j1)-s2s(j1))/z2
       endif
c
c
      lk=l3
      if (i.eq.j) then
       lk=min(l3,Ll(l1)-Ll(l3)+l2)
      endif
c
      do 706 l4=1,lk
c
       f2=1.D0
       t1=Q(l4)-r(Nuc(j),l4)
       t2=-Dw(j1,l4)
       tna=t1*d0p-t2*d1p
c
       if (l4.eq.l1) then
        tna=tna+(pj0p-pj1p)/z2
       endif
c
       if (l4.eq.l2) then
        tna=tna+(pi0p-pi1p)/z2
       endif
c
       if (l4.eq.l3) then
        f2=sq3
        tna=tna+(d0s-d1s)/z2
       endif
c
       cc=ccoef/(f1*f2)
       term=-cc*pc(j1)*tna
c
       l12=l1*(l1-1)/2+l2
       l34=l3*(l3-1)/2+l4
       ii=i+l12-1
       jj=j+l34-1
c
       k=ii+((M2-jj)*(jj-1))/2
       RMM(M11+k-1)=RMM(M11+k-1)+term
       E1s=E1s+RMM(k)*term
       FQQ(J1)=FQQ(j1)-TERM*RMM(K)/PC(J1)
       EAC(J1)=EAC(J1)+RMM(K)*term
c
 706  continue
 703  continue
 700  continue
c       write(*,*)'7.E1s ',E1s,rmm(k),term
c

c      IF(NPAS.NE.1)THEN
c        IF(MOD((IT-NIN),NPAS).EQ.0)THEN
c        write(*,*) 'E1s=',E1s
c        ENDIF
c      ENDIF

c      do i=natom+1,natom+natsol*nsol
c      write(*,*)' fq INTSOL ',i,pc(i),-fqq(i)
c      enddo
c      pause
C      do 3 i=1,ns
C      do 3 j=1,i
C      k=i+((M2-j)*(j-1))/2
C      coco= caca(k)-RMM(k)
C3     continue
999    continue
c      write(*,*) 'Esol  posta',E1s
      return
      end
c-------------------------------------------------------------------
