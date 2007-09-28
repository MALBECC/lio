c-------------------------------------------------------------------
c Calculation of electrical potential in a grid
c
c calls elec.f for evaluating the potential using the
c final density
c The grid is chosen according to the prescriptions
c given in : J. Comp. Chem 11 361 (1990)
c
c A cube of points containing the molecule (9.0 a.u) from
c all atoms is generated. All points inside a predefined
c Van der Waals radius of an atom are discarded, as well as
c the points outside the 9.00 a.u maximum radius
c
c After generating the grid, the electrostatical potential
c in all points of the grid is evaluated by calling the
c subroutine elec
c
c Point charges are evaluated by using a least-squares method
c
c
c values obtained could be printed or used in order
c to fit a set of nuclear centered point charges
c 
c loop over all basis functions
c-------------------------------------------------------------------
      subroutine charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,RMM,map,Q)
c
      implicit real*8 (a-h,o-z)
      INCLUDE 'param'
      logical NORM,OUT
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
      dimension r(nt,3),nshell(0:3)
      dimension RMM(*),xi(3),vdw(54),map(natom)
c
      COMMON /TABLE/ STR(880,0:21)
c data in angstroms
      data vdw /1.50,1.50,1.50,2.00,2.00,2.00,2.00,2.00,2.00,2.00,
     > 2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,
     > 2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,
     > 2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,
     > 2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50/
c
c
      do i=1,54
       vdw(i)=vdw(i)/0.529177
      enddo
c
      MM=M*(M+1)/2
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, F also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
      M11=M9+MMd
c
c auxiliar quantities
c
c determination of the number of classes
c  caso particular HCl - H2O --------------------
      natom1=natom+1
      Q=Q+1.0
      write(*,*) natom1,r(natom1,1),r(natom1,2),r(natom1,3),Q
c -----------------------------------------------
       nclass=1
       do n=1,natom1
        if (map(n).gt.nclass) then
          nclass=map(n)
        endif
       enddo

       do kk=M3,M9
        RMM(kk)=0.0D0
       enddo
c
c finds minimum and maximum for X,Y,Z in order to define the box
       Xmin=0.0D0
       Xmax=0.0D0
       Ymin=0.0D0
       Ymax=0.0D0
       Zmin=0.0D0
       Zmax=0.0D0

       do n=1,natom1
        if (r(n,1).gt.Xmax) then
         Xmax=r(n,1)
        endif
        if (r(n,1).lt.Xmin) then
         Xmin=r(n,1)
        endif
c
        if (r(n,2).gt.Ymax) then
         Ymax=r(n,2)
        endif
        if (r(n,2).lt.Ymin) then
         Ymin=r(n,2)
        endif
c
        if (r(n,3).gt.Zmax) then
         Zmax=r(n,3)
        endif
        if (r(n,3).lt.Zmin) then
         Zmin=r(n,3)
        endif
c 
        enddo
c
       Xmin=Xmin-19.00D0
       Xmax=Xmax+19.00D0
       Ymin=Ymin-19.00D0
       Ymax=Ymax+19.00D0
       Zmin=Zmin-19.00D0
       Zmax=Zmax+19.00D0
c
c GRID GENERATION ------------------------------------------------
c
       Npoint=0
       Ndis=0
       Ndis1=0
       do 10  x=Xmin,Xmax,2.00D0
       do 10  y=Ymin,Ymax,2.00D0
       do 10  z=Zmin,Zmax,2.00D0
c
       xi(1)=x
       xi(2)=y
       xi(3)=z
c
       OUT=.false.
       iout=0
c
       do 111 n=1,natom1
       dd=(xi(1)-r(n,1))**2+(xi(2)-r(n,2))**2+(xi(3)-r(n,3))**2
       dd=sqrt(dd)
       if (dd.lt.vdw(Iz(n))) then
        OUT=.true.
       endif
c
       if (dd.gt.19.0D0) then
        iout=iout+1
       endif
 111   continue
c
       if (iout.eq.natom1) then
        Ndis1=Ndis1+1
       endif

       if ((OUT).or.iout.eq.natom1) then
        Ndis=Ndis+1
        goto 99
       endif
c
       Npoint=Npoint+1
       Nx=3*(Npoint-1)
       RMM(M11+Nx)=x
       RMM(M11+Nx+1)=y
       RMM(M11+Nx+2)=z
c
 99    continue
 10    continue
c
       M15=M11+3*Npoint+1
c CALL ELEC TO COMPUTE ELECTROSTATIC POTENTIAL IN ALL POINTS OF THE GRID
       call elec(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,RMM,Npoint)
c
c----------------------------------------------------------
       do 101  np=1,Npoint
c
       do 2 n=1,nclass
        RMM(M3+n-1)=0.0D0
 2     continue
c
       Nx=3*(np-1)
       xi(1)=RMM(M11+Nx)
       xi(2)=RMM(M11+Nx+1)
       xi(3)=RMM(M11+Nx+2)
c
c tmp weight
       tmp=1.0D0
c
c find out last term contribution in order to obtain constrained charges
c
       Flast=0.0D0
       Nlast=0
       do 7 n=1,natom1
       if (map(n).eq.nclass) then
        dd=(xi(1)-r(n,1))**2+(xi(2)-r(n,2))**2+(xi(3)-r(n,3))**2
        dd=sqrt(dd)
        Flast=Flast+1.0D0/dd
        Nlast=Nlast+1
       endif
 7     continue
       Flast=Flast/Nlast
c
       V=0.0D0
       do 11 n=1,natom1
       dd=(xi(1)-r(n,1))**2+(xi(2)-r(n,2))**2+(xi(3)-r(n,3))**2
       dd=sqrt(dd)
       V=V-Iz(n)/dd
c fitting function : sum over all atoms of the same class
       RMM(M3+map(n)-1)=RMM(M3+map(n)-1)+1.D0/dd-Flast
  11   continue
c

       V=V+RMM(M15+np)+Flast*Q
c
c   LEAST SQUARES
#ifdef essl
c ESSL requires lower packed form
       kk=M7-1
       do 12 j2=1,nclass-1
        tt=RMM(M3+j2-1)*tmp
        RMM(M5+j2-1)=RMM(M5+j2-1)-V*tt
c 
        do 12 i2=j2,nclass-1
        kk=kk+1
        RMM(kk)=RMM(kk)+RMM(M3+i2-1)*tt
 12     continue
#endif
c
#ifdef pack
c LAPACK requires upper packed form
       kk=M7-1
       do 13 j2=1,nclass-1
        tt=RMM(M3+j2-1)*tmp
        RMM(M5+j2-1)=RMM(M5+j2-1)-V*tt
c 
        do 13 i2=1,j2
        kk=kk+1
        RMM(kk)=RMM(kk)+RMM(M3+i2-1)*tt
 13     continue
#endif
c
 101    continue
c
c ESSL OPTION
#ifdef essl
      CALL DPPF(RMM(M7),nclass-1,1)
      CALL DPPS(RMM(M7),nclass-1,RMM(M5),1)
#endif
c
c LAPACK OPTION
#ifdef pack
      call dppco(RMM(M7),nclass-1,rcond,RMM(M9),info)
      call dppsl(RMM(M7),nclass-1,RMM(M5))
#endif
c
      write(*,*) 'CALCULATED POINT CHARGES (A.U)'
      write(*,770)
      
c
       qq=0.0
      do n=1,natom1
       if (map(n).ne.0) then
       qq=qq+RMM(M5+map(n)-1)
       endif
      enddo
c
      qq=Q-qq
      RMM(M5+nclass-1)=qq/Nlast
c
      do n=1,natom1
       if (map(n).ne.0) then
       write(*,760) n,Iz(n),RMM(M5+map(n)-1)
       else
       write(*,760) n,Iz(n),0.0D0
       endif
      enddo
c
c CALCULATION OF FITTED DIPOLE MOMENT
c
c Calculation of center of charge
       xcc=0.0D0
       ycc=0.0D0
       zcc=0.0D0
c
       do n=1,natom1
       xcc=xcc+r(n,1)*abs(RMM(M5+map(n)-1))
       ycc=ycc+r(n,2)*abs(RMM(M5+map(n)-1))
       zcc=zcc+r(n,3)*abs(RMM(M5+map(n)-1))
       enddo
c
       if (Q.ne.0.0D0) then
       xcc=xcc/Q
       ycc=ycc/Q
       zcc=zcc/Q
       endif
c
c      write(*,*) xcc,ycc,zcc
       dx=0.0
       dy=0.0
       dz=0.0
      do n=1,natom1
       dx=dx+RMM(M5+map(n)-1)*(r(n,1)-xcc)
       dy=dy+RMM(M5+map(n)-1)*(r(n,2)-ycc)
       dz=dz+RMM(M5+map(n)-1)*(r(n,3)-zcc)
      enddo
c

      dx=dx*2.54
      dy=dy*2.54
      dz=dz*2.54
      d=sqrt(dx**2+dy**2+dz**2)
c
      write(*,*)
      write(*,*) 'FITTED DIPOLE MOMENT (DEBYES)'
      write(*,900) dx,dy,dz,d
c     return
c     write(*,*) Xmin,Xmax
c     write(*,*) Ymin,Ymax
c     write(*,*) Zmin,Zmax
c     write(*,*) Npoint,Ndis,Ndis1

 760  format(I3,9x,I3,4x,F11.6)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'CHARGE')
 900  format(3(F10.4,2x),2x,F10.4)
      end
