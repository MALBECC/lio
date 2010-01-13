c MD Subroutine ----------------------------------
c at each nuclear configuration, a SCF calculation
c time steps : normal for flexible molecules
c
c
c calls SCF subroutine for calculating energy at each geometry ,
c then calls int1G int3G and intSG, that calculate gradients of
c 1 electron part, 2 electron part and overlap respectively
c After that, a move is done, and the process is continued
c
c all specifications given in the namelist SCF are necessary
c
c Dario Estrin, 1992
c---------------------------------------------------
       subroutine MD2(MEMO,NORM,natom,Iz,r,v,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,name1,ikk,
     >    nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
      implicit real*8 (a-h,o-z)
      integer nopt,iconst,igrid,iforce
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write,SVD
      logical SHFT,OPEN,GRAD,integ,field,sol,free,MEMO
      INCLUDE 'param'
      dimension r(nt,3),nshelld(0:3),nshell(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
      dimension Em(ntq+nss),Rm(ntq+nss),pc(nss),alpha(nss)
c
      common /sol1/ Nsol,natsol,alpha,Em,Rm,pc,sol,free
c
      character*17 name5
      character*17 name6
      character*12 name1
c
c nuclear MD
      dimension v(nt,3),f(nt,3),v1(nt,3),f1(nt,3),r1(nt,3)
      dimension Pm(nt)
      dimension RMM(*)
c
c auxiliars
c X scratch space in matrices
      dimension X(M,3*M),XX(Md,Md)
c
c
      COMMON /TABLE/ STR(880,0:21)
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
c     common /HF/ nopt,OPEN,NMAX,NCO,ATRHO,VCINP,DIRECT,
c    >             IDAMP,EXTR,SHFT,SHI,GOLD,told,write,Nunp
c     common /dyn/ h,nsteps,pmc,target,Pm,scale,scale2
      common /dyn/ h,nsteps,Pm,Nscale
      common /Sys/ SVD,iconst
      common /index/ index(ng)
      common /force/ iforce
      common/Ngeom/ ngeo
      common /ENum/ GRAD
      common /cav/ a0,epsilon,field
c------------------------------------------------------------------
c 
      name5= name1(1:ikk)//'.traj'
      name6= name1(1:ikk)//'.eng'
      open(unit=19,file=name5,form='unformatted')
      open(unit=25,file=name6)
c
c Pointers
c
      write(25,125) 
      NCOa=NCO
      NCOb=NCO+Nunp
      iforce=0
c
      Nel=NCOa+NCOb
      GRAD=.false.
      MM=M*(M+1)/2 
      MM2=M**2
      MMd=Md*(Md+1)/2
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
c now H
      M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
cLeast squares
      M17=M15+MM
c
c vectors of MO alpha
      M18=M17+MMd
c vectors of MO beta
      M18b=M18+M*NCOa
c
      ntom=natom
      if (sol.and.free) then
        ntom=ntom + Nsol*natsol
      endif
c
        if (field) then
        g0=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
        endif
c
      Nlib=(3*ntom-6)
c
      if (ntom.eq.2) then
        Nlib=3*natom-5
      endif
c First step : call SCF subroutine
c
c first call , needs also E , not only gradients -----------
c
      if (.not.OPEN) then
      call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      else
      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      endif
c
c now gradients
c
c
      call int1G(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En,f)
c
c 
      call int3G(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,Exc,f,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write,true)
c
      if (field) then
        g1=g0
        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,Nel,
     >       ux,uy,uz)
        call dipg(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >               Nel,g1,ux,uy,uz,f)
      endif
c
c----------------------------
c classical solvent case ----
        if (sol) then
         call mmsolG(natom,Nsol,natsol,Iz,pc,r,Em,Rm,f)
         call intsolG(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,f)
        endif
c---------------------------
      call intSG(NORM,natom,r,Nuc,M,Md,ncont,nshell,c,a,RMM,f)
c
c ----------------------------------------------------------
c
      E=E+Exc
      ENK=0.0D0
       write(*,480) E
        do i=1,ntom
        write(*,580) i,f(i,1),f(i,2),f(i,3)
        ENK=ENK+(v(i,1)**2+v(i,2)**2+v(i,3)**2)*Pm(i)
        enddo
c
       ENK=1822.7926*ENK/2.D0
       Etot=E+ENK
c--- velocity form of verlet algorithm ---
c
c fac= 1/(2.419E-2 * 1822.7926) : factor to pass from time step
c in femptoseconds and form the mass in m C12/12 in a.u
c
c to pass from fs to a.u of time
      h=h/2.419E-02
c corrected for passing mass in units of mC12/12 to a.u
      dt=h**2
      dt2=.5d0*dt/1822.7926
      dth=.5d0*h/1822.7926
c
c 'Big' Loop over steps -------------------------------------------
      do 50 istep=1,nsteps
c
c
      do 65 i=1,ntom
       r1(i,1)=r(i,1)+h*v(i,1)-dt2*f(i,1)/Pm(i)
       r1(i,2)=r(i,2)+h*v(i,2)-dt2*f(i,2)/Pm(i)
       r1(i,3)=r(i,3)+h*v(i,3)-dt2*f(i,3)/Pm(i)
 65    continue
c
c call forces again
      if (.not.OPEN) then
      call SCF(MEMO,NORM,natom,Iz,r1,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      else
      call SCFop(MEMO,NORM,natom,Iz,r1,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      endif
c
c now gradients
c
c
      call int1G(NORM,natom,r1,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En,f1)
c
      call int3G(NORM,natom,Iz,r1,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,Exc,f1,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write,true)
c
      if (field) then
        g1=g0
        call dip(NORM,Iz,natom,r1,Nuc,M,ncont,nshell,c,a,RMM,Nel,
     >       ux,uy,uz)
        call dipg(NORM,Iz,natom,r1,Nuc,M,ncont,nshell,c,a,RMM,
     >               Nel,g1,ux,uy,uz,f1)
      endif
c
c----------------------------
c classical solvent case ----
        if (sol) then
         call mmsolG(natom,Nsol,natsol,Iz,pc,r1,Em,Rm,f1)
         call intsolG(NORM,natom,Nsol,natsol,r1,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,f1)
        endif
c---------------------------
c
      call intSG(NORM,natom,r1,Nuc,M,Md,ncont,nshell,c,a,RMM,f1)
      E=E+Exc
      write(*,480) E
c
c nuclear velocities updates
c
c
      ENK=0.D0
      do 80 i=1,ntom
       v1(i,1)=v(i,1)-dth*(f(i,1)+f1(i,1))/Pm(i)
       v1(i,2)=v(i,2)-dth*(f(i,2)+f1(i,2))/Pm(i)
       v1(i,3)=v(i,3)-dth*(f(i,3)+f1(i,3))/Pm(i)
 80   ENK=ENK+(v1(i,1)**2+v1(i,2)**2+v1(i,3)**2)*Pm(i)
      ENK=1822.7926*ENK/2.D0
      T=63158.2*ENK/Nlib
c
       fac=1.0D0
       fac0=1.D0
c test in order to scale correctly (E conservation) ----------------
c      Edif=Etot-E
c     if (mod(istep,Nscale).eq.0.and.Edif.gt.1.0D-06) then
c      fac0=(Etot-E)/ENK
c      fac=sqrt(fac0)
c      endif
c--------------------------------------------------------------------
        ENK=ENK*fac0
c     <<move "new" arrays to "old" arrays for next step>>
c
       ss=0.D0
       do 95 i=1,ntom
        r(i,1)=r1(i,1)
        r(i,2)=r1(i,2)
        r(i,3)=r1(i,3)
        v(i,1)=fac*v1(i,1)
        v(i,2)=fac*v1(i,2)
        v(i,3)=fac*v1(i,3)
        f(i,1)=f1(i,1)
        f(i,2)=f1(i,2)
        f(i,3)=f1(i,3)
c       write(*,*) i,f(i,1),f(i,2),f(i,3)
 95     ss=ss+f(i,1)**2+f(i,2)**2+f(i,3)**2
c
         ss=sqrt(ss)
c
c       write(*,450) istep,E,ENK+E
c       write(19,450) istep,E,ss,ENK
c       write(*,450) istep,E,E+ENK
        write(25,450) istep,E,E+ENK,T
c
c also writes down geometry
c
c


      do 117 i=1,ntom
       write(19) r(i,1),r(i,2),r(i,3)
  117  continue
c
      do 127 i=1,ntom
       write(19) v(i,1),v(i,2),v(i,3)
  127  continue
c
c       write(12,*) istep,E+ENK
c  
 50     continue
c
c Outputs final  MO ---------------------
      if (OPEN) then
c alpha
      kk=M18-1
      do 221 n=1,NCOa
      do 221 l=1,M
       kk=kk+1
 221   X(index(l),M+n)=RMM(kk)
c
      do 226 l=1,M
 226   write(2,400) (X(l,M+n),n=1,NCOa)
c
      kk=M18b-1
      do 222 n=1,NCOb
      do 222 l=1,M
       kk=kk+1
 222   X(index(l),M+n)=RMM(kk)
c
      do 227 l=1,M
 227   write(2,400) (X(l,M+n),n=1,NCOb)
c
      else
      do 220 l=1,M
      do 220 n=1,NCO
 220   X(index(l),M+n)=X(l,M2+n)
c
      do 225 l=1,M
 225   write(2,400) (X(l,M+n),n=1,NCO)
c
      endif
c
c  -------------------------                                            
c also writes down geometry
c
c
      do 17 i=1,ntom
       write(2,500) Iz(i),r(i,1),r(i,2),r(i,3)
  17  continue
c
      do 27 i=1,ntom
       write(2,600) v(i,1),v(i,2),v(i,3)
  27  continue
c
c
 400  format(4(E14.7E2,2x))
 450  format(I4,2x,2(F15.7,2x),E16.8E3)
 500  format(i3,2x,F9.6,2x,F9.6,2x,F9.6)
 600  format(F9.6,2x,F9.6,2x,F9.6)
 650  format(F9.6,2x,F9.6,2x,F9.6,I2)
 480  format ('SCF ENERGY = ',F14.7)
 580  format (I3,2x,3(F12.5,2x))
 125  format('# STEP',4x,'POT. ENERGY',6x,'TOT. ENERGY',6x,'TEMP.')
       return
       end
c
