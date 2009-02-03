c SCF subroutine ----------------------------------
c DIRECT VERSION
c calls all integrals generator subroutines : 1 el integrals,
c 2 el integrals, exchange fitting , so it gets S matrix, F matrix
c and P matrix in lower storage mode ( symmetric matrices)
c
c Dario Estrin, 1992
c---------------------------------------------------
       subroutine SCFOP(MEMO,NORM,natom,Iz,r,Nuc,M,ncont
     >    ,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     > nopt,pc,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1,
     > FQQ,Q,IT,ITEL,NIN,IPR1,E1s,EAC,
     > ux,uy,uz,NPAS)
c
      implicit real*8 (a-h,o-z)
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write1
      logical OPEN,SVD,SHFT,GRAD,BSSE,integ,field,sol,free
      logical exter,MEMO
      integer nopt,iconst,igrid,igrid2
      INCLUDE 'param'
      dimension r(nt,3),nshelld(0:3),nshell(0:3),q(ntq)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(ngd),ncontd(ngd)
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng),Iz(nt)

      dimension Em(ntq+nss),Rm(ntq+nss),pc(nt),alpha(nss)
      dimension pci(nt,4),ighost(ntq)

      dimension xi(3),Ef(3)

c auxiliars
c X scratch space in matrices
c      dimension X(ngDyn,ng3),XX(ngdDyn,ngdDyn),RMM(ng2)
      dimension X(M,3*M),XX(Md,Md),RMM(*)
      DIMENSION FQQ(NATOM+NSOL*NATSOL),EAC(NATOM+NSOL*NATSOL) 

      COMMON /TABLE/ STR(880,0:21)
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
c     common /HF/ nopt,OPEN,NMAX,NCO,ATRHO,VCINP,DIRECT,
c    >             IDAMP,EXTR,SHFT,SHI,GOLD,told,write1,Nunp
      common /Nc/ Ndens
      common /Sys/ SVD,iconst
      common /cav/ a0,epsilon,field,exter,Fx,Fy,Fz

c
      common /index/ index(ng)
      common /coef/ af(ngd)
      common /coef2/ B(ngd,3)
      common /bsse/ BSSE,ighost
      common/Ngeom/ ngeo
      common /ENum/ GRAD
      common /propt/ idip,ipop,ispin,icharge,map(ntq)

      common /sol1/ Nsol,natsol,alpha,Em,Rm,sol,free
c------------------------------------------------------------------
c
c Pointers
c

      E=0.0D0
      E1=0.0D0
      En=0.0D0
      E2=0.0D0
      Es=0.0D0
      idip=1
      ngeo=ngeo+1


      sq2=sqrt(2.D0)
      MM=M*(M+1)/2 
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M
      M2a=3*M
c
      NCOa=NCO
      NCOb=NCO+Nunp
      Nel=2*NCO+Nunp

c      write(*,*)'SCFop: write1',write1

c------------------------------------------------
c
c first P
      M1=1
c now F alpha
      M3=M1+MM
c now S, F beta also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues closed shell and alpha spin  in open shell case)
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
c Least squares
      M17=M15+MM
c vectors of MO alpha
      M18=M17+MMd
c vectors of MO beta
      M18b=M18+M*NCOa
c weights (in case of using option )
      M19=M18b+M*NCOb
c new Fock matrix alpha
      M20=M19+natom*50*Nang
c new Fock matrix beta
      M21=M20+MM
c eigenvalues (beta spin in open shell case_
      M22=M21+MM
c
* RAM storage of two-electron integrals (if MEMO=T)
      M23 = M22 +  M  
c
c------------------------------------------------
c Initializations/Defaults
c
      IF(MOD((IT-NIN),IPR1).EQ.0)THEN
      write(*,*) ' SCF OPEN SHELL CALCULATION  ',itel
      ENDIF
c
c
      good=1.00D0
      niter=0
      D1=1.D0
      D2=1.D0
      DAMP0=GOLD
      DAMP=DAMP0
     
c
       Qc=0.0D0
      do 99 i=1,natom
       Qc=Qc+Iz(i)
   99   continue
       Qc=Qc-Nel
      Qc2=Qc**2
c
c H H core, 1 electron matrix elements
c
c
      call int1(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En)
c
c -- SOLVENT CASE --------------------------------------
      if (sol) then
      call intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,E1s,FQQ,IT,ITEL,NIN,
     >            IPR1,EAC,NPAS)


c      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
      endif
c
         E1=0.D0
        do 302 k=1,MM
  302     E1=E1+RMM(k)*RMM(M11+k-1)
c
c
c------------------------------------------------------
c Diagonalization of S matrix, after this is not needed anymore
c
c -----
c ESSL OPTION ----------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X,M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
       call dspev('V','L',M,RMM(M5),RMM(M13),X,M,RMM(M15),info)
#endif
c-----------------------------------------------------------
c 
c X transformation matrix , canonical orthogonalization
c LINEAR DEPENDENCY ELIMINATION
c
        do 10 j=1,M
          if (RMM(M13+j-1).lt.1.0D-06) then
          write(*,*) 'LINEAR DEPENDENCY DETECTED'
         do 11 i=1,M
   11      X(i,j)=0.0D0
          else
         do 12 i=1,M
   12      X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))
         endif
   10      continue
c
          do i=1,M
          do j=1,M
           X(i,M2a+j)=X(i,j)
          enddo
          enddo
c CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
c
        if ((.not.ATRHO).and.(.not.VCINP).and.ngeo.eq.1) then
        do 220 i=1,M
        do 220 j=1,M
         X(i,M+j)=0.D0
         do 222 k=1,j
         X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+j+(M2-k)*(k-1)/2-1)
  222     continue
c
         do 223 k=j+1,M
  223     X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+k+(M2-j)*(j-1)/2-1)
c
  220     continue
c
         kk=0
        do 230 j=1,M
        do 230 i=j,M
         kk=kk+1
         RMM(M5+kk-1)=0.D0
        do 232 k=1,M
  232      RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
  230     continue
c
c diagonalization now
c
        do 249 i=1,M
         RMM(M15+i-1)=0.D0
  249      RMM(M13+i-1)=0.D0
c
c ESSL OPTION ---------------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
       call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
c-----------------------------------------------------------
       do 250 i=1,M
       do 250 j=1,M
c
        X(i,M2+j)=0.D0
       do 252 k=1,M
  252    X(i,M2+j)=X(i,M2+j)+X(i,k)*X(k,M+j)
  250    continue
c
c Density Matrix
c alpha and beta coefficients set equal
       kk=0
c
       do 261 k=1,NCOa
       do 261 i=1,M
       kk=kk+1
  261     RMM(M18+kk-1)=X(i,M2+k)
c
       kk=0
       do 262 k=1,NCOb
       do 262 i=1,M
       kk=kk+1
  262     RMM(M18b+kk-1)=X(i,M2+k)
c
      kk=0
      do 330 j=1,M
      do 330 i=j,M
      kk=kk+1
      RMM(kk)=0.D0
c
      do 337 k=1,NCOa
       RMM(kk)=RMM(kk)+X(i,M2+k)*X(j,M2+k)
  337  continue
c
      do 331 k=1,NCOb
       RMM(kk)=RMM(kk)+X(i,M2+k)*X(j,M2+k)
  331  continue
c
       if (i.ne.j) then
       RMM(kk)=2.0D0*RMM(kk)
       endif
c
  330  continue
c
      endif
c End of Starting guess (No MO , AO known)-------------------------------
c

      call intt2(NORM,natom,r,Nucd,M,Md,ncontd,
     >          nshelld,cd,ad,RMM,XX)
**
      if (MEMO) then
         call int3mem(NORM,natom,r,Nuc,M,M23,ncont,nshell,c,a,Nucd,Md,
     >                 ncontd,nshelld,RMM,cd,ad,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
      endif
****       
c---------------------------------------------------------------------
c CASE  SAVING BASIS FUNCTIONS ON DISK
c ONLY IN LEAST SQUARES
      if (.not.integ) then
      if (dens) then
      write(*,*) 'in nwrite'
      call nwrite(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,NCOa,
     >           NCOb,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)
c
      else
      write(*,*) 'in write'
      call write(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,
     >    NCOa,NCOb,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)

      endif
      endif
c-------------------------------------------------------------------
c
      do 999 while (good.ge.told)
c
c
      if (niter.ge.NMAX) then
       write(*,*) 'NO CONVERGENCE AT ',NMAX,' ITERATIONS'
       goto 995
      endif
c
c
c int3N is called close to convergence, it uses previous values
c for the fitting coefficient vector af
c
      niter=niter+1
      if (MEMO) then
         call int3lu(NORM,natom,Iz,r,Nuc,M,M23,ncont,nshell,c,a,Nucd,Md,
     >               ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
      else         

      if (good.gt.3.0D0*told) then
      call int3(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
c
c------------------
      else
      call int3N(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
c     write(*,*) 'INTN'
      endif
      endif
c
c-------------------------------------------------------
c
      E1=0.0D0
c
c REACTION FIELD  --------------------------------------------
c
c       if (field) then
c       call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
c    >     Nel,ux,uy,uz)
c       ux=ux/2.54
c       uy=uy/2.54
c       uz=uz/2.54
c       g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
c       call dip2(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,g,ux,uy,uz)
c
c       E1=-0.50D0*g*(ux**2+uy**2+uz**2)
c       endif
c
c REACTION FIELD CASE --------------------------------------------
c
        if (field) then
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,Nel,
     >       ux,uy,uz)
c
        if (exter) then
         g=1.0D0
         fac=2.54D0
        else
         g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
         Fx=ux/2.54D0
         Fy=uy/2.54D0
         Fz=uz/2.54D0
         fac=(2.00D0*2.54D0)
        endif
c
        call dip2(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,g,Fx,Fy,Fz,
     >    nopt,OPEN)
c
        E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/fac-
     >    0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
        endif
c----------------------------------------------------------------
c E1 includes the solvent electrostatic matrix elements (sol T)
c
        do 303 k=1,MM
  303     E1=E1+RMM(k)*RMM(M11+k-1)
c
c ---  debugging ------
c      write(*,*) E1+E2
c      return
c ---------------------
c
c
c now, we know S matrix, and F matrix, and E for a given P
c 1) diagonalize S, get X=U s^(-1/2)
c 2) get U F Ut
c 3) diagonalize F
c 4) get vec ( coeff) ---->  P new
c 5) iterate
c call diagonalization routine for S , get after U s^(-1/2)
c where U matrix with eigenvectors of S , and s is vector with
c eigenvalues
c
c here in RMM(M5) it is stored the new Fock matrix
c test damping on Fock matrix
c
c      DAMP=gold
      if (niter.eq.1) then
      DAMP=0.0D0
      endif
c
c Damping Alpha Matrix
c
      do 135 k=1,MM
       kk=M5+k-1
       kk2=M20+k-1
  135   RMM(kk)=(RMM(kk)+DAMP*RMM(kk2))/(1.D0+DAMP)
c
c Damping  Beta Matrix
c
      do 136 k=1,MM
       kk=M3+k-1
       kk2=M21+k-1
  136   RMM(kk)=(RMM(kk)+DAMP*RMM(kk2))/(1.D0+DAMP)
c
c the newly constructed damped matrix is stored, for next iteration
c in RMM(M20) and RMM(M21) for alpha and beta spin respectively.
c
      do 137 k=1,MM
       kk=M5+k-1
       kk2=M20+k-1
       RMM(kk2)=RMM(kk)
c
       kk=M3+k-1
       kk2=M21+k-1
  137   RMM(kk2)=RMM(kk)
c
c alpha case -------
c
        do 20 i=1,M
        do 20 j=1,M
         X(i,M+j)=0.D0
         do 25 k=1,j
         X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M5+j+(M2-k)*(k-1)/2-1)
   25     continue
c
         do 26 k=j+1,M
   26     X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M5+k+(M2-j)*(j-1)/2-1)
c
   20     continue
c
         kk=0
        do 30 j=1,M
        do 30 i=j,M
         kk=kk+1
         RMM(M5+kk-1)=0.D0
        do 32 k=1,M
   32      RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
   30     continue
c
c now F contains transformed F
c
c diagonalization now
c
        do 47 i=1,M
         RMM(M15+i-1)=0.D0
   47      RMM(M13+i-1)=0.D0
c
        if (SHFT) then
        if (niter.ge.2) then
c adition of level shifts
c constant to diagonal (virtual) elements
       do i=NCOa+1,M
        ii=i+(i-1)*(M2-i)/2
        RMM(M5+ii-1)=RMM(M5+ii-1)+shi
       enddo
       endif
       endif
c
c ESSL OPTION ------------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
c
#ifdef pack
       call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
c-----------------------------------------------------------
c diagonalization now
c
c
c new coefficients
c
       do 50 i=1,M
       do 50 j=1,M
c
        X(i,M2+j)=0.D0
       do 52 k=1,M
   52    X(i,M2+j)=X(i,M2+j)+X(i,k)*X(k,M+j)
   50    continue
c
c
       kk=0
       do 61 k=1,NCOa
       do 61 i=1,M
        kk=kk+1
   61     RMM(M18+kk-1)=X(i,M2+k)
c
c xxxxxxx aca poner que escriba ------------------------------
      if ((good.le.5.0D0*told.and.nopt.eq.0).or.niter.eq.nmax) then
      do 4200 l=1,M
      do 4200 n=1,M
 4200  X(index(l),M+n)=X(l,M2+n)

        write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, SPIN ALPHA'
c this option in order to look orbitals
c      do n=1,NCO+5
c this other if one is interested in fargment  populations
       do n=1,M
        write(29,850) n,RMM(M13+n-1)
        write(29,400) (X(l,M+n),l=1,M)
       enddo
      endif
c--------------------------------------------------------------
       if (SHFT) then
c      if (niter.ge.2) then
c Level Shifting
      do 141 i=1,M
      do 141 j=1,M
       X(i,j)=X(i,M2+j)
  141  continue
c     endif
      endif
c
c Beta case -------------------
c
c
        do 21 i=1,M
        do 21 j=1,M
         X(i,M+j)=0.D0
         do 22 k=1,j
         X(i,M+j)=X(i,M+j)+X(k,M2a+i)*RMM(M3+j+(M2-k)*(k-1)/2-1)
   22     continue
c
         do 23 k=j+1,M
   23     X(i,M+j)=X(i,M+j)+X(k,M2a+i)*RMM(M3+k+(M2-j)*(j-1)/2-1)
c
   21     continue
c
         kk=0
        do 35 j=1,M
        do 35 i=j,M
         kk=kk+1
         RMM(M3+kk-1)=0.D0
        do 37 k=1,M
   37      RMM(M3+kk-1)=RMM(M3+kk-1)+X(i,M+k)*X(k,M2a+j)
   35     continue
c now F contains transformed F
c
c diagonalization now
c
        do 49 i=1,M
         RMM(M15+i-1)=0.D0
   49      RMM(M22+i-1)=0.D0
c
        if (SHFT) then
        if (niter.ge.2) then
c adition of level shifts
c constant to diagonal (virtual) elements
       do i=NCOb+1,M
        ii=i+(i-1)*(M2-i)/2
        RMM(M3+ii-1)=RMM(M3+ii-1)+shi
       enddo
       endif
       endif
c
c
c ESSL OPTION------------------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M3),RMM(M22),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
       call dspev('V','L',M,RMM(M3),RMM(M22),X(1,M+1),M,RMM(M15),info)
#endif
c-----------------------------------------------------------
c diagonalization now
c
c
c new coefficients
c
       do 51 i=1,M
       do 51 j=1,M
c
        X(i,M2+j)=0.D0
       do 53 k=1,M
   53    X(i,M2+j)=X(i,M2+j)+X(i,M2a+k)*X(k,M+j)
   51    continue
c
c xxxxxxx  aca poner que escriba -------------------------------
      if ((good.le.5.0D0*told.and.nopt.eq.0).or.niter.eq.nmax) then
      do 4201 l=1,M
      do 4201 n=1,M
 4201  X(index(l),M+n)=X(l,M2+n)
c
        write(29,*)
        write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, SPIN BETA'
c this option is for looking orbitals
c      do n=1,NCO+5
c this other for fragment populations
       do n=1,M
        write(29,850) n,RMM(M22+n-1)
        write(29,400) (X(l,M+n),l=1,M)
       enddo
        rewind(29)
       endif
c------------------------------------------------------------
c
       kk=0
       do 67 k=1,NCOb
       do 67 i=1,M
        kk=kk+1
   67     RMM(M18b+kk-1)=X(i,M2+k)
c
      if (SHFT) then
c     if (niter.ge.2) then
c Level Shifting
      do 142 i=1,M
      do 142 j=1,M
       X(i,M2a+j)=X(i,M2+j)
  142  continue
      endif
c     endif
c
c-----------------------------------------
c Construction of new density matrix and comparison with old one
c
      kk=0
      good=0.0D0
      do 131 j=1,M
      do 131 i=j,M
      kk=kk+1
      tmp=RMM(kk)
      RMM(kk)=0.D0
c
      do 139 k=1,NCOa
       k0=M18+M*(k-1)-1
       ki=k0+i
       kj=k0+j
       RMM(kk)=RMM(kk)+RMM(ki)*RMM(kj)
c
  139  continue
c
      do 149 k=1,NCOb
       k0=M18b+M*(k-1)-1
       ki=k0+i
       kj=k0+j
       RMM(kk)=RMM(kk)+RMM(ki)*RMM(kj)
c
  149  continue
c
       if (i.ne.j) then
       RMM(kk)=2.0D0*RMM(kk)
       endif
c
        del=RMM(kk)-tmp
        if (i.ne.j) then
         del=del*sq2
        endif
        good=good+del**2
  131  continue
c
       good=sqrt(good)/float(M)
c
c--- Damping factor update - 
        DAMP=DAMP0
        if (IDAMP.EQ.1) then
        DAMP=DAMP0
        if (abs(D1).lt.1.D-5) then
         fac=dmax1(0.90D0,abs(D1/D2))
         fac=dmin1(fac,1.1D0)
        DAMP=DAMP0*fac
        endif
c
       E=E1+E2+En
c
       if (sol) then
        E=E+Es
       endif
c
        D2=D1
        D1=(E-E0)
c
        E0=E
        DAMP0=DAMP
        endif
c
       E=E1+E2+En
       if (sol) then
        E=E+Es
       endif

c
c      if (nopt.le.2) then
c       write(*,300) niter,DAMP,E
c      endif
c
c      write(*,*) 'Coulomb E',E2-Ex,Ex
       if (write1) then
c
      open(unit=3,file='restart')
c outputs final  MO ---------------------
c alpha
      do 320 l=1,M
      do 320 n=1,NCOa
       kk=M18+(l-1)+M*(n-1)
  320   X(index(l),M+n)=RMM(kk)
c
      do 325 l=1,M
  325   write(3,400) (X(l,M+n),n=1,NCOa)
c
c beta
      do 323 l=1,M
      do 323 n=1,NCOb
       kk=M18b+(l-1)+M*(n-1)
  323   X(index(l),M+n)=RMM(kk)
c
      do 327 l=1,M
  327   write(3,400) (X(l,M+n),n=1,NCOb)
c
       write(3,*) niter,E
      close(3)
      endif
c
  999   continue
  995   continue
c
c -- SOLVENT CASE --------------------------------------
      if (sol) then
      call intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,ncont,nshell,
     > c,a,pc,RMM,E1s,FQQ,IT,ITEL,NIN,IPR1,EAC,NPAS)
    

c      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)

c--- E1s NO esta sumado aca al Esolvente!!(Esta sumado en main)
c      Es=Es+E1s
      endif


c---------------------------------------------------------
       if (GRAD) then

       if (nopt.eq.0) then
       IF(MOD((IT-NIN),IPR1).EQ.0)THEN
       write(*,*)
       write(*,600)
       write(*,610)
       write(*,620) E1,E2-Ex,En
       ENDIF
       endif
c         if (sol) then
c          write(*,615)
c          write(*,625) Es
c         endif
        call exchnumop(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >                M18,NCOa,NCOb,Exc,nopt,IT,ITEL,NIN,IPR1)

       E=E+Exc-Ex

c---------------------------------------------------------

          IF(MOD((IT-NIN),IPR1).EQ.0)THEN
       write(*,450) E
          ENDIF

       else
       E=E-Ex
       endif
c
c calculation of energy weighted density matrix
c
c test -----------
c     write(*,*) 'Eig alpha'
c     do i=1,NCOa+10
c      write(*,*) i,RMM(M13+i-1)
c     enddo
c
c     write(*,*) 'Eig beta'
c     do i=1,NCOb+10
c      write(*,*) i,RMM(M22+i-1)
c     enddo
c
      kk=M15-1
      do 307 j=1,M
      do 307 i=j,M
      kk=kk+1
      RMM(kk)=0.D0
c
      do 309 k=1,NCOa
       k0=M18+M*(k-1)-1
       ki=k0+i
       kj=k0+j
       RMM(kk)=RMM(kk)-RMM(M13+k-1)*RMM(ki)*RMM(kj)
  309  continue
c
      do 409 k=1,NCOb
       k0=M18b+M*(k-1)-1
       ki=k0+i
       kj=k0+j
       RMM(kk)=RMM(kk)-RMM(M22+k-1)*RMM(ki)*RMM(kj)
  409  continue
c
       if (i.ne.j) then
        RMM(kk)=2.0D0*RMM(kk)
       endif
c
  307   continue
c
c-----------------------------------------------------------------
c PROPERTIES CALCULATION 
c calculates dipole moment
c
       if (idip.eq.1) then
      call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     > Nel,ux,uy,uz)
      u=sqrt(ux**2+uy**2+uz**2)

      IF(MOD((IT-NIN),IPR1).EQ.0)THEN
      write(*,*)
      write(*,*) 'DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
      write(*,900) ux,uy,uz,u
      write(*,*)
      ENDIF

c u in Debyes
c
       endif

c Calculates Mulliken Populations
       if (ipop.eq.1) then
        call int1(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En)
c
c
        do n=1,natom
         q(n)=Iz(n)
        enddo
c
        do i=1,M
c
        do j=1,i-1
         kk=i+(M2-j)*(j-1)/2
         t0=RMM(kk)*RMM(M5+kk-1)/2.D0
         q(Nuc(i))=q(Nuc(i))-t0
        enddo
c
         kk=i+(M2-i)*(i-1)/2
         t0=RMM(kk)*RMM(M5+kk-1)
         q(Nuc(i))=q(Nuc(i))-t0
c
         do j=i+1,M
         kk=j+(M2-i)*(i-1)/2
         t0=RMM(kk)*RMM(M5+kk-1)/2.D0
         q(Nuc(i))=q(Nuc(i))-t0
         enddo
        enddo
c

         IF(MOD((IT-NIN),IPR1).EQ.0)THEN
         write(*,*) 'MULLIKEN POPULATION ANALYSIS  STEP No ',ITEL
         write(*,770)

        do n=1,natom
         write(*,760) n,Iz(n),q(n)
        enddo
c
         write(*,*)
         ENDIF
c
c UNPAIRED SPIN POPULATION
c M7 spin alpha density matrix, M11 spin beta density matrix
      kk=M7-1
      kk1=M11-1
      do 507 j=1,M
      do 507 i=j,M
      kk=kk+1
      kk1=kk1+1
      RMM(kk)=0.0D0
      RMM(kk1)=0.0D0
c
      do 509 k=1,NCOa
       k0=M18+M*(k-1)-1
       ki=k0+i
       kj=k0+j
       RMM(kk)=RMM(kk) + RMM(ki)*RMM(kj)
  509  continue
c
      do 519 k=1,NCOb
       k0=M18b+M*(k-1)-1
       ki=k0+i
       kj=k0+j
       RMM(kk1)=RMM(kk1)  + RMM(ki)*RMM(kj)
  519  continue
c
       if (i.ne.j) then
        RMM(kk)=2.0D0*RMM(kk)
        RMM(kk1)=2.0D0*RMM(kk1)
       endif
c
  507   continue
       do n=1,natom
         q(n)=0.0D0
        enddo
c
        do i=1,M
c
        do j=1,i-1
         kk=i+(M2-j)*(j-1)/2
         kka=kk+M7-1
         kkb=kk+M11-1
         
         t0=(RMM(kka)-RMM(kkb))*RMM(M5+kk-1)/2.D0
         q(Nuc(i))=q(Nuc(i))-t0
        enddo
c
         kk=i+(M2-i)*(i-1)/2
         kka=kk+M7-1
         kkb=kk+M11-1
         t0=(RMM(kka)-RMM(kkb))*RMM(M5+kk-1)
         q(Nuc(i))=q(Nuc(i))-t0
c
         do j=i+1,M
         kk=j+(M2-i)*(i-1)/2
         kka=kk+M7-1
         kkb=kk+M11-1
c
         t0=(RMM(kka)-RMM(kkb))*RMM(M5+kk-1)/2.D0
         q(Nuc(i))=q(Nuc(i))-t0
         enddo
        enddo
c
        IF(MOD((IT-NIN),IPR1).EQ.0)THEN
         write(*,*) 'UNPAIRED SPIN MULLIKEN POPULATION ANALYSIS ',ITEL
         write(*,770)

        do n=1,natom
         write(*,760) n,Iz(n),q(n)
        enddo
        ENDIF

        endif
c
c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
c
        if (icharge.eq.1) then
          Q1=-(2*NCO+Nunp)
         do n=1,natom
          Q1=Q1+Iz(n)
         enddo
c         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c     *                 c,a,RMM,map,Q1)
        endif
c-----------------------------------------------------
      do 420 l=1,M
      do 420 n=1,NCOa
       kk=M18+(l-1)+M*(n-1)
  420   X(index(l),M+n)=RMM(kk)
c
      do 225 l=1,M
  225   write(2,400) (X(l,M+n),n=1,NCOa)
c
c-------------------------------------------------
       do 520 l=1,M
       do 520 n=1,NCOb
       kk=M18b+(l-1)+M*(n-1)
  520   X(index(l),M+n)=RMM(kk)
c
      do 425 l=1,M
  425   write(2,400) (X(l,M+n),n=1,NCOb)
c-------------------------------------------------
c
c

  500  format('SCF TIME ',I6,' sec')
  450  format ('SCF ENERGY = ',F14.7)
  400  format(4(E14.7E2,2x))
  600  format('  ENERGY CONTRIBUTIONS IN A.U.')
  610  format(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
  615  format(2x,'SOLVENT')
  620  format(F14.7,4x,F14.7,4x,F14.7)
  625  format(F14.7)
  760  format(I3,9x,I3,6x,F10.4)
  770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
  300  format(I3,E14.6,2x,F14.7)
  850  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
  900  format(3(F10.4,2x),2x,F10.4)
c---- DEBUGGINGS
c      write(*,*) 'Exc, integrated and calculated',Exc,Ex
c      write(*,*) 'Coulomb energy',E2-Ex
c
       return
       end
C  -------------------------                                            
