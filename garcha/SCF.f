c SCF subroutine ----------------------------------
c DIRECT VERSION
c Calls all integrals generator subroutines : 1 el integrals,
c 2 el integrals, exchange fitting , so it gets S matrix, F matrix
c and P matrix in lower storage mode ( symmetric matrices)
c
c Dario Estrin, 1992
c---------------------------------------------------
       subroutine SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
c
      implicit real*8 (a-h,o-z)
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write1,just_int3n
      logical OPEN,SVD,SHFT,GRAD,BSSE,integ,field,sol,free
      logical exter,MEMO
      INCLUDE 'param'
      dimension r(nt,3),nshelld(0:3),nshell(0:3),q(ntq)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
c
      dimension Em(ntq+nss),pc(nss),alpha(nss)
      dimension pci(nt,4),ighost(ntq)
c
c
       dimension RMM(*),xi(3),Ef(3)
c
c auxiliars
c X scratch space in matrices
      dimension X(M,3*M),XX(Md,Md)
c
c test
      dimension aux(ng),ds(nt)
c
      COMMON /TABLE/ STR(880,0:21)
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
c     common /HF/ OPEN,NMAX,NCO,ATRHO,VCINP,DIRECT
c    >   ,IDAMP,EXTR,SHFT,SHI,GOLD,told,write1,Nunp 
      common /Nc/ Ndens
      common /Sys/ SVD,iconst
      common /cav/ a0,epsilon,field,exter,Fx,Fy,Fz
c
      common /index/ index(ng)
      common /coef/ af(ngd)
      common /coef2/ B(ngd,2)
      common /bsse/ BSSE,ighost
      common/Ngeom/ ngeo
      common /ENum/ GRAD
      common /propt/ idip,ipop,ispin,icharge,map(ntq)
      common /intg1/ e_(50,3),wang(50),Nr(0:54)
      common /intg2/ e_2(116,3),wang2(116),Nr2(0:54),e3(194,3),
     > wang3(194)
			common /radii/ Rm(0:54)
c
      common /sol1/ Nsol,natsol,alpha,Em,pc,sol,free

      call timer_start('SCF')

      just_int3n = false


c------------------------------------------------------------------
c
c Pointers
c
c chequeo -----------
c
      Ndens=1
c---------------------
      E=0.0D0
      E1=0.0D0
      En=0.0D0
      E2=0.0D0
      idip=1
      ngeo=ngeo+1

      sq2=sqrt(2.D0)
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
c now H
      M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
c Least squares
      M17=M15+MM
c vectors of MO
      M18=M17+MMd
c weights (in case of using option )
      M19=M18+M*NCO
c
* RAM storage of two-electron integrals (if MEMO=T)
      M20 = M19 + natom*50*Nang   
c
      Nel=2*NCO+Nunp
c Initializations/Defaults
c
      write(*,*) ' SCF CALCULATION  '
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
c -------------------------------------------------------------
c H H core, 1 electron matrix elements
c
c
c
      call int1(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En)
c  
c
c -- SOLVENT CASE --------------------------------------
      if (sol) then
      call intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,pc,RMM,E1s)
c prueba----------
c      xi(1)=0.0
c      xi(2)=0.0
c      xi(3)=0.0
c      call efield(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c    >            c,a,RMM,xi,V,Ef)
c      write(*,*) Ef(1),Ef(2),Ef(3)
c end prueba------
      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
      endif
c test ---------------------------------------------------------
c to check if MO are normalized
c       kki=0
c       kkj=0
c       TE=0.0D0
c       do i=1,NCO
c        ss=0.0D0
c        do ii=1,M
c        do jj=1,M
c        
c          kki=M18+M*(i-1)+ii-1
c          kkj=M18+M*(i-1)+jj-1
c          if (ii.ge.jj) then
c          kk=ii+((M2-jj)*(jj-1))/2
c          else
c          kk=jj+((M2-ii)*(ii-1))/2
c          endif
c        ss=ss+RMM(kki)*RMM(kkj)*RMM(M5+kk-1)
c        enddo
c        enddo
c        write(*,*) i,ss
c        TE=TE+ss
c       enddo
c
c       write(*,*) 2.0D0*TE
c       pause
         E1=0.D0
        do 302 k=1,MM
 302     E1=E1+RMM(k)*RMM(M11+k-1)
c
c Diagonalization of S matrix, after this is not needed anymore
c
c ESSL OPTION ------------------------------------------
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
c ESSL OPTION
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c LAPACK OPTION -----------------------------------------
#ifdef pack
c
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
c
       kk=0
c
       do 261 k=1,NCO
       do 261 i=1,M
       kk=kk+1
 261     RMM(M18+kk-1)=X(i,M2+k)
c
      kk=0
      do 330 j=1,M
      do 330 i=j,M
      kk=kk+1
      RMM(kk)=0.D0
c
      if(i.eq.j) then
       ff=2.D0
      else
       ff=4.D0
      endif
c
      do 330 k=1,NCO
       RMM(kk)=RMM(kk)+ff*X(i,M2+k)*X(j,M2+k)
 330  continue
c
c
      endif
c End of Starting guess (No MO , AO known)-------------------------------
c

      call int22(NORM,natom,r,Nucd,M,Md,ncontd,
     >          nshelld,cd,ad,RMM,XX)
c     
**
      if (MEMO) then
         call int3mem(NORM,natom,r,Nuc,M,M20,ncont,nshell,c,a,Nucd,Md,
     >                 ncontd,nshelld,RMM,cd,ad,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
      endif
****        
c---------------------------------------------------------------------
c Now, damping is performed on the density matrix
c The first 4 iterations ( it may be changed, if necessary)
c when the density is evaluated on the grid, the density
c matrix is used ( slower), after that it is calculated
c using the vectors . Since the vectors are not damped,
c only at the end of the SCF, the density matrix and the
c vectors are 'coherent'
c---------------------------------------------------------------
c CASE  SAVING BASIS FUNCTIONS ON DISK
c ONLY IN LEAST SQUARES
      if (.not.integ) then
      if (dens) then
      write(*,*) 'in nwrite'
      call nwrite(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,NCO,
     >           NCO,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)
c
      else
      write(*,*) 'in write'
      call write(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,
     >           NCO,NCO,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)

      endif
      endif
c LEVEL SHIFT CASE, contruction of initial vectors ------------------
c
      if (SHFT) then
c
         if (MEMO) then
            call int3lu(NORM,natom,Iz,r,Nuc,M,M20,ncont,nshell,c,a,Nucd,
     >                  Md,ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
         else  
c
      call int3(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
c
      endif

        do 720 i=1,M
        do 720 j=1,M
         X(i,M+j)=0.D0
         do 722 k=1,j
         X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M5+j+(M2-k)*(k-1)/2-1)
  722     continue
c
         do 723 k=j+1,M
  723     X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M5+k+(M2-j)*(j-1)/2-1)
c
  720     continue
c
         kk=0
        do 730 j=1,M
        do 730 i=j,M
         kk=kk+1
         RMM(M5+kk-1)=0.D0
        do 732 k=1,M
 732      RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
 730     continue
c
c diagonalization now
c
        do 749 i=1,M
         RMM(M15+i-1)=0.
 749      RMM(M13+i-1)=0.
c
c ESSL OPTION --------------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
c
       call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
c-----------------------------------------------------------
c
c new coefficients
c
       do 750 i=1,M
       do 750 j=1,M
c
        X(i,M2+j)=0.D0
       do 752 k=1,M
  752    X(i,M2+j)=X(i,M2+j)+X(i,k)*X(k,M+j)
  750    continue
c
      do 741 i=1,M
      do 741 j=1,M
       X(i,j)=X(i,M2+j)
 741  continue
c
c Density Matrix, Old one ?
      kk=0
      do 830 j=1,M
      do 830 i=j,M
      kk=kk+1
      RMM(M3+kk-1)=0.D0
c
      if(i.eq.j) then
       ff=2.D0
      else
       ff=4.D0
      endif
c
      do 830 k=1,NCO
       RMM(M3+kk-1)=RMM(M3+kk-1)+ff*X(i,M2+k)*X(j,M2+k)
 830  continue
c
      do 840 k=1,MM
 840   RMM(k)=RMM(M3+k-1)
c
      endif

c
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c
      do 999 while (good.ge.told)
c
      if (niter.ge.NMAX) then
       write(*,*) 'NO CONVERGENCE AT ',NMAX,' ITERATIONS'
       goto 995
      endif
c
c prueba----------
c      xi(1)=0.0
c      xi(2)=0.0
c      xi(3)=0.0
c      call efield(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c    >            c,a,RMM,xi,V,Ef)
c      write(*,*) Ef(1),Ef(2),Ef(3)
c end prueba------
c
c int3N is called close to convergence, it uses previous values
c for the fitting coefficient vector af
c
      niter=niter+1
      if (MEMO) then
         call int3lu(NORM,natom,Iz,r,Nuc,M,M20,ncont,nshell,c,a,Nucd,Md,
     >               ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
      else                            
      if ((.not.just_int3n).and.(good.gt.3.0D0*told)) then
c
      call int3(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
c
c------------------
      else
      just_int3n = true
      call int3N(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,XX,E2,Ex,
     >     nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
c     write(*,*) 'INTN'
      endif
      endif
c-------------------------------------------------------
      E1=0.0D0
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
         fac=(2.54D0*2.00D0)
        endif
c
        call dip2(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,g,Fx,Fy,Fz,
     >     nopt,OPEN)
c
        E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/fac -
     >    0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
        endif
c----------------------------------------------------------------
c E1 includes solvent 1 electron contributions
        do 303 k=1,MM
 303     E1=E1+RMM(k)*RMM(M11+k-1)
c
c ---  debugging ------
c      write(*,*) E1+E2
c      return
c ---------------------
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
c
      do 135 k=1,MM
       kk=M5+k-1
       kk2=M3+k-1
 135   RMM(kk)=(RMM(kk)+DAMP*RMM(kk2))/(1.D0+DAMP)
c
c the newly constructed damped matrix is stored, for next iteration
c in RMM(M3)
c
      do 137 k=1,MM
       kk=M5+k-1
       kk2=M3+k-1
 137   RMM(kk2)=RMM(kk)
c
c
        do 20 i=1,M
        do 20 j=1,M
         X(i,M+j)=0.D0
         do 22 k=1,j
         X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M5+j+(M2-k)*(k-1)/2-1)
  22     continue
c
         do 23 k=j+1,M
  23     X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M5+k+(M2-j)*(j-1)/2-1)
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
        do 49 i=1,M
         RMM(M15+i-1)=0.D0
 49      RMM(M13+i-1)=0.D0
c
c---- LEVEL SHIFT 
c
       if (SHFT) then
c       shi=shi*0.99
c adition of level shifts
c constant to diagonal (virtual) elements
       do i=NCO+1,M
        ii=i+(i-1)*(M2-i)/2
        RMM(M5+ii-1)=RMM(M5+ii-1)+shi
       enddo
       endif
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
c
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
c --- For the first iteration, damping on density matrix 
c Important for the case of strating guess of AO
c

c
       kk=0
       do 61 k=1,NCO
       do 61 i=1,M
        kk=kk+1
 61     RMM(M18+kk-1)=X(i,M2+k)
c
c Construction of new density matrix and comparison with old one
      kk=0
      good=0.0D0
      do 130 j=1,M
      do 130 i=j,M
      kk=kk+1
      tmp=RMM(kk)
      RMM(kk)=0.D0
c
      if(i.eq.j) then
       ff=2.D0
      else
       ff=4.D0
      endif
c
      do 139 k=1,NCO
       RMM(kk)=RMM(kk)+ff*X(i,M2+k)*X(j,M2+k)
 139  continue
        del=RMM(kk)-tmp
        if (i.ne.j) then
         del=del*sq2
        endif
        good=good+del**2
 130  continue
c
       good=sqrt(good)/float(M)
      if (SHFT) then
c Level Shifting
      do 141 i=1,M
      do 141 j=1,M
       X(i,j)=X(i,M2+j)
 141  continue
      endif
c
c
      if (SHFT) then
      DAMP=0.0D0
      endif
c
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
c
       if (nopt.ne.3) then
       write(*,300) niter,DAMP,E
       endif
c
c      write(*,*) 'Coulomb E',E2-Ex,Ex
       if (write1) then
c
      open(unit=3,file='restart')
c outputs final  MO ---------------------
      do 320 l=1,M
      do 320 n=1,NCO
 320   X(index(l),M+n)=X(l,M2+n)
c
       write(3,*) niter,E
c     
      do 325 l=1,M
 325   write(3,400) (X(l,M+n),n=1,NCO)
c
      close(3)
      endif
c
 999   continue
 995   continue
c
c -- SOLVENT CASE --------------------------------------
      if (sol) then
      call intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,pc,RMM,E1s)
      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
      Es=Es+E1s
      endif
c--------------------------------------------------------------
       if (GRAD) then
        if(nopt.eq.0) then
         write(*,*)
         write(*,600)
         write(*,610)
         write(*,620) E1,E2-Ex,En
         if (sol) then
          write(*,615)
          write(*,625) Es
         endif
       endif
			 write(*,*) 'ultimo paso SCF'
       call timer_start('exchnum')
#ifdef GPU
#ifdef ULTIMA_CPU
       call exchnum(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >              M18,NCO,Exc,nopt)
#else
			 if (igrid.ne.igrid2) then
       call gpu_new_grid(igrid)
			 endif
       call gpu_solve_groups(1, Exc, 0)
#endif
#else
#ifdef ULTIMA_GPU
			 if (igrid.ne.igrid2) then
       call gpu_new_grid(igrid)
			 endif
       call gpu_solve_groups(1, Exc, 0)
#else      
      call exchnum(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >              M18,NCO,Exc,nopt)
			write(*,*) 'total final:',Exc
#endif       
#endif
       call timer_stop('exchnum')      
       E=E+Exc-Ex
c
c--------------------------------------------------------------
       write(*,*)
       write(*,450) E
       else
       E=E-Ex
       endif
c calculation of energy weighted density matrix
c
      kk=0
      do 307 j=1,M
      do 307 i=j,M
      kk=kk+1
      RMM(M15+kk-1)=0.D0
c
      if(i.eq.j) then
       ff=2.D0
      else
       ff=4.D0
      endif
c
      do 309 k=1,NCO
       RMM(M15+kk-1)=RMM(M15+kk-1)-RMM(M13+k-1)*ff*X(i,M2+k)*X(j,M2+k)
 309  continue
c
 307   continue
c
c
c
      if (nopt.eq.0) then
c
c PROPERTIES CALCULATION
c calculates dipole moment
c
       if (idip.eq.1) then
        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >       Nel,ux,uy,uz)
      u=sqrt(ux**2+uy**2+uz**2)
c
      write(*,*)
      write(*,*) 'DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
      write(*,900) ux,uy,uz,u
      write(*,*)
c u in Debyes
       endif
c
c calculates Mulliken poputations
       if (ipop.eq.1) then
        call int1(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En)
c
c--------------------------------------------------------------
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
         write(*,*) 'MULLIKEN POPULATION ANALYSIS'
         write(*,770)

        do n=1,natom
         write(*,760) n,Iz(n),q(n)
        enddo
c
        write(*,*)
        endif
c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
c
        if (icharge.eq.1) then
          Q1=-(2*NCO+Nunp)
         do n=1,natom
          Q1=Q1+Iz(n)
         enddo
         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,RMM,map,Q1)
        endif
c
c--------------------------------------------------------------
c outputs final  MO ---------------------
      do 420 l=1,M
c     do 420 n=1,NCO+3
      do 420 n=1,M
 420   X(index(l),M+n)=X(l,M2+n)
c
      do 225 l=1,M
 225   write(2,400) (X(l,M+n),n=1,NCO)
c-------------------------------------------------
c writes down MO coefficients and orbital energies
        write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
       do n=1,NCO
        write(29,850) n,RMM(M13+n-1)
        write(29,400) (X(l,M+n),l=1,M)
       enddo
       do n=NCO+1,M
        write(29,851) n,RMM(M13+n-1)
        write(29,400) (X(l,M+n),l=1,M)
       enddo
       close(29)
c
c-------------------------------------------------
      endif
c
c

 500  format('SCF TIME ',I6,' sec')
 450  format ('SCF ENERGY = ',F19.12)
 400  format(4(E14.7E2,2x))
 300  format(I3,E14.6,2x,F14.7)
 600  format('  ENERGY CONTRIBUTIONS IN A.U.')
 610  format(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
 615  format(2x,'SOLVENT')
 620  format(F14.7,4x,F14.7,4x,F14.7)
 625  format(F14.7)
 760  format(I3,9x,I3,6x,F10.4)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
 850  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
 851  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7,
     >    '(NON OCC.)')
 900  format(3(F10.4,2x),2x,F10.4)
 777  format(4(F8.4,2x))
 776  format (3(F8.4,2x))
c---- DEBUGGINGS
c      write(*,*) 'Exc, integrated and calculated',Exc,Ex
c      write(*,*) 'Coulomb energy',E2-Ex
c
      call timer_stop('SCF')
       return
       end
C  -------------------------                                            
