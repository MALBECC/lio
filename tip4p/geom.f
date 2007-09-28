c Geometry Optimization Subroutine ----------------------------------
c Vibrational Analysis
c Thermodynamical Properties
c
c calls SCF subroutine for calculating energy at each geometry ,
c then calls int1G int3G and intSG, that calculate gradients of
c 1 electron part, 2 electron part and overlap respectively
c After that, a move is done, and the process is continued
c
c all specifications given in the namelist SCF are necessary
c
c FTOL, tolerance for the difference in the results of 2
c consecutive line minimizations
c
c imin 0 DFP method with line searches
c imin 1 DFP line minimization
c imin 2 CG line minimization
cc ibrent : used for # of frozen atoms (always last ones)
c inmod vibrational analysis
c thermo thermodynamical properties
c
c Dario Estrin, 1992
c---------------------------------------------------
       subroutine geom(MEMO,FTOL,imin,ibrent,NORM,
     >  natom,Iz,r,Nuc,M,ncont,
     >  nshell,c,a,Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,Pm,delta,
     >  inmod,name1,ikk,thermo,TEMP,sigma,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
      implicit real*8 (a-h,o-z)
      character*20 strng
      character*17 name5
      character*12 name1
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write,MEMO
      logical OPEN,SVD,SHFT,GRAD,BSSE,integ,thereis
      logical field,sol,free
      integer nopt,iconst,igrid,imin,ibrent,thermo,igrid2
      INCLUDE 'param'
      parameter(nt3=nt*3,nx=nt3*2)
      dimension r(nt,3),nshelld(0:3),nshell(0:3),q(ntq)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
      dimension Pm(nt),vib(3*nt),Ei(3),Xi(20)
c
      dimension f(nt,3)
      dimension RMM(*)
*
      dimension xH(nt3,nt3),f1(nt,3),f2(nt,3),aux(nx)
      dimension xWW(nt3)
      dimension xU(nt3,3),xUU(nt3,3)
*
c auxiliars
c X scratch space in matrices
      dimension X(M,3*M),XX(Md,Md)
      dimension Em(ntq+nss),Rm(ntq+nss),pc(nss),alpha(nss)
c
c
      COMMON /TABLE/ STR(880,0:21)
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
      common /Sys/ SVD,iconst
c     common /HF/ nopt,OPEN,NMAX,NCO,ATRHO,VCINP,DIRECT,
c    >             IDAMP,EXTR,SHFT,SHI,GOLD,told,write,Nunp
c
      common /index/ index(ng)
      common /coef/ af(ngd)
      common /coef2/ B(ngd,2)
      common /Ngeom/ ngeo
      common /ENum/ GRAD
      common /propt/ idip,ipop
      common /cav/ a0,epsilon,field
      common /sol1/ Nsol,natsol,alpha,Em,Rm,pc,sol,free
c------------------------------------------------------------------
c
c Pointers
c
      MM=M*(M+1)/2 
      MM2=M**2
      M2=2*M
      MMd=Md*(Md+1)/2
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
*
      Nel=2*NCO+Nunp
      if(inmod.eq.1) goto 451
*
c
c First step : Call SCF subroutine
c
c minimization (RECOMENDED option)
        call dfp2(MEMO,FTOL,ITER,FRET,NORM,natom,Iz,r,Nuc,M,ncont,
     >   nshell,c,a,Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,xH,ibrent,
     >  nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
c
      write(*,*) iter,fret
      write(*,*) 'Number of geometries',ngeo
c
c PROPERTIES CALCULATIONS
       if (idip.eq.1) then
        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >       Nel,ux,uy,uz)
       endif
c
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
         write(*,*) 'MULLIKEN POPULATION ANALYSIS'
         write(*,770)

        do n=1,natom
         write(*,760) n,Iz(n),q(n)
        enddo
c
        write(*,*)
        endif
c--------------------------------------------------------------
c Outputs final  MO ---------------------
      if (OPEN) then
       NCOa=NCO
       NCOb=NCO+Nunp
       M18=M17+MMd
       M18b=M18+M*NCOa
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
c
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
      if (sol.and.free) then
       ntom=natom+Nsol*natsol
      else
       ntom=natom
      endif
c
      do 10 i=1,ntom
c      write(2,500) Iz(i),Pm(i),r(i,1),r(i,2),r(i,3)
       write(2,500) Iz(i),r(i,1),r(i,2),r(i,3)
  10  continue
c
c
 400  format(4(E14.7E2,2x))
 500  format(i3,2x,F9.6,2x,F9.6,2x,F9.6)
 760  format(I3,9x,I3,6x,F10.4)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
*
      if(inmod.eq.0)  return
*
*     Normal modes; finite-differences second derivatives calculations,
*                   masses in atomic units, final results in cm**-1
*     Intensities; finite-differences dipole first derivatives
*                  final results in km/mol (hopefully)
*
 451  continue
c
      if (sol.and.free) then
       ntom=natom+Nsol*natsol
      else
       ntom=natom
      endif
c
      numcalc = 0
      numcalctot = 6*ntom
      nat3 = ntom*3
      do i = 1,nt3
        do j = 1,3
          xU(i,j) = 0.D0
*          xUU(i,j) = 0.D0
        enddo
      enddo
      do i = 1,nt
       do j = 1,nt
         xH(i,j) = 0.D0 
       enddo 
      enddo
      del2 = delta/2.D0
      istart1 = 1
      istart2 = 1
****
      name5 = name1(1:ikk)//'.rest'                   !  RESTART
      inquire(file=name5, exist=thereis)
      if (thereis) then
         write(*,*)'THIS IS A RESTART RUN'   
         open(99,file=name5,form='unformatted')
         read(99) numcalc2
         do i = 1,numcalc2
            read(99) (xH(i,j),j=1,nat3)
            read(99) xU(i,1),xU(i,2),xU(i,3)
         enddo
      istart1 = numcalc2/3 +1
      istart2 = numcalc2 - (numcalc2/3)*3 + 1
      numcalc = numcalc2*2
      close(99)
      endif                                           !  END RESTART
****
      do 201 i = istart1,ntom
       sqmi = dsqrt(Pm(i))
       udenom = delta*sqmi
       do 51 k = istart2,3
        r(i,k) = r(i,k) + del2 
        GRAD=.false.
        numcalc = numcalc + 1
        write(*,663) numcalc,numcalctot
 663  format('Calculation no. ',i3,' of ',i3)
c
        if (field) then
        g0=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
        endif
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,FP,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,FP,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
c
        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux1,uy1,uz1)
        call int1G(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En,f1)
        call int3G(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >             Nucd,Md,ncontd,nshelld,cd,ad,RMM,E2,f1,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c continuum solvent case
      if (field) then
        g1=g0
c
        call dipg(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >               Nel,g1,ux1,uy1,uz1,f1)
      endif
c----------------------------
c classical solvent case ----
        if (sol) then
         call mmsolG(natom,Nsol,natsol,Iz,pc,r,Em,Rm,f1)
         call intsolG(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,f1)
        endif
c---------------------------
c
        call intSG(NORM,natom,r,Nuc,M,Md,ncont,nshell,c,a,RMM,f1)
c test
*
        r(i,k) = r(i,k) - delta 
        GRAD=.false. 
        numcalc = numcalc + 1
        write(*,663) numcalc,numcalctot
c
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,FP,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,FP,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
c
        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux2,uy2,uz2)
        call int1G(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En,f2)
        call int3G(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >             Nucd,Md,ncontd,nshelld,cd,ad,RMM,E2,f2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
*
c continuum solvent case ----
      if (field) then
        g1=g0
        call dipg(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >               Nel,g1,ux2,uy2,uz2,f2)
      endif
c
c----------------------------
c classical solvent case ----
        if (sol) then
         call mmsolG(natom,Nsol,natsol,Iz,pc,r,Em,Rm,f2)
         call intsolG(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,f2)
        endif
c---------------------------
        call intSG(NORM,natom,r,Nuc,M,Md,ncont,nshell,c,a,RMM,f2)
        ii = i*3 - (3-k)     
        xU(ii,1) = (ux1-ux2)/udenom
        xU(ii,2) = (uy1-uy2)/udenom
        xU(ii,3) = (uz1-uz2)/udenom
c     xH  Force-constant matrix
        do 101 j = 1,ntom
           hdenom = udenom*dsqrt(Pm(j))
           jj = j*3 - 2
           xH(ii,jj) = (f1(j,1) - f2(j,1))/hdenom
           xH(ii,jj+1) = (f1(j,2) - f2(j,2))/hdenom
           xH(ii,jj+2) = (f1(j,3) - f2(j,3))/hdenom
 101    continue
        r(i,k) = r(i,k) + del2 
*****
      open(99,file=name5,form='unformatted')      ! WRITING RESTART FILE
      numcalc2 = numcalc/2
      write(99) numcalc2
      do isf = 1,numcalc2
         write(99) (xH(isf,jsf),jsf=1,nat3)
         write(99) xU(isf,1),xU(isf,2),xU(isf,3)
      enddo
      close(99)                                   ! END WRITING
*****
        
  51   continue
      istart2 = 1
 201  continue
*
* Writing file .nmod
      write(45,*)'Force-constant matrix  (Hessian with respect to ',
     >           'mass-weighted coordinates)'
      do i = 1,nat3
         write(45,*)(xH(i,j),j=1,nat3)
      enddo
      write(45,*)'Dipole moment derivatives'
      do i = 1,nat3
         write(45,*)(xU(i,j),j=1,3)
      enddo 
*
*
*     upper-packed storage mode
      do 489 i = 1,nat3
       do j = i,nat3
         RMM(i + j*(j-1)/2) = (xH(i,j) + xH(j,i))*0.5D0
       enddo
 489  continue
      nnx = nat3*2
c ESSL option
#ifdef essl
      call dspev(21,RMM,xWW,XX,Md,nat3,aux,nx)
#endif
c LAPACK OPTION
#ifdef pack
      call dspev('V','U',nat3,RMM,xWW,XX,Md,aux,info)
#endif
*
*  CHECKS
*
*  Check if the system is linear
c Calculation of principal inertia moments
c
c
      Xcm=0.0D0
      Ycm=0.0D0
      Zcm=0.0D0
      Rcm=0.0D0
c
      do i=1,ntom
       Rcm=Rcm+Pm(i)
       Xcm=Xcm+Pm(i)*r(i,1)
       Ycm=Ycm+Pm(i)*r(i,2)
       Zcm=Zcm+Pm(i)*r(i,3)
      enddo
c
      Xcm=Xcm/Rcm
      Ycm=Ycm/Rcm
      Zcm=Zcm/Rcm
c
      do i=1,10
       Xi(i)=0.0D0
      enddo
c
      do i=1,ntom
       x1=r(i,1)-Xcm
       x2=x1**2
       y1=r(i,2)-Ycm
       y2=y1**2
       z1=r(i,3)-Zcm
       z2=z1**2
       r2=x2+y2+z2
c
       Xi(1)=Xi(1)+pm(i)*(r2-x2)
       Xi(4)=Xi(4)+pm(i)*(r2-y2)
       Xi(6)=Xi(6)+pm(i)*(r2-z2)
c
       Xi(2)=Xi(2)-pm(i)*x1*y1
       Xi(3)=Xi(3)-pm(i)*x1*z1
       Xi(5)=Xi(5)-pm(i)*y1*z1
c
      enddo
c
c ESSL OPTION ----------------------------------------
#ifdef essl
      call DSPEV(1,Xi,Ei,Wi,3,3,Xi(7),6)
#endif
c
c EISPACK OPTION -----------------------------------------------
#ifdef pack
c
       call dspev('N','L',3,Xi,Ei,Wi,3,Xi(7),info)
#endif
c-----------------------------------------------------------
c Ei(1), Ei(2) and Ei(3) three principal inertia moments
c if molecule is linear one should be zero and the other two equal
c in atoms the 3 of them should be 0
*
      if(abs(Ei(1)).lt.1.0D-06) then
        if (abs(Ei(2)).lt.1.0D-06) then
         ntrasrot =3
        else
         ntrasrot = 5
        endif
      else
        ntrasrot = 6
      endif
      nnull = nt3 - nat3 + ntrasrot +1
c
* Checks if the eigenvector is normalized and if not, normalizes it 
      itemp = ntrasrot + 1
*      do i = nat3,itemp,-1
      do i = nat3,1,-1
         vecsum = 0.D0
         do j = 1,nat3
            vecsum = vecsum + XX(j,i)**2
         enddo
         check = dabs(vecsum - 1.D0)
         if (check.gt.1.D-6) then
            write(*,*)'check.gt.1.D-6',vecsum
            xnorm = dsqrt(vecsum)
            do j = 1,nat3
               XX(j,i) = XX(j,i)/xnorm
cTEST
            write(17,*) XX(j,i),j
            enddo
         endif
      enddo
*  END CHECKS
*
      zpecm = 0.D0
      do i = 1,nt3
         do j = 1,3
            xUU(i,j) = 0.D0
         enddo
      enddo
*      do i = nat3,itemp,-1
      do i = nat3,1,-1
         do j = 1,nat3
            autovet = XX(j,i)
            xUU(i,1) = xUU(i,1) + autovet*xU(j,1)
            xUU(i,2) = xUU(i,2) + autovet*xU(j,2)
            xUU(i,3) = xUU(i,3) + autovet*xU(j,3)
         enddo
      enddo
*      do i = nat3,itemp,-1 
      do i = nat3,1,-1			          ! restoring
         do j = 1,nat3				  ! cartesian
            indm = (j-1)/3 + 1			  ! coordinates
            XX(j,i) = XX(j,i)/dsqrt(Pm(indm))	  ! (no more 
         enddo					  !  mass-weighted) 
      enddo					  ! 
*
*  OUTPUTS :   unit 44 = .vib      unit 46 = .freq      
*
*  Conversion Factors :
*                      const1 converts frequency to cm**-1 
*                      const2 converts intensity to km/mol
*                      const3 converts z.p. energy to kcal/mol
      const1 = 5141.1182009496200901D0
      const2 = 150.89915D0
      const3 = 1.430D-3
c
      const4= 1.43*2.0D0*4.186/(8.31451*TEMP)
      const5=TEMP*0.011545543
      const6=-1.164856775
      const7=2.072364943
c
      write(44,662)ntom
 662  format('Number of Atoms and Molecular Geometry   (a.u.)',/,i2)
      do i=1,ntom
       write(44,500) Iz(i),r(i,1),r(i,2),r(i,3)
      enddo
      write(44,661) nat3 - ntrasrot
 661  format(/,'Total Number of Normal Modes',/,i3)
*
      write(46,*)'FREQUENCES (cm**-1) and INTENSITIES (kcal/mol)'
      do i = nat3,1,-1
           xUU(i,1) = xUU(i,1)*xUU(i,1)
           xUU(i,2) = xUU(i,2)*xUU(i,2)
           xUU(i,3) = xUU(i,3)*xUU(i,3)
           if (xWW(i).lt.0.D0) then
              write(46,*)'Purely imaginary eigenvalue'
              xWW(i) = dabs(xWW(i))
           endif
           freq = dsqrt(xWW(i))
           freqcm = freq * const1
           xintensity = const2 * (xUU(i,1)+xUU(i,2)+xUU(i,3))
           write(46,*) freqcm, xintensity 
           if(i.ge.itemp) then
              zpecm = zpecm + freqcm
              vib(i)=freqcm*const4
           endif
c
           ni = nat3 - i + 1
           if(i.eq.(itemp)) write(46,669)
           write(44,664) ni,freqcm, xintensity 
           write(44,667) (XX(iku,i),iku=1,nat3)
      enddo
      zpecm = zpecm * const3 
      write(46,*)'Zero-point energy (kcal/mol)  ',zpecm
      write(44,668) zpecm
*
c THERMO OPTION -------------------------------------------------
c ALL THIS INFORMATION CORRESPONDS TO THE FOLLOWING MODEL:
c 1) RIGID ROTOR
c 2) HARMONIC UNCOUPLED OSCILLATORS
c 3) NO INTERNAL ROTATION
c 4) NO SYMMETRY FACTORS
c
c ENTROPY IN  IN CAL/K/MOL - ENTHAPLY CAL/MOL - HEAT CAPACITY CAL/K/MOL
c
      if (thermo.gt.0) then
c
c Translational entropy and enthalpy
c enthalpy : energy + RT, Pressure considered 1 atm
      Prs=1.0D0
      St=const6+1.50D0*log(Rcm)+2.50D0*log(TEMP)-log(Prs)
      St=St*8.31451/4.186
c   
      Et=2.50D0*8.31451/4.186*TEMP
      Ct=1.50D0*8.31451/4.186
c
c Rotational entropy and energy
c
       if (abs(Ei(1)).lt.1.0D-06) then
        if (abs(Ei(2)).gt.1.0D-04) then
c linear molecules
         Er=8.31451/4.186*TEMP
c        Sr=const7+0.50D0*log(const5**3*Ei(2)*Ei(3))
         Sr=log(const5*Ei(2)/sigma)+1.0D0
         Cr=8.31451/4.186
        else
c atoms
         Er=0.0D0
         Sr=0.0D0
         Cr=0.0D0
        endif
       else
c non-linear molecules
       Er=1.50D0*8.31451/4.186*TEMP
       Sr=const7+0.50D0*log(const5**3*Ei(1)*Ei(2)*Ei(3))-log(sigma)
       Cr=1.50D0*8.31451/4.186
       endif
c
      Sr=Sr*8.31451/4.186
c Vibrational entropy
c
      Sv=0.0D0
      Ev=0.0D0
      Cv=0.0D0
       do i=itemp,3*natom
        if (vib(i).gt.100.D0) then
       rterm=0.0D0
       rterm2=0.0D0
      else
       rt1=exp(vib(i))
       rterm=vib(i)/(rt1-1.0D0)
       rterm2=rt1/(rt1-1.0D0)**2
      endif
      Sv=Sv+rterm-log(1.0D0-exp(-vib(i)))
c first expression includes zero point energy
c     Ev=Ev+0.50D0*vib(i)+rterm
      Ev=Ev+rterm
      Cv=Cv+vib(i)**2*rterm2
      enddo
      Sv=Sv*8.31451/4.186
      Ev=Ev*8.31451/4.186*TEMP
      Cv=Cv*8.31451/4.186
c
      write(44,700) TEMP,Prs
      write(44,*)
      write(44,701)
      write(44,702)
c     
      write(44,703) Et,Ct,St
      write(44,704) Er,Cr,Sr
      write(44,705) Ev,Cv,Sv
c
      endif
c
 664  format(/,' Normal Mode N.      Frequency  (cm**-1)',
     >       '    IR Intensity  (Km/mol)  ',/,4x,i3,15x,f10.0,12x,f10.0)
 667  format(3(f10.5,3x))
 668  format(/,'Zero-point Energy (kcal/mol)  ',f10.6)
 669  format(/,'Eigenvalues of rotations and translations  ',
     >       '(they should be zero;',/,'if not, they can be purely',
     >       ' imaginary eigenvalues)')
*
 700  format(/,'THERMODYNAMIC PROPERTIES AT TEMPERATURE = ',f4.0,'K',3x,
     > 'PRESSURE= ',f5.2,' ATM')
 701  format(/,'CONTRIBUTION',5x,'ENTHALPY',5x,'HEAT CAPACITY',
     >  5x,'ENTROPY')
 702  format(17x,'CAL/MOL',6x,'CAL/MOL/K',9x,'CAL/MOL/K')
 703  format('TRANSLATIONAL',3x,f10.4,4x,f8.4,10x,f8.4)
 704  format('ROTATIONAL',6x,f10.4,4x,f8.4,10x,f8.4)
 705  format('VIBRATIONAL',5x,f10.4,4x,f8.4,10x,f8.4)
c
      close(44)
      close(45)
      close(46)
      strng='rm '//name5
      call system(strng)
*
      write(*,*) 'Time for normal modes ', (itimefinal-itimeinit)/100
      return
      end
