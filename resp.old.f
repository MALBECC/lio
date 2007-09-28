c This subroutine is used for the evaluation of the electrical
c properties.
c The dipole moment, polarizability
c
c calls SCF subroutine for calculating energy at
c different intensities of the applied field
c Fi, -Fi, 2Fi, -2Fi
c with i =x,y,z
c this correspond to 12 calculations
c (A standard SCF calculations, with reaction field,
c  performed by calling dip2, but with the applied
c field instead of the reaction field)
c The electrical parameters are evaluated using
c the energy expansion, or the dipole expansion
c See J. Phys. Chem. 98, 2545 1994 and J. Phys. Chem. 97,
c 1158 1993 for details
c
c In the case of using the dipole expansion, after each
c energy calculations, dip should de called
c
c F intensity of applied field
c u dipole moment vector
c alpha diagonal terms of polarizability tensor
c beta diagonal terms of first order hiporpelarizability tensor
c gamma diagonal terms of second order hiperpolarizability tensor
c
c all specifications given in the namelist SCF are necessary
c
c Dario Estrin, January 1995
c---------------------------------------------------
       subroutine resp(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,
     >  nshell,c,a,Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,F,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
      implicit real*8 (a-h,o-z)
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write
      logical OPEN,SVD,SHFT,GRAD,BSSE,integ,field,sol,free
      logical exter,MEMO
      integer nopt,iconst,igrid,igrid2
      INCLUDE 'param'
      dimension r(nt,3),nshelld(0:3),nshell(0:3),q(ntq)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
      dimension RMM(*)
*
*
c auxiliars
c X scratch space in matrices
      dimension X(M,3*M),XX(Md,Md)
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
      common /cav/ a0,epsilon,field,exter,Fx,Fy,Fz
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
*
c
      Fx=0.0D0
      Fy=0.0D0
      Fz=0.0D0
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E0,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E0,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux,uy,uz)
c
c Calculation of x components ------------------------
c
      Fx=F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E11,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E11,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux11,uy11,uz11)
c
        Fx=-F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em11,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em11,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxm11,uym11,uzm11)

c
        Fx=2.0D0*F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux2,uy2,uz2)
c
        Fx=-2.0D0*F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxm2,uym2,uzm2)
c
c Calculations of properties
c
c USING ENERGY EXPANSION
       uex=(-0.6666666666D0*(E11-Em11)+0.083333333333D0*(E2-Em2))/F
       alphaex=(2.5D0*E0-1.3333333333D0*(E11+Em11)+0.083333333333D0*
     >         (E2+Em2))/F**2
       btexxx= (0.5D0*(E11-Em11)-0.25D0*(E2-Em2))/F**3
       gmex=(-E0+0.666666666666D0*(E11+Em11)-0.1666666666666D0*
     >         (E2+Em2))/F**4     
c
c USING DIPOLE EXPANSION
       alphaux=(0.66666666666D0*(ux11-uxm11)-0.083333333333D0*
     >         (ux2-uxm2))/(F*2.54D0)
       alpuyx=(0.66666666666D0*(uy11-uym11)-0.083333333333D0*
     >         (uy2-uym2))/(F*2.54D0)
       alpuzx=(0.66666666666D0*(uz11-uzm11)-0.083333333333D0*
     >         (uz2-uzm2))/(F*2.54D0)
       btuxxx=(-1.25D0*ux+0.6666666666666D0*(ux11+uxm11)-  
     >         0.0416666666D0*(ux2+uxm2))/(F**2*2.54D0)
       btuyxx=(0.333333333333D0*(uy2+uym2-uy11-uym11))/(F**2*2.54D0)
       btuzxx=(0.333333333333D0*(uz2+uzm2-uz11-uzm11))/(F**2*2.54D0)
       gmux=(-0.1666666666666D0*(ux11-uxm11)+0.08333333333333D0*
     >         (ux2-uxm2))/(F**3*2.54D0)
c
c Calculation of y components ------------------------
c
      Fx=0.0D0
      Fy=F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E21,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E21,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux21,uy21,uz21)
c
        Fy=-F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em21,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em21,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)

        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxm21,uym21,uzm21)
c
        Fy=2.0D0*F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux2,uy2,uz2)
c
        Fy=-2.0D0*F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxm2,uym2,uzm2)
c
c Calculations of properties
c
c USING ENERGY EXPANSION
       uey=(-0.6666666666D0*(E21-Em21)+0.083333333333D0*(E2-Em2))/F
       alphaey=(2.5D0*E0-1.3333333333D0*(E21+Em21)+0.083333333333D0*
     >         (E2+Em2))/F**2
       bteyyy= (0.5D0*(E21-Em21)-0.25D0*(E2-Em2))/F**3
       gmey=(-E0+0.666666666666D0*(E21+Em21)-0.1666666666666D0*
     >         (E2+Em2))/F**4
c
c USING DIPOLE EXPANSION
       alphauy=(0.66666666666D0*(uy21-uym21)-0.083333333333D0*
     >         (uy2-uym2))/(F*2.54D0)
       alpuxy=(0.66666666666D0*(ux21-uxm21)-0.083333333333D0*
     >         (ux2-uxm2))/(F*2.54D0)
       alpuzy=(0.66666666666D0*(uz21-uzm21)-0.083333333333D0*
     >         (uz2-uzm2))/(F*2.54D0)
       btuyyy=(-1.25D0*uy+0.6666666666666D0*(uy21+uym21)-
     >         0.04166666666666D0*(uy2+uym2))/(F**2*2.54D0)
       btuxyy=(0.3333333333333D0*(ux2+uxm2-ux21-uxm21))/(F**2*2.54D0)
       btuzyy=(0.3333333333333D0*(uz2+uzm2-uz21-uzm21))/(F**2*2.54D0)
       gmuy=(-0.1666666666666D0*(uy21-uym21)+0.08333333333333D0*
     >         (uy2-uym2))/(F**3*2.54D0)
c
c------------------------------------------------------------
c Calculation of z components
c
      Fy=0.0D0
      Fz=F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E31,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E31,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux31,uy31,uz31)
c
        Fz=-F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em31,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em31,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxm31,uym31,uzm31)

c
        Fz=2.0D0*F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,ux2,uy2,uz2)
c
        Fz=-2.0D0*F
        if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Em2,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxm2,uym2,uzm2)
c
c Calculations of properties
c
c USING ENERGY EXPANSION
       uez=(-0.6666666666D0*(E31-Em31)+0.083333333333D0*(E2-Em2))/F
       alphaez=(2.5D0*E0-1.3333333333D0*(E31+Em31)+0.083333333333D0*
     >         (E2+Em2))/F**2
       btezzz= (0.5D0*(E31-Em31)-0.25D0*(E2-Em2))/F**3
       gmez=(-E0+0.666666666666D0*(E31+Em31)-0.1666666666666D0*
     >         (E2+Em2))/F**4
c
c USING DIPOLE EXPANSION
       alphauz=(0.66666666666D0*(uz31-uzm31)-0.083333333333D0*
     >         (uz2-uzm2))/(F*2.54D0)
       alpuxz=(0.66666666666D0*(ux31-uxm31)-0.083333333333D0*
     >         (ux2-uxm2))/(F*2.54D0)
       alpuyz=(0.66666666666D0*(uy31-uym31)-0.083333333333D0*
     >         (uy2-uym2))/(F*2.54D0)
       btuzzz=(-1.25D0*uz+0.6666666666666D0*(uz31+uzm31)-
     >         0.04166666666666666D0*(uz2+uzm2))/(F**2*2.54D0)
       btuxzz=(0.33333333333333D0*(ux2+uxm2-ux31-uxm31))/(F**2*2.54D0)
       btuyzz=(0.33333333333333D0*(uy2+uym2-uy31-uym31))/(F**2*2.54D0)
       gmuz=(-0.1666666666666D0*(uz31-uzm31)+0.08333333333333D0*
     >         (uz2-uzm2))/(F**3*2.54D0)
c
c Calculation of non diagonal terms ------------------------------
c
c si se quiere solo la polarizabilidad saltearse este pedazo!!!!!
c
c
       goto 876
       Fx=F
       Fy=F
       Fz=0
                 if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ea)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ea)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxa,uya,uza)
c
       Fx=F
       Fy=-F
       Fz=0
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eb)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eb)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxb,uyb,uzb)
c
       Fx=-F
       Fy=F
       Fz=0
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ec)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ec)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxc,uyc,uzc)
c
       Fx=-F
       Fy=-F       
       Fz=0
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ed)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ed)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxd,uyd,uzd)
c
       Fx=2D0*F
       Fy=2D0*F
       Fz=0
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ee)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ee)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxe,uye,uze)
c
       Fx=2D0*F
       Fy=-2D0*F
       Fz=0
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ef)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ef)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxf,uyf,uzf)
c
       Fx=-2D0*F
       Fy=2D0*F
       Fz=0
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eg)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eg)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxg,uyg,uzg)
c
       Fx=-2D0*F
       Fy=-2D0*F
       Fz=0
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eh)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eh)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxh,uyh,uzh)
c
c Calculation of properties
c
c USING ENERGY EXPANSION
         alpexy=(0.02083333333D0*(Ee-Ef-Eg+Eh)-0.33333333333D0*
     >           (Ea-Eb-Ec+Ed))/F**2
         alpeyx=(0.02083333333D0*(Ee-Eg-Ef+Eh)-0.33333333333D0*
     >           (Ea-Ec-Eb+Ed))/F**2
         btexyy=(0.5D0*(Ed-Ea+Ec-Eb)+E11-Em11)/F**3
         bteyxx=(0.5D0*(Ed-Ea+Eb-Ec)+E21-Em21)/F**3
         gmexxyy=(-4.0D0*E0-Ea-Ed-Eb-Ec+2D0*(E11+Em11+E21+Em21))/F**4
         gmeyyxx=(-4.0D0*E0-Ea-Ed-Ec-Eb+2D0*(E21+Em21+E11+Em11))/F**4
c
c USING DIPOLE EXPANSION
         gmuxxyy=(0.5D0*(uxa-uxc+uxb-uxd)-ux11+uxm11)/(F**3*2.54D0)
         gmuyyxx=(0.5D0*(uya-uyb+uyc-uyd)-uy21+uym21)/(F**3*2.54D0)
c--------------------------------------------------------------------------
       Fx=0
       Fy=F
       Fz=F
                 if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ea)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ea)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxa,uya,uza)
c
       Fx=0
       Fy=F
       Fz=-F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eb)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eb)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxb,uyb,uzb)
c
       Fx=0
       Fy=-F
       Fz=F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ec)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ec)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxc,uyc,uzc)
c
       Fx=0
       Fy=-F
       Fz=-F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ed)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ed)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxd,uyd,uzd)
c
       Fx=0
       Fy=2D0*F
       Fz=2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ee)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ee)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxe,uye,uze)
c
       Fx=0
       Fy=2D0*F
       Fz=-2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ef)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ef)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxf,uyf,uzf)
c
       Fx=0
       Fy=-2D0*F
       Fz=2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eg)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eg)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxg,uyg,uzg)
c
       Fx=0
       Fy=-2D0*F
       Fz=-2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eh)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eh)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxh,uyh,uzh)
c
c USING ENERGY EXPANSION
         alpeyz=(0.02083333333D0*(Ee-Ef-Eg+Eh)-0.33333333333D0*
     >           (Ea-Eb-Ec+Ed))/F**2
         alpeyz=(0.02083333333D0*(Ee-Eg-Ef+Eh)-0.33333333333D0*
     >           (Ea-Ec-Eb+Ed))/F**2
         bteyzz=(0.5D0*(Ed-Ea+Ec-Eb)+E21-Em21)/F**3
         btezyy=(0.5D0*(Ed-Ea+Eb-Ec)+E31-Em31)/F**3
         gmeyyzz=(-4.0D0*E0-Ea-Ed-Eb-Ec+2D0*(E21+Em21+E31+Em31))/F**4
         gmezzyy=(-4.0D0*E0-Ea-Ed-Ec-Eb+2D0*(E31+Em31+E21+Em21))/F**4
c
c USING DIPOLE EXPANSION
         gmuyyzz=(0.5D0*(uya-uyc+uyb-uyd)-uy21+uym21)/(F**3*2.54D0)
         gmuzzyy=(0.5D0*(uza-uzb+uzc-uzd)-uz31+uzm31)/(F**3*2.54D0)
c ----------------------------------------------------------------
       Fx=F
       Fy=0
       Fz=F
                 if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ea)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ea)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxa,uya,uza)
c
       Fx=-F
       Fy=0
       Fz=F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eb)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eb)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxb,uyb,uzb)
c
       Fx=F
       Fy=0
       Fz=-F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ec)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ec)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxc,uyc,uzc)
c
       Fx=-F
       Fy=0
       Fz=-F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ed)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ed)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxd,uyd,uzd)
c
       Fx=2D0*F
       Fy=0
       Fz=2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ee)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ee)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxe,uye,uze)
c
       Fx=2D0*F
       Fy=0
       Fz=-2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ef)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ef)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxf,uyf,uzf)
c
       Fx=-2D0*F
       Fy=0
       Fz=2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eg)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eg)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxg,uyg,uzg)
c
       Fx=-2D0*F
       Fy=0
       Fz=-2D0*F
                if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eh)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Eh)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxh,uyh,uzh)
c
c Calculation of properties
c 
c
c USING ENERGY EXPANSION
         alpezx=(0.02083333333D0*(Ee-Ef-Eg+Eh)-0.33333333333D0*
     >           (Ea-Eb-Ec+Ed))/F**2
         alpexz=(0.02083333333D0*(Ee-Eg-Ef+Eh)-0.33333333333D0*
     >           (Ea-Ec-Eb+Ed))/F**2
         btezxx=(0.5D0*(Ed-Ea+Ec-Eb)+E31-Em31)/F**3
         btexzz=(0.5D0*(Ed-Ea+Eb-Ec)+E11-Em11)/F**3
         gmezzxx=(-4.0D0*E0-Ea-Ed-Eb-Ec+2D0*(E31+Em31+E11+Em11))/F**4
         gmexxzz=(-4.0D0*E0-Ea-Ed-Ec-Eb+2D0*(E11+Em11+E31+Em31))/F**4
c
c USING DIPOLE EXPANSION
         gmuzzxx=(0.5D0*(uza-uzc+uzb-uzd)-uz31+uzm31)/(F**3*2.54D0)
         gmuxxzz=(0.5D0*(uxa-uxb+uxc-uxd)-ux11+uxm11)/(F**3*2.54D0)
c ---------------------------------------------------------------
c Calculation of others terms
c
       Fx=F
       Fy=F
       Fz=F
                 if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew1)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew1)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxw1,uyw1,uzw1)
c
       Fx=-F
       Fy=-F
       Fz=-F
                 if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew3)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew3)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxw3,uyw3,uzw3)
c
       Fx=2.0D0*F
       Fy=2.0D0*F
       Fz=2.0D0*F
                 if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew2)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew2)
        endif
         call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >           Nel,uxw2,uyw2,uzw2)
c
       Fx=-2.0D0*F
       Fy=-2.0D0*F
       Fz=-2.0D0*F
                 if (OPEN) then
        call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew4)
        else
        call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >           Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,Ew4)
       endif
        call dip(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,RMM, 
     >           Nel,uxw4,uyw4,uzw4)
c
c Calculation of properties
        suma=(0.5D0*(2.0D0*Ew3-2.0D0*Ew1-Ew4+Ew2))/F**3
c -------------------------------------------------------------
  876    continue
c Calculation of observables
c
         u=(uex**2+uey**2+uez**2)**(0.50D0)
         alphae=0.333333333333D0*(alphaex+alphaey+alphaez)
         alphau=0.333333333333D0*(alphaux+alphauy+alphauz)
         betaex=btexxx+btexyy+btexzz
         betaey=bteyyy+bteyxx+bteyzz
         betaez=btezzz+btezxx+btezyy
         betae=(betaex**2+betaey**2+betaez**2)**0.5D0 
         betaux=btuxxx+btuxyy+btuxzz
         betauy=btuyyy+btuyxx+btuyzz
         betauz=btuzzz+btuzxx+btuzyy
         betau=(betaux**2+betauy**2+betauz**2)**0.5D0 
         gammae=0.2D0*(gmex+gmez+gmey+2D0*(gmezzxx+
     >          gmezzyy+gmeyyxx))
         gammau=0.2D0*(gmux+gmuz+gmuy+2D0*(gmuzzxx+
     >          gmuzzyy+gmuyyxx))
         sumbet=suma-(betaex+betaey+betaez)
c ---------------------------------------------------------------
       write(*,*) 'ELECTRICAL RESPONSE CALCULATION'
       write(*,*) 'DIAGONAL TERMS XYZ'
       write(*,*) 
       write (*,*) 'RESULTS OBTAINED USING ENERGY EXPANSION'
       write(*,700) uex,uey,uez
       write(*,800) alphaex,alphaey,alphaez
       write(*,900) btexxx,bteyyy,btezzz
       write(*,950) gmex,gmey,gmez
       write(*,*) 
       write(*,*) 'RESULTS OBTAINED USING DIPOLE EXPANSION'
       write(*,800) alphaux,alphauy,alphauz
       write(*,900) btuxxx,btuyyy,btuzzz
       write(*,950) gmux,gmuy,gmuz
       write(*,*)
       write(*,*)  'NON DIAGONAL TERMS'
       write(*,*)
       write(*,*) 'RESULTS OBTAINED USING ENERGY EXPANSION' 
       write(*,*) 'Alpha xy   Alpha xz'
       write(*,910) alpexy,alpexz
       write(*,*) 'Alpha yx   Alpha yz'
       write(*,910) alpeyx,alpeyz
       write(*,*) 'Alpha zx   Alpha zy'
       write(*,910) alpezx,alpezy
       write(*,*) 'Beta xyy   Betaxzz'
       write(*,910) btexyy,btexzz
       write(*,*) 'Beta yxx   Betayzz'
       write(*,910) bteyxx,bteyzz   
       write(*,*) 'Beta zxx   Betazyy'
       write(*,910) btezxx,btezyy 
       write(*,*) 'Gammaxxyy       Gammaxxzz'
       write(*,*) gmexxyy,gmexxzz
       write(*,*) 'Gammayyxx       Gammayyzz'
       write(*,*) gmeyyxx,gmeyyzz
       write(*,*) 'Gammazzxx       Gammazzyy'
       write(*,*) gmezzxx,gmezzyy
       write(*,*)
       write(*,*) 'RESULTS OBTAINED USING DIPOLE EXPANSION' 
       write(*,*) 'Alpha xy   Alpha xz'
       write(*,910) alpuxy,alpuxz
       write(*,*) 'Alpha yx   Alpha yz'
       write(*,910) alpuyx,alpuyz
       write(*,*) 'Alpha zx   Alpha zy'
       write(*,910) alpuzx,alpuzy
       write(*,*) 'Beta xyy   Betaxzz'
       write(*,910) btuxyy,btuxzz
       write(*,*) 'Beta yxx   Betayzz' 
       write(*,910) btuyxx,btuyzz   
       write(*,*) 'Beta zxx   Betazyy' 
       write(*,910) btuzxx,btuzyy
       write(*,*) 'Gammaxxyy       Gammaxxzz'
       write(*,*) gmuxxyy,gmuxxzz
       write(*,*) 'Gammayyxx       Gammayyzz'
       write(*,*) gmuyyxx,gmuyyzz
       write(*,*) 'Gammazzxx       Gammazzyy'
       write(*,*) gmuzzxx,gmuzzyy
       write(*,*) 'OBSERVABLES'
       write(*,*) 'RESULTS OBTAINED USING ENERGY EXPANSION'
       write(*,940) alphae
       write(*,*) 'Beta x   Beta y   Beta z'
       write(*,820) betaex,betaey,betaez
       write(*,960) Sumbet
       write(*,930) betae
       write(*,920) gammae
       write(*,*) 'RESULTS OBTAINED USING DIPOLE EXPANSION'
       write(*,940) alphau
       write(*,*) 'Beta x   Beta y   Beta z'
       write(*,820) betaux,betauy,betauz
       write(*,930) betau
       write(*,920) gammau
c---------
 700   format('DIPOLE MOMENT =',3(1x,F7.4))
 800   format('POLARIZABILITY =',3(1x,F7.4))
 820   format(3(1x,F15.6))
 900   format('HIPERPOLARIZABILITY 1st =',3(1x,F15.4))
 910   format(2(1x,F15.6))
 920   format('GAMMA =',(1x,F20.5))
 930   format('BETA =',(1x,F20.5)) 
 940   format('ALPHA =',(1x,F7.4))
 950   format('HIPERPOLARIZABILITY 2nd =',3(1x,F9.1))  
 960   format('SUMA NON DIAGONAL BETA =',(1x,F9.4))
      return
      end




