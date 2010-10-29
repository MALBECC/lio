c     ===========================================
c BROYDEN-FLETCHER-GOLDFARB-SHANNO SUBROUTINE,
c using search, for line minimizations, no bracketing,
c searchs for 0 of gradient
c
c----------------------------------------------------------

      SUBROUTINE dfp2(MEMO,GTOL,ITER,FRET,NORM,natom,Iz,r,Nuc,M
     > ,ncont,nshell,c,a,Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,HESSIN,
     >   nfrozen,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param'
c  STPMX 0.05 original
      PARAMETER (NMAX=3*nt,STPMX=0.01D0,ITMAX=200,EPS=3.E-08)
      PARAMETER (TOLX=4.D0*EPS)
      DIMENSION P(NMAX),HESSIN(NMAX,NMAX),XI(NMAX),G(NMAX),DG(NMAX),
     > HDG(NMAX),f(NMAX/3,3),PNEW(NMAX)
c
c
c
      logical SHFT,NORM,ATRHO,VCINP,DIRECT,EXTR,dens,SVD
      logical OPEN,GRAD,integ,write,field,sol,free,MEMO
      integer nopt,iconst,igrid,iforce,nfrozen
      dimension r(nt,3),nshelld(0:3),nshell(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
      dimension RMM(*),X(M,3*M),XX(Md,Md)
      dimension Em(ntq+nss),Rm(ntq+nss),pc(nss),alpha(nss)
      character*2 atname

c
      COMMON /TABLE/ STR(880,0:21)
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
c     common /HF/ nopt,OPEN,NMX,NCO,ATRHO,VCINP,DIRECT,
c    >             IDAMP,EXTR,SHFT,SHI,GOLD,told,write,Nunp
      common /Sys/ SVD,iconst
      common /coef2/ B(ngd,2)
      common /coef/ af(ngd)
      common /force/ iforce
      common /ENum/ GRAD
      common /cav/ a0,epsilon,field
      common /sol1/ Nsol,natsol,alpha,Em,Rm,pc,sol,free
c
      iforce=0
c
      if (field) then
       g0=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
      endif
c
      zero=0.0D0
      Nel=2*NCO+Nunp
      write(*,*) 'NEW DAVIDSON-FLETCHER-POWELL OPTIMIZATION'
      write(*,*)
c
      ntom=natom-nfrozen
c
      if (sol.and.free)then
       ntom=ntom+Nsol*natsol
      endif
c
      N=3*ntom
c
      k=0
      do 20 i=1,ntom
        k=k+1
        P(k)=r(i,1)
        k=k+1
        P(k)=r(i,2)
        k=k+1
        P(k)=r(i,3)
 20   continue
c
c
c
c first call , needs also E , not only gradients -----------
c

#ifdef G2G
			write(*,*) 'primera carga de posiciones (con igrid2)'
			call g2g_reload_atom_positions(igrid2)
#endif

      GRAD=.false.
      if (OPEN) then
      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,FP,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
       else
      call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,FP,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      endif
c
c now gradients 
c
c
#ifdef G2G
			write(*,*) 'cambio de grilla para fuerza+energia (igrid)'
      call g2g_new_grid(igrid)
#endif

      call int1G(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En,f)
c
c
      call int3G(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,Exc,f,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write,.true.)
c
      call intSG(NORM,natom,r,Nuc,M,Md,ncont,nshell,c,a,RMM,f)
c
c reaction field case -------
      if (field) then
        g1=g0
        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,Nel,
     >       ux,uy,uz)
c
        call dipg(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
     >               Nel,g1,ux,uy,uz,f)
      endif
c----------------------------
c classical solvent case ----
        if (sol) then
         call mmsolG(natom,Nsol,natsol,Iz,pc,r,Em,Rm,f)
         call intsolG(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     >            ncont,nshell,c,a,pc,RMM,f)
        endif
c---------------------------
      FP=FP+Exc
      write(*,480) FP
c ----------------------------------------------------------
      write(*,*)
c
      write(*,*) 'ENERGY GRADIENTS, IN A.U'
      k=0
      ss=0.D0
      ss1=0.0D0
      do 29 i=1,ntom
        ss=ss+f(i,1)**2+f(i,2)**2+f(i,3)**2
        ss1=ss1+r(i,1)**2+r(i,2)**2+r(i,3)**2
        k=k+1
        G(k)=f(i,1)
        XI(k)=-f(i,1)
        k=k+1
        G(k)=f(i,2)
        XI(k)=-f(i,2)
        k=k+1 
        G(k)=f(i,3)
        XI(k)=-f(i,3)
        write(*,500) i,f(i,1),f(i,2),f(i,3)
 29    continue
c
      if (ibrent.eq.2) then
      do 31 i=ntom+1,natom
        write(*,500) i,zero,zero,zero
 31   continue
      endif
c
      stpmax=STPMX*dmax1(dsqrt(ss1),dfloat(N))
      write(*,450) sqrt(ss)
      write(*,*)
      if (ss.lt.1.D-08) then
       write(*,*) 'INITIAL GEOMETRY ALREADY OPTIMIZED'
       return
      endif
c
      DO 12 I=1,N
        DO 11 J=1,N
          HESSIN(I,J)=0.D0
11      CONTINUE
        HESSIN(I,I)=1.D0
12    CONTINUE
c
      DO 27 ITS=1,ITMAX
        ITER=ITS
c
        CALL LSEARCH(MEMO,N,P,FP,G,XI,PNEW,
     >               FRET,STPMAX,CHECK,NORM,natom,Iz,
     >  r,Nuc,M,ncont,nshell,c,a,Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,
     >  XX,nfrozen,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
        FP=FRET
       open(unit=69,file='opt.xyz')
       A0=0.5291771
        write(69,*) int(N/3)
        write(69,*)
       atname='??' 
        do 121 i=1,N,3
       ii=(i/3)+1
      if(Iz(ii).eq.1) atname='H'
      if(Iz(ii).eq.6) atname='C'
      if(Iz(ii).eq.7) atname='N'
      if(Iz(ii).eq.8) atname='O'
      if(Iz(ii).eq.3) atname='Li'
       write(69,*)  atname,PNEW(i)*A0,PNEW(i+1)*A0,PNEW(i+2)*A0
 121     write(*,501) PNEW(i),PNEW(i+1),PNEW(i+2)
c
        if (ibrent.eq.2) then
        do 122 i=natom-nfrozen+1,natom
         write(*,501) r(i,1),r(i,2),r(i,3)
 122    continue
        endif
c
        DO 13 I=1,N
         XI(I)=PNEW(I)-P(I)
         P(I)=PNEW(I)
 13     CONTINUE
c
        k=0
        do 30 i=1,ntom
        k=k+1
        r(i,1)=P(k)
        k=k+1
        r(i,2)=P(k)
        k=k+1
        r(i,3)=P(k)
30      continue
c
        test=0.0D0
c
      DO 14 i=1,N
        temp=abs(XI(i))/dmax1(abs(P(i)),1.D0)
        if (temp.gt.test) test=temp
 14   CONTINUE
c
      IF (test.lt.TOLX) return
c
c
        DO 15 I=1,N
          DG(I)=G(I)
15      CONTINUE

#ifdef G2G
			write(*,*) 'actualizacion de posiciones por movimiento'
			call g2g_reload_atom_positions(igrid2)
#endif
c
      GRAD=.false.
      if (OPEN) then
      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
       else
      call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)

      endif
c
c now gradients
c
c
#ifdef G2G
			write(*,*) 'cambio de grilla para fuerza+energia (igrid)'
      call g2g_new_grid(igrid)
#endif


      call int1G(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En,f)
c
      call int3G(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,Exc,f,
     > nopt,OPEN,NMX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write,.true.)
c
      call intSG(NORM,natom,r,Nuc,M,Md,ncont,nshell,c,a,RMM,f)
c reaction field case ------
      if (field) then
        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,Nel,
     >       ux,uy,uz)
c
        g1=g0
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
c
      E=E+Exc
      write(*,480) E
c ----------------------------------------------------------
c
c
        k=0
        ss=0.0D0
        do 140 i=1,ntom
        ss=ss+f(i,1)**2+f(i,2)**2+f(i,3)**2
        k=k+1
        G(k)=f(i,1)
        k=k+1
        G(k)=f(i,2)
        k=k+1 
        G(k)=f(i,3)
c test
c        write(*,500) i,f(i,1),f(i,2),f(i,3)       
140   continue
c
      test=0.0D0
c

      den=dmax1(fret,1.D0)
      write(*,450) sqrt(ss)
c
        DO 16 I=1,N
         temp=abs(g(i))*dmax1(abs(P(i)),1.D0)/den
         if (temp.gt.test) test=temp
 16     CONTINUE
c
       if (test.lt.gtol) return
c
        DO 17 i=1,N
         DG(i)=G(i)-DG(i)
 17     CONTINUE
c
        DO 19 I=1,N
          HDG(I)=0.D0
          DO 18 J=1,N
18           HDG(I)=HDG(I)+HESSIN(I,J)*DG(J)
19        CONTINUE
        FAC=0.D0
        FAE=0.D0
        SUMDG=0.0D0
        SUMXI=0.0D0
c
        DO 21 I=1,N
          FAC=FAC+DG(I)*XI(I)
          FAE=FAE+DG(I)*HDG(I)
          SUMDG=SUMDG+DG(I)**2
          SUMXI=SUMXI+XI(I)**2
21      CONTINUE
c
        if (FAC**2.gt.EPS*SUMDG*SUMXI) then
        FAC=1.D0/FAC
        FAD=1.D0/FAE
c
        DO 22 I=1,N
          DG(I)=FAC*XI(I)-FAD*HDG(I)
22      CONTINUE
c
        DO 24 I=1,N
          DO 23 J=1,N
            HESSIN(I,J)=HESSIN(I,J)+FAC*XI(I)*XI(J)
     *        -FAD*HDG(I)*HDG(J)+FAE*DG(I)*DG(J)
23        CONTINUE
24      CONTINUE
        ENDIF
c
        DO 26 I=1,N
          XI(I)=0.D0
          DO 25 J=1,N
            XI(I)=XI(I)-HESSIN(I,J)*G(J)
25        CONTINUE
26      CONTINUE
27    CONTINUE
c      PAUSE 'too many iterations in DFPMIN'
      write(*,*) 'too many iterations in DFPMIN'
      return
      
c
 450  format ('NORM OF GRADIENT ',F17.10)
 480  format ('SCF ENERGY = ',F17.10)
 500  format (I3,2x,3(F17.10,2x))
 501  format (3(F12.6,2x))
      RETURN
      END
c
c     ===========================================================
