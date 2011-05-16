c-------------------------------------------------------------------
c Integrals subroutine -Third part
c 2 e integrals, 3 index : wavefunction and density fitting functions
c All of them are calculated
c using the Obara-Saika recursive method.
c
c
c loop over all basis functions
c now the basis is supposed to be ordered according to the type,
c all s, then all p, then all d, .....
c inside each type, are ordered in shells
c px,py,pz , dx2,dxy,dyy,dzx,dzy,dzz, .....
c
c ns ... marker for end of s
c np ... marker for end of p
c nd ... marker for end of d
c
c r(Nuc(i),j) j component of position of nucleus i , j=1,3
c Input : G ,F,  standard basis and density basis
c F comes, computed the 1 electron part, and here the
c Coulomb part is added, without storing the integrals
c Output: F updated with Coulomb part, also Coulomb energy
c F also updated with exchange correlation part, also energy
c is updated
c this subroutine calls the fitting for exchange-correlation
c-----------------------------------------------------------------
      subroutine int3lu(NORM,natom,Iz,r,Nuc,M,Mmem,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,E2,Ex,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
       use latom
c
c
      implicit real*8 (a-h,o-z)
      logical NORM,dens,OPEN,SVD,ATRHO,integ
      logical VCINP,DIRECT,EXTR,SHFT,write
      integer nopt,iconst,igrid,igrid2,ndens_local
      INCLUDE 'param'
      parameter(pi52=34.9868366552497108D0,pi=3.14159265358979312D0)
      dimension r(nt,3),nshelld(0:3),nshell(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
c
      dimension Q(3),W(3),Rc(ngd),af(ngd),FF(ngd),P(ngd)
      dimension d(ntq,ntq),Jx(ng),Ll(3)
      dimension RMM(*),X(Md,Md)
*****
*****
c scratch space
c
c auxiliars
      dimension B(ngd,3),aux(ngd)
c
*check      COMMON /TABLE/ STR(880,0:21)
      common /Sys/ SVD,iconst
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
      common /coef/ af
      common /coef2/ B
c     common /index/ iii(ng),iid(ng)
      common /Nc/ Ndens
      common /intg1/ e_(50,3),wang(50)
      common /intg2/ e_2(116,3),wang2(116),Nr(0:54),e3(194,3),wang3(194)

c
c------------------------------------------------------------------
c now 16 loops for all combinations, first 2 correspond to
c wavefunction basis, the third correspond to the density fit
c Rc(k) is constructed adding t(i,j,k)*P(i,j)
c cf(k) , variationally obtained fitting coefficient, is
c obtained by adding R(i)*G-1(i,k)
c if the t(i,j,k) were not stored, then in order to evaluate
c the corresponding part of the Fock matrix, they should be
c calculated again.
c V(i,j) obtained by adding af(k) * t(i,j,k)
c
c
c------------------------------------------------------------------
c # of calls
      Ncall=Ncall+1
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
      pi32=pi**1.50000000000000000D0
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md
c  pointers
c
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
c
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, also F later
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues ), also space used in least-squares
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
c least squares
      M17=M15+MM
c vectors of MO
      M18=M17+MMd
c
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
*  Mmem is the pointer to the RAM memory fo two-e integrals
*  Mmem = M20 in CLOSED SHELL case,  Mmem = M23 in OPEN SHELL case
c
         Mmem1=Mmem + M*(M+1)/2*Md+1
c
        NCOa=NCO
        NCOb=NCO+Nunp
c
c end ------------------------------------------------
      if (NORM) then
      sq3=dsqrt(3.D0)
      else
      sq3=1.D0
      endif

      call timer_start('principio int3lu')
c
      do 1 l=1,3
 1     Ll(l)=l*(l-1)/2
c
      do 6 k=1,Md
 6     Rc(k)=0.D0
c
******
c        write(*,*) 'cosas',M,Md,kknums,MM
        do kk=1,kknums
          iikk=(kk-1)*Md
             do k=1,Md
          Rc(k)=Rc(k)+RMM(kkind(kk))*cool(iikk+k)
                 
c               write(88,*) cool(iikk+k), iikk+k
         enddo
         enddo


c         do  kk = 1,MM
c           do  k = 1,Md
c            Rc(k) = Rc(k) + RMM(kk)*cool((kk-1)*Md+k)
c         enddo
c      enddo

       call timer_stop('principio int3lu')
*
*
c
c------------------------------------------------
c     write(*,*) 'Nro integ',ix
c--- calculation of variational coefficients
c
c calculation of fitting coefficients
c
c Constraint that integrated fitted density = N electrons
c
cSVD PART  ------------------------------------------------
c
      if (SVD) then
c
      MMp=Md*(Md+1)/2
      do 199 k=1,MMp
 199   RMM(M9+k-1)=0.0D0
      do 208 k=1,Md
       af(k)=0.0D0
 208   RMM(M9+k-1)=Rc(k)
c
       k1=0
       do 116 j=1,Md
       do 116 i=j,Md
c
       k1=k1+1
c
        X(i,j)=RMM(M7+k1-1)
        X(j,i)=X(i,j)
 116   continue
c
      M10=M9+Md
      M12=M10+Md
      Md3=3*Md

      call timer_start('dgelss')
c ESSL OPTION ------------------------
#ifdef essl
      CALL DGESVF(2,X,Md,RMM(M9),Md,1,RMM(M10),
     >             Md,Md,RMM(M12),Md3)
      imax=idamax(Md,RMM(M10),1)
      ss=RMM(M10+imax-1)
      tau=0.1192D-14*Md*ss

      CALL DGESVS(X,Md,RMM(M9),Md,1,RMM(M10),af,Md,Md,Md,tau)
#endif
c---------------------------------------
c LAPACK OPTION
#ifdef pack
      do i=1,Md
        af(i)=Rc(i)
      enddo
c
      Md5=5*Md
      rcond=1.0D-06
      call dgelss(Md,Md,1,X,Md,af,Md,RMM(M9),rcond,irank,RMM(M10),
     >            Md5,info)
c
c
#endif
      call timer_stop('dgelss')
c
c END SVD PART --
c
c if SVD.eq.false, then goes to Normal equation method, with or without
c constraint
c
      else
c---------------------------------------
c NORMAL EQUATION PART
      if (iconst.eq.1) then
c P : integrals of fitting functions --------------
c
      do 294 k=1,Md
 294   P(k)=0.0D0
c
      do 295 k=1,nsd
       do 295 nk=1,ncontd(k)
c
 295   P(k)=P(k)+cd(k,nk)/dsqrt(ad(k,nk)**3)
c
c p functions
c all integrals 0.0
c
c d functions ----
      do 297 k=nsd+npd+1,Md,6
       do 297 nk=1,ncontd(k)
c
        t0=cd(k,nk)/(dsqrt(ad(k,nk)**3)*2.0D0*ad(k,nk))
       do 297 l1=1,3
       do 297 l2=1,l1
        kk=k+Ll(l1)+l2-1
c
        if (l1.eq.l2) then
        P(kk)=P(kk)+t0/sq3
        endif
c
 297    continue
c
        do 298 k=1,Md
 298     P(k)=P(k)*pi32
c
c--------------------------------------
      do 300 m1=1,Md
       FF(m1)=0.0D0
       do 301 k=1,m1-1
 301   FF(m1)=FF(m1)+P(k)*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
       do 302 k=m1,Md
 302   FF(m1)=FF(m1)+P(k)*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
 300  continue
c
      r0=0.0D0
      r1=0.0D0
      do 900 m1=1,Md
       r0=r0+FF(m1)*Rc(m1)
       r1=r1+FF(m1)*P(m1)
 900  continue
c
      Nel=2*NCO+Nunp
      bda=(Nel-r0)/r1
c
      ss9=0.0D0
      do 200 m1=1,Md
      af(m1)=0.D0
      do 201 k=1,m1-1
 201  af(m1)=af(m1)+(Rc(k)+bda*P(k))*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
      do 202 k=m1,Md
 202  af(m1)=af(m1)+(Rc(k)+bda*P(k))*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
       ss9=ss9+af(m1)*P(m1)
 200  continue

c tests
c     do i=1,Md
c      write(20,*) i,Rc(i),af(i)
c     enddo
c     write(20,*) Nel,bda,ss9
c
c no constraint applied
      else
      do 1200 m1=1,Md
      af(m1)=0.D0
      do 1201 k=1,m1-1
 1201  af(m1)=af(m1)+Rc(k)*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
      do 1202 k=m1,Md
 1202  af(m1)=af(m1)+Rc(k)*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
 1200  continue

      endif
      endif
c----------------------------------------------------------------
c
c Initialization of Fock matrix elements
c
c
      do 215 k=1,MM
 215   RMM(M5+k-1)=RMM(M11+k-1)
c
      if (OPEN) then
      do 216 k=1,MM
 216   RMM(M3+k-1)=RMM(M11+k-1)
      endif
c----------------------------------------------------------------
c
c
c here af is ready
c
c Ndens says if the density subroutines should evaluate
c the density using the Density Matrix or the vectors
c Since Damping is applied on Density Matrix, at the beggining
c of the SCF it is more convenient to use the D.M.
c
c     if (Ncall.ge.2) then
c      Ndens=2
c     endif
c
c call fit for exchange correlation routine
      if (integ) then
c
       do i=1,Md
        B(i,2)=0.0D0
        B(i,1)=0.0D0
        B(i,3)=0.0D0
       enddo
c
       else
c
c ---- Least squares options, true recalculating functions at every
c      iteration, and false, saving on disk
c
      if (dens) then
      call exch(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,
     > Md,ncontd,nshelld,Nucd,cd,ad,RMM,NCOa,NCOb,M17,B)
      else
c
      call exch2(OPEN,Iz,natom,RMM,nshelld,M,Md,M17,NCOa,NCOb,B)
c
      endif
      endif
c------------------------------------------------------------------
c
c
      Ex=0.D0
c end of Least squares option ----------
c
c

      Ea=0.D0
      Eb=0.D0
c
      do 610 m1=1,Md
       Ex=Ex+B(m1,1)*Rc(m1)
       Ea=Ea+af(m1)*Rc(m1)
      do 611 k=1,m1
 611   Eb=Eb+af(k)*af(m1)*RMM(M7+m1+(2*Md-k)*(k-1)/2-1)
      do 612 k=m1+1,Md
 612   Eb=Eb+af(k)*af(m1)*RMM(M7+k+(2*Md-m1)*(m1-1)/2-1)
 610  continue
c
c
c------------------------------------------------------------------
c Calculation of all integrals again, for constructing the
c Fock matrix
c Previously energy was computed, and coefficients for
c the fit were generated
c------------------------------------------------------------------
c
      if (OPEN) then
       do k=1,Md
        aux(k)=af(k)+B(k,3)
       enddo
       else
       do k=1,Md
        aux(k)=0.0D0
       enddo
       endif

c
        do 217 k=1,Md
  217    af(k)=af(k)+B(k,2)
c
****
****

         call timer_start('int3lu')
          if (open) then       
         do kk=1,kknums
c             iikk=(kkind(kk)-1)*Md
              iikk=(kk-1)*Md
             do k=1,Md
        RMM(M5+kkind(kk)-1)=RMM(M5+kkind(kk)-1)+af(k)*cool(iikk+k)
        RMM(M3+kkind(kk)-1)=RMM(M3+kkind(kk)-1)+aux(k)*cool(iikk+k)
         enddo
         enddo
         else


          do kk=1,kknums
c            iikk=(kkind(kk)-1)*Md
               iikk=(kk-1)*Md
              do k=1,Md
              RMM(M5+kkind(kk)-1)=RMM(M5+kkind(kk)-1)+af(k)*cool(iikk+k)
c              write(88,*) cool(iikk+k),kk,k
              enddo
          enddo

         endif


         call timer_stop('int3lu')



****
****
c
c Numerical integration for obtaining the exchange-correlation part
c of Fock matrix and also for the exchange-correlation energy
c
      do 317 k=1,Md
 317   af(k)=af(k)-B(k,2)
c
       if (integ) then

       NCOa=NCO
       NCOb=NCO+Nunp
       write(*,*) 'exchnum int3lu'
c       write(957,*) 'int3lu'
c       call timer_start('exchfock')
c        write(*,*) 'que pasa 2?'
#ifdef G2G
       call g2g_solve_groups(0,Ex,0)
#else
       call EXCHFOCK(OPEN,NORM,natom,Iz,Nuc,ncont,nshell,a,c,r,
     >        M,M18,NCOa,NCOb,RMM,Ex)
       write(*,*) 'energia cpu',Ex
c      do kk=1,m
c        do jj=kk,m
c          write(*,*) 'rmm output',RMM(i)
c        enddo
c      enddo

#endif
c       call timer_stop('exchfock')
      
       Ndens=Ndens+1
       endif
c
c
c
      E2=Ea-Eb/2.D0+Ex
c
      return
      end
