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
      module subm_int3lu; contains
      subroutine int3lu(E2)

      use garcha_mod, only: RMM, ll, af, X, B, ngd, md, integ, M
     >                    , kknumd, kknums, MEMO, Nang, natom, NCO
     >                    , NORM, Nunp, OPEN, pi32, nshell, nshelld
     >                    , SVD, cool, cools, kkind, kkinds, ncontd
     >                    , ad, cd
      implicit none
c
      real*8, intent(inout) :: E2
      real*8 :: Q(3), W(3), Rc(Md), FF(Md), P(Md), Jx(M)
      real*8 :: aux(ngd)

      integer :: M1, M2, M3, M5, M7, M9
      integer :: M10, M11, M12, M13, M15, M17, M18, M18b, M19
      integer :: M20, M21, M22, M23
      integer :: MM, MMd, MMp, Md2, Md3, Md5
      integer :: NCOa, NCOb, Ndens, Nel, iconst, irank, info
      integer :: nd, ndd, nk, np, npd, ns, nsd
      integer :: l, l1, l2, i, j, k, k1, kk, kkk
      integer :: iikk

      real*8  :: sq3, r0, r1, rcond, rr, bda
      real*8  :: t0, ss9, Ex, Ea, Eb


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
c      Ncall=Ncall+1
      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      M2=2*M
c
      nsd=nshelld(0)
      npd=nshelld(1)
      ndd=nshelld(2)
      Md2=2*Md

c  pointers
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
        NCOa=NCO
        NCOb=NCO+Nunp

       iconst=0
       Ndens = 1
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
c         Mmem1=Mmem + M*(M+1)/2*Md+1
c end ------------------------------------------------

      if (NORM) then
        sq3=dsqrt(3.D0)
      else
        sq3=1.D0
      endif

      if (MEMO) then
        call g2g_timer_start('principio int3lu')

        do 1 l=1,3
 1        Ll(l)=l*(l-1)/2

        do 6 k=1,Md
 6        Rc(k)=0.D0

        do kk=1,kknumd !ns*(ns+1)/2+1,ns*(ns+1)/2+ns*np!kknumd
          iikk=(kk-1)*Md
          do k=1,Md
            Rc(k)=Rc(k)+RMM(kkind(kk))*cool(iikk+k)
          enddo
        enddo

        do kk=1,kknums
          iikk=(kk-1)*Md
          do k=1,Md
            Rc(k)=Rc(k)+RMM(kkinds(kk))*cools(iikk+k)
          enddo
        enddo
                
c------------------------------------------------
c calculation of variational coefficients
c calculation of fitting coefficients
c Constraint that integrated fitted density = N electrons
c------------------------------------------------
c       SVD Part
        if (SVD) then
          MMp=Md*(Md+1)/2
          do 199 k=1,MMp
 199        RMM(M9+k-1)=0.0D0
          do 208 k=1,Md
            af(k)=0.0D0
 208      RMM(M9+k-1)=Rc(k)

          k1=0
          do 116 j=1,Md
          do 116 i=j,Md
            k1=k1+1
            X(i,j)=RMM(M7+k1-1)
            X(j,i)=X(i,j)
 116      continue

          M10=M9+Md
          M12=M10+Md
          Md3=3*Md
          call g2g_timer_start('dgelss')

c        ESSL OPTION
#ifdef essl
          CALL DGESVF(2,X,Md,RMM(M9),Md,1,RMM(M10),Md,Md,RMM(M12),Md3)
          imax=idamax(Md,RMM(M10),1)
          ss=RMM(M10+imax-1)
          tau=0.1192D-14*Md*ss
          CALL DGESVS(X,Md,RMM(M9),Md,1,RMM(M10),af,Md,Md,Md,tau)
#endif
c         LAPACK OPTION
#ifdef pack
          do i=1,Md
            af(i)=Rc(i)
          enddo
          Md5=5*Md
          rcond=1.0D-06
          call dgelss(Md,Md,1,X,Md,af,Md,RMM(M9),rcond,irank,RMM(M10),
     >                Md5,info)
#endif
          call g2g_timer_stop('dgelss')
!         END SVD

!     if SVD.eq.false, then goes to Normal equation method, with or without constraint
        else

!       NORMAL EQUATION PART
          if (iconst.eq.1) then ! Constrain applied

!           P : integrals of fitting functions
            do 294 k=1,Md
 294          P(k)=0.0D0

            do 295 k=1,nsd
            do 295 nk=1,ncontd(k)
 295          P(k)=P(k)+cd(k,nk)/dsqrt(ad(k,nk)**3)

c           p and d functions
            do 297 k=nsd+npd+1,Md,6
            do 297 nk=1,ncontd(k)
              t0=cd(k,nk)/(dsqrt(ad(k,nk)**3)*2.0D0*ad(k,nk))
              do 297 l1=1,3
              do 297 l2=1,l1
                kk=k+Ll(l1)+l2-1
                if (l1.eq.l2) then
                  P(kk)=P(kk)+t0/sq3
                endif
 297        continue

            do 298 k=1,Md
 298          P(k)=P(k)*pi32

c--------------------------------------
            do 300 m1=1,Md
              FF(m1)=0.0D0
              do 301 k=1,m1-1
 301          FF(m1)=FF(m1)+P(k)*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
              do 302 k=m1,Md
 302          FF(m1)=FF(m1)+P(k)*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
 300        continue

            r0=0.0D0
            r1=0.0D0
            do 900 m1=1,Md
              r0=r0+FF(m1)*Rc(m1)
              r1=r1+FF(m1)*P(m1)
 900        continue

            Nel=2*NCO+Nunp
            bda=(Nel-r0)/r1

            ss9=0.0D0
            do 200 m1=1,Md
              af(m1)=0.D0
              do 201 k=1,m1-1
 201       af(m1)=af(m1)+(Rc(k)+bda*P(k))*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
              do 202 k=m1,Md
 202       af(m1)=af(m1)+(Rc(k)+bda*P(k))*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
              ss9=ss9+af(m1)*P(m1)
 200        continue

          else ! No constraint applied

            do 1200 m1=1,Md
               af(m1)=0.D0
            do 1201 k=1,m1-1
 1201          af(m1)=af(m1)+Rc(k)*RMM(M9+m1+(2*Md-k)*(k-1)/2-1)
            do 1202 k=m1,Md
 1202          af(m1)=af(m1)+Rc(k)*RMM(M9+k+(2*Md-m1)*(m1-1)/2-1)
 1200       continue
          endif
        endif


c-------------------------------------------------------
c       Initialization of Fock matrix elements
        do 215 k=1,MM
 215    RMM(M5+k-1)=RMM(M11+k-1)

        if (OPEN) then
        do 216 k=1,MM
 216      RMM(M3+k-1)=RMM(M11+k-1)
        endif

c Here af is ready
c Ndens says if the density subroutines should evaluate
c the density using the Density Matrix or the vectors
c Since Damping is applied on Density Matrix, at the beggining
c of the SCF it is more convenient to use the D.M.
c call fit for exchange correlation routine
        if (integ) then
          do i=1,Md
            B(i,1)=0.0D0
            B(i,2)=0.0D0
            B(i,3)=0.0D0
          enddo
        else
          stop
        endif

        Ex=0.D0
        Ea=0.D0
        Eb=0.D0

        do 610 m1=1,Md
          Ex=Ex+B(m1,1)*Rc(m1)
          Ea=Ea+af(m1)*Rc(m1)
        do 611 k=1,m1
 611      Eb=Eb+af(k)*af(m1)*RMM(M7+m1+(2*Md-k)*(k-1)/2-1)
        do 612 k=m1+1,Md
 612      Eb=Eb+af(k)*af(m1)*RMM(M7+k+(2*Md-m1)*(m1-1)/2-1)
 610    continue

c Calculation of all integrals again, for constructing the
c Fock matrix
c Previously energy was computed, and coefficients for
c the fit were generated
        if (OPEN) then
          do k=1,Md
            aux(k)=af(k)+B(k,3)
          enddo
        else
          do k=1,Md
            aux(k)=0.0D0
          enddo
        endif

        do 217 k=1,Md
  217     af(k)=af(k)+B(k,2)

        call g2g_timer_stop('principio int3lu')
        call g2g_timer_start('int3lu')
        if (open) then ! Does the single and double precision
          do kk=1,kknumd
            iikk=(kk-1)*Md
            do k=1,Md
           RMM(M5+kkind(kk)-1)=RMM(M5+kkind(kk)-1)+af(k)*cool(iikk+k)
           RMM(M3+kkind(kk)-1)=RMM(M3+kkind(kk)-1)+aux(k)*cool(iikk+k)
            enddo
          enddo
          do kk=1,kknums
            iikk=(kk-1)*Md
            do k=1,Md
        RMM(M5+kkinds(kk)-1)=RMM(M5+kkinds(kk)-1)+af(k)*cools(iikk+k)
        RMM(M3+kkinds(kk)-1)=RMM(M3+kkinds(kk)-1)+aux(k)*cools(iikk+k)
            enddo
          enddo
        else
          do kk=1,kknumd
            iikk=(kk-1)*Md
            kkk=M5+kkind(kk)-1
            do k=1,Md
              RMM(kkk)=RMM(kkk)+af(k)*cool(iikk+k)
            enddo
          enddo
          do kk=1,kknums
            iikk=(kk-1)*Md
            rr=0.
            do k=1,Md
              rr=rr+af(k)*cools(iikk+k)
            enddo
            RMM(M5+kkinds(kk)-1)=RMM(M5+kkinds(kk)-1)+rr
          enddo
        endif
        call g2g_timer_stop('int3lu')
      else
           do k=1, MM
             if (OPEN) then
               RMM(M5+k-1)=RMM(M11+k-1)
               RMM(M3+k-1)=RMM(M11+k-1)
             else
               RMM(M5+k-1)=RMM(M11+k-1)
             endif
           enddo
           Ea=0.D0
           call aint_coulomb_fock(Ea)
           Eb=0.D0

           do 6101 m1=1,Md
             do 6111 k=1,m1
 6111           Eb = Eb +af(k)*af(m1)*RMM(M7+m1+(2*Md-k)*(k-1)/2-1)
             do 6121 k=m1+1,Md
 6121           Eb = Eb +af(k)*af(m1)*RMM(M7+k+(2*Md-m1)*(m1-1)/2-1)
 6101      continue
      endif

 
c Numerical integration for obtaining the exchange-correlation part
c of Fock matrix and also for the exchange-correlation energy

      do 317 k=1,Md
        af(k)=af(k)-B(k,2)
  317 end do

      if (integ) then
        NCOa=NCO
        NCOb=NCO+Nunp
        Ndens=Ndens+1
      endif
      E2=Ea-Eb/2.D0

      return
      end subroutine int3lu
      end module subm_int3lu
