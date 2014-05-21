c SCF subroutine ----------------------------------
c DIRECT VERSION
c calls all integrals generator subroutines : 1 el integrals,
c 2 el integrals, exchange fitting , so it gets S matrix, F matrix
c and P matrix in lower storage mode ( symmetric matrices)
c
c Dario Estrin, 1992
c---------------------------------------------------
      subroutine SCFOP(E,dipxyz)
      use garcha_mod
      REAL*8:: En,E2,E,Es,Ex,Exc

      call g2g_timer_start('SCF')
      write(*,*) '======>>>> INGRESO A SCFop <<<<=========='
c------------------------------------------------------------------
c
c Pointers
c
      E=0.0D0
      E1=0.0D0
      En=0.0D0
      E2=0.0D0
      Es=0.0D0

c  QUE ES ngeo ?????      
      ngeo=ngeo+1
c Number of electrons
      Nel=2*NCO+Nunp
c first P
      MM=M*(M+1)/2 
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M
      M2a=3*M
c
c Number of OM up
      NCOa=NCO
c Number of OM down
      NCOb=NCO+Nunp
c------------------------------------------------
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
c W ( eigenvalues of alpha spin  in open shell case)
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
c eigenvalues (beta spin in open shell case)
      M22=M21+MM
c
      M23 = M22 +  M  

c RAM storage of two-electron integrals (if MEMO=T)
c      M20=M19+natom*50*Nang
c------------------------------------------------
c Initializations/Defaults
c
      good=1.00D0
      niter=0
      D1=1.D0
      D2=1.D0
      DAMP0=GOLD
      DAMP=DAMP0
c
c      Qc=0.0D0
c      do i=1,natom
c        Qc=Qc+Iz(i)
c      enddo
c      Qc=Qc-Nel
c      Qc2=Qc**2
C----------------------------------------
c Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano
c        write(*,*) 'que pasa?'
  
      do i=1,natom
        natomc(i)=0
        do j=1,natom
          d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >        (r(i,3)-r(j,3))**2
          zij=atmin(i)+atmin(j)
          ti=atmin(i)/zij
          tj=atmin(j)/zij
          alf=atmin(i)*tj
          rexp=alf*d(i,j)
          if (rexp.lt.rmax) then
            natomc(i)=natomc(i)+1
            jatc(natomc(i),i)=j
          endif 
        enddo
      enddo

       do iij=nshell(0),1,-1
         nnps(nuc(iij))=iij
       enddo
       do iik=nshell(0)+nshell(1),nshell(0)+1,-1
         nnpp(nuc(iik))=iik
       enddo

       do iikk=M,nshell(0)+nshell(1)+1,-1
         nnpd(nuc(iikk))=iikk
       enddo
      
      call g2g_reload_atom_positions(igrid2)
c
c H H core, 1 electron matrix elements
c
      call int1(En)
      
c
c -- SOLVENT CASE --------------------------------------
c      if (sol) then
c      call intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,ncont,nshell,
c     >            c,a,pc,RMM,E1s)
c      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
c      endif
c-------------------------------------------------------
c E1 Energia monoelectronica
c
      E1=0.D0
      do k=1,MM
        E1=E1+RMM(k)*RMM(M11+k-1)
      enddo
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
c X transformation matrix , canonical orthogonalization
c LINEAR DEPENDENCY ELIMINATION

c
      do j=1,M
        if (RMM(M13+j-1).lt.1.0D-06) then
          write(*,*) 'LINEAR DEPENDENCY DETECTED ACA !!!!'
          do i=1,M
            X(i,j)=0.0D0
          enddo
        else
          do i=1,M
            X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))
          enddo
        endif
      enddo
c QUE ES ESTO ????
      do i=1,M
        do j=1,M
          X(i,M2a+j)=X(i,j)
        enddo
      enddo
c
c ======>>>>>> CASE OF NO STARTING GUESS PROVIDED,  <<<<<=========
c   1 E FOCK MATRIX USED
c
      if((.not.ATRHO).and.(.not.VCINP).and.primera) then
        primera=.false.
c QUE HACE ACA ?????
        do i=1,M
          do j=1,M
            X(i,M+j)=0.D0
            do k=1,j
              X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+j+(M2-k)*(k-1)/2-1)
            enddo
c
            do k=j+1,M
              X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+k+(M2-j)*(j-1)/2-1)
            enddo
c
          enddo
        enddo

c QUE HACE ACA ?????
        kk=0
        do j=1,M
          do i=j,M
            kk=kk+1
            RMM(M5+kk-1)=0.D0
            do k=1,M
              RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
            enddo
          enddo
        enddo
        
c
c diagonalization now
c
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M13+i-1)=0.D0
        enddo
c
c ESSL OPTION
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c LAPACK OPTION -----------------------------------------
#ifdef pack
c
c QUE VaRIABleS ENTRAN Y CUALES SALEN ?????
c
        call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
c-----------------------------------------------------------
        do i=1,M
          do j=1,M
            X(i,M2+j)=0.D0
            do k=1,M
c calculando  C=XC   
              X(i,M2+j)=X(i,M2+j) + X(i,k)*X(k,M+j)
            enddo
          enddo
        enddo
c
c Density Matrix
c alpha and beta coefficients set equal
c
        kk=0
        do k=1,NCOa
          do i=1,M
            kk=kk+1
            RMM(M18+kk-1)=X(i,M2+k)
          enddo
        enddo
c
        kk=0
        do k=1,NCOb
          do i=1,M
            kk=kk+1
            RMM(M18b+kk-1)=X(i,M2+k)
          enddo
        enddo
c
        kk=0
        do j=1,M
          do i=j,M
c
            kk=kk+1
            RMM(kk)=0.D0
            rhoalpha(kk)=0.D0
            rhobeta(kk)=0.D0
c
            do k=1,NCOa
              RMM(kk)=RMM(kk)+X(i,M2+k)*X(j,M2+k)
              rhoalpha(kk)=rhoalpha(kk)+X(i,M2+k)*X(j,M2+k)
            enddo
c
            do k=1,NCOb
              RMM(kk)=RMM(kk)+X(i,M2+k)*X(j,M2+k)
              rhobeta(kk)=rhobeta(kk)+X(i,M2+k)*X(j,M2+k)
            enddo
c
            if (i.ne.j) then
              RMM(kk)=2.0D0*RMM(kk)
              rhoalpha(kk)=2.0D0*rhoalpha(kk)
              rhobeta(kk)=2.0D0*rhobeta(kk)
            endif
c
          enddo
        enddo
c
c
#ifdef PRINT_MATRICES
c------ IMPRIMIENDO DENSIDADES ---------------------------------
        kk=0
        do i=1,M
          do j=i,M
            kk=kk+1 
            write(*,'(I,X,I,X,I,X,F8.5,X,F8.5,X,F8.5)') 
     <       kk,i,j,RMM(kk),rhoalpha(kk),rhobeta(kk)
          enddo
        enddo
#endif
      endif
c
c End of Starting guess (No MO , AO known)-------------------------------
c------------------------------------------------------------------------
c
      call int22()
c
**
      if (MEMO) then
         call g2g_timer_start('int3mem')
         call int3mem() 
         call int3mems()
         call g2g_timer_stop('int3mem')
      endif
****
c---------------------------------------------------------------------
c CASE  SAVING BASIS FUNCTIONS ON DISK
c ONLY IN LEAST SQUARES
c      if (.not.integ) then
c      if (dens) then
c      write(*,*) 'in nwrite'
c      call nwrite(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,NCOa,
c     >           NCOb,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)
c
c      else
c      write(*,*) 'in write'
c      call write(OPEN,NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,
c     >    NCOa,NCOb,Nucd,Md,ncontd,nshelld,cd,ad,M17,RMM)

c      endif
c      endif
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c
c      write(*,*) 'empiezo el loop=========>>>>>>',NMAX

      do 999 while (good.ge.told.and.niter.le.NMAX)
        call g2g_timer_start('Total iter')
        niter=niter+1

c        kk=0
c        do i=1,M
c          do j=i,M
c            kk=kk+1
c            write(*,'(I,X,I,X,I,X,F8.5,X,F8.5,X,F8.5)')
c     <       kk,i,j,RMM(kk),rhoalpha(kk),rhobeta(kk)
c          enddo
c        enddo


c ====>>>>> CALL INT3LU <<<<<=======
        call int3lu(E2)

c        do i=1,mm
c          write(84,*) 'rmm input',RMM(m5+i-1),RMM(m3+i-1)
c        enddo


        call g2g_solve_groups(0,Ex,0)

c        do i=1,mm
c          write(83,*) 'rmm output',RMM(m5+i-1),RMM(m3+i-1)
c       enddo

c        write(*,*) 'Ex=',Ex
c-------------------------------------------------------
        E1=0.0D0
c
c REACTION FIELD CASE --------------------------------------------
c
        call g2g_timer_start('actualiza rmm')
c----------------------------------------------------------------
c E1 includes the solvent electrostatic matrix elements (sol T)
        do k=1,MM
          E1=E1+RMM(k)*RMM(M11+k-1)
        enddo
c
c now, we know S matrix, and F matrix, and E for a given P
c 1) diagonalize S, get X=U s^(-1/2)
c 2) get U F Ut
c 3) diagonalize F
c 4) get vec ( coeff) ---->  P new
c 5) iterate
c
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
        do k=1,MM
          RMM(M5+k-1)=(RMM(M5+k-1) + DAMP*RMM(M20+k-1))/(1.D0+DAMP)
        enddo
c
c Damping  Beta Matrix
c
        do k=1,MM
          RMM(M3+k-1)=(RMM(M3+k-1) + DAMP*RMM(M21+k-1))/(1.D0+DAMP)
        enddo
c
c the newly constructed damped matrix is stored, for next iteration
c in RMM(M20) and RMM(M21) for alpha and beta spin respectively.
c
         write(*,*)
         do k=1,MM
           RMM(M20+k-1)=RMM(M5+k-1)
           RMM(M21+k-1)=RMM(M3+k-1)
c           write(*,*) RMM(M5+k-1),RMM(M3+k-1)
         enddo

c
c alpha case -------
c
         do i=1,M
           do j=1,M
             X(i,M+j)=0.D0
             do k=1,j
               X(i,M+j)=X(i,M+j) + X(k,i)*RMM(M5+j+(M2-k)*(k-1)/2-1)
             enddo
             do k=j+1,M
               X(i,M+j)=X(i,M+j) + X(k,i)*RMM(M5+k+(M2-j)*(j-1)/2-1)
             enddo
           enddo
         enddo
c
         kk=0
         do j=1,M
           do i=j,M
             kk=kk+1
             RMM(M5+kk-1)=0.D0
             do k=1,M
               RMM(M5+kk-1)=RMM(M5+kk-1) + X(i,M+k)*X(k,j)
             enddo
           enddo
         enddo
c
c now F contains transformed F
c diagonalization now
c
         do i=1,M
           RMM(M15+i-1)=0.D0
           RMM(M13+i-1)=0.D0
         enddo
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
        call g2g_timer_start('dspev')
         
c         write(*,*)"alpha antes"
c         do k=1,MM
c           write(*,*) RMM(M5+k-1),RMM(M3+k-1)
c         enddo


c ESSL OPTION ------------------------------------------------
#ifdef essl
        call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#endif
c
c LAPACK OPTION -----------------------------------------
#ifdef pack
        call dspev('V','L',M,RMM(M5),RMM(M13),X(1,M+1),M,RMM(M15),info)
#endif
c         write(*,*)"alpha despues"
c         do k=1,MM
c           write(*,*) RMM(M5+k-1),RMM(M3+k-1)
c         enddo


        call g2g_timer_stop('dspev')
        call g2g_timer_start('coeff')
c-----------------------------------------------------------
c
c diagonalization now
c
c
c new coefficients
c
        do i=1,M
          do j=1,M
c
            X(i,M2+j)=0.D0
            do k=1,M
              X(i,M2+j)=X(i,M2+j) + X(i,k)*X(k,M+j)
            enddo
          enddo
        enddo
c
c
        kk=0
        do k=1,NCOa
          do i=1,M
            kk=kk+1
            RMM(M18+kk-1)=X(i,M2+k)
          enddo
        enddo
c
c xxxxxxx aca poner que escriba ------------------------------
        if ((good.le.5.0D0*told.and.nopt.eq.0).or.niter.eq.nmax) then
          do l=1,M
            do n=1,M
              X(indexii(l),M+n)=X(l,M2+n)
            enddo
          enddo

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
c          if (niter.ge.2) then
c Level Shifting
            do i=1,M
              do j=1,M
                X(i,j)=X(i,M2+j)
              enddo
            enddo
c         endif
        endif
c
c Beta case -------------------
c
c
        do i=1,M
          do j=1,M
            X(i,M+j)=0.D0
            do k=1,j
              X(i,M+j)=X(i,M+j)+X(k,M2a+i)*RMM(M3+j+(M2-k)*(k-1)/2-1)
            enddo
c
            do k=j+1,M
              X(i,M+j)=X(i,M+j)+X(k,M2a+i)*RMM(M3+k+(M2-j)*(j-1)/2-1)
            enddo
          enddo
        enddo
c
        kk=0
        do j=1,M
          do i=j,M
            kk=kk+1
            RMM(M3+kk-1)=0.D0
            do k=1,M
              RMM(M3+kk-1)=RMM(M3+kk-1)+X(i,M+k)*X(k,M2a+j)
            enddo
          enddo
        enddo
c now F contains transformed F
c
c diagonalization now
c
        do i=1,M
          RMM(M15+i-1)=0.D0
          RMM(M22+i-1)=0.D0
        enddo
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
c         write(*,*)"beta antes"
c         do k=1,MM
c           write(*,*) RMM(M5+k-1),RMM(M3+k-1)
c         enddo

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
c         write(*,*)"beta despues"
c         do k=1,MM
c           write(*,*) RMM(M5+k-1),RMM(M3+k-1)
c         enddo

c
c new coefficients
c
       do i=1,M
         do j=1,M
           X(i,M2+j)=0.D0
           do k=1,M
             X(i,M2+j)=X(i,M2+j)+X(i,M2a+k)*X(k,M+j)
           enddo
         enddo
       enddo

c
c xxxxxxx  aca poner que escriba -------------------------------
       if ((good.le.5.0D0*told.and.nopt.eq.0).or.niter.eq.nmax) then
         do l=1,M
           do n=1,M
             X(indexii(l),M+n)=X(l,M2+n)
           enddo
         enddo
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
       do k=1,NCOb
         do i=1,M
           kk=kk+1
           RMM(M18b+kk-1)=X(i,M2+k)
c           write(*,*) X(i,M2+k)
         enddo
       enddo
c
       if (SHFT) then
c     if (niter.ge.2) then
c Level Shifting
          do i=1,M
            do j=1,M
              X(i,M2a+j)=X(i,M2+j) 
            enddo
          enddo
        endif
c     endif

      call g2g_timer_stop('coeff') 

c
c-----------------------------------------
c Construction of new density matrix and comparison with old one
c

      kk=0
      good=0.0D0
      do j=1,M
        do i=j,M
          kk=kk+1
          tmp=RMM(kk)
          RMM(kk)=0.D0
          rhoalpha(kk)=0.D0
          rhobeta(kk)=0.D0
c
          do k=1,NCOa
            k0=M18+M*(k-1)-1
            ki=k0+i
            kj=k0+j
            RMM(kk)=RMM(kk)+RMM(ki)*RMM(kj)
c rhoalpha(M*(M+1)/2)
            rhoalpha(kk)=rhoalpha(kk)+RMM(ki)*RMM(kj)
          enddo
c
          do k=1,NCOb
            k0=M18b+M*(k-1)-1
            ki=k0+i
            kj=k0+j
            RMM(kk)=RMM(kk)+RMM(ki)*RMM(kj)
c rhobeta(M*(M+1)/2) 
            rhobeta(kk)=rhobeta(kk)+RMM(ki)*RMM(kj)
          enddo
c
          if (i.ne.j) then
            RMM(kk)=2.0D0*RMM(kk)
            rhoalpha(kk)=2.0D0*rhoalpha(kk)
            rhobeta(kk)=2.0D0*rhobeta(kk)
          endif
c
          del=RMM(kk)-tmp
          if (i.ne.j) then
            del=del*sqrt(2.D0)
          endif
          good=good+del**2
        enddo
      enddo
c
#ifdef PRINT_MATRICES
c------ IMPRIMIENDO DENSIDADES ---------------------------------
        kk=0 
        do j=1,M
          do i=j,M
            kk=kk+1 
            write(*,'(I,X,I,X,I,X,E12.5,X,E12.5,X,E12.5,X,E12.5)') 
     <       kk,j,i,RMM(kk),rhoalpha(kk),rhobeta(kk),
     <       (rhoalpha(kk)+rhobeta(kk))
          enddo
        enddo
#endif
        good=sqrt(good)/float(M)
c
c--- Damping factor update - 
        DAMP=DAMP0
c        IDAMP=0
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
c          if (sol) then
            E=E+Es
c          endif
c
          D2=D1
          D1=(E-E0)
c
          E0=E
          DAMP0=DAMP
        endif
c
        E=E1+E2+En+Ex
c        if (sol) then
        E=E+Es
c        endif

         write(*,300) niter,DAMP,E
c         write(*,*) E1,E2,En,Ex
         write(*,"(A,X,F13.6,X,F13.6,X,F13.6,X,F13.6,X,F13.6)")
     >         'En,E1,E2,Ex,Etotal',En,E1,E2,Ex,E1+E2+En+Ex

c
        call g2g_timer_stop('Total iter')

c      if (nopt.le.2) then
c       write(*,300) niter,DAMP,E
c      endif
c
c      write(*,*) 'Coulomb E',E2-Ex,Ex
c        if (write1) then
c
c          open(unit=3,file='restart')
c outputs final  MO ---------------------
c alpha
c          dol=1,M
c            do n=1,NCOa
c              kk=M18+(l-1)+M*(n-1)
c              X(indexii(l),M+n)=RMM(kk)
c            enddo
c          enddo
c
c          do l=1,M
c            write(3,400) (X(l,M+n),n=1,NCOa)
c          enddo
c
c beta
c          do l=1,M
c            do n=1,NCOb
c              kk=M18b+(l-1)+M*(n-1)
c              X(indexii(l),M+n)=RMM(kk)
c            enddo
c          enddo
c
c          do l=1,M
c            write(3,400) (X(l,M+n),n=1,NCOb)
c          enddo
c
c          write(3,*) niter,E
c          close(3)
c        endif
c
 999  continue
c 995   continue
c -- SOLVENT CASE --------------------------------------
c       if (sol) then
c         call intsol(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,ncont,nshell,
c     >            c,a,pc,RMM,E1s)
c         call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
c         Es=Es+E1s
c       endif
c
c       if (GRAD) then
c
c         if (nopt.eq.0) then
c           write(*,*)
c           write(*,600)
c           write(*,610)
c           write(*,620) E1,E2-Ex,En
c         endif
c         if (sol) then
c           write(*,615)
c           write(*,625) Es
c         endif
c         call exchnumop(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
c     >                   M18,NCOa,NCOb,Exc,nopt)
c
#ifdef G2G
#ifdef ULTIMA_CPU
       call exchnum(NORM,natom,r,Iz,Nuc,M,ncont,nshell,c,a,RMM,
     >              M18,NCO,Exc,nopt)
#else
       call g2g_new_grid(igrid)

       call g2g_solve_groups(1, Exc, 0)

c       write(*,*) 'g2g-Exc',Exc
#endif
#else
#ifdef ULTIMA_G2G
c       call g2g_new_grid(igrid)
c       call g2g_solve_groups(1, Exc, 0)
#else
#endif       
#endif
        E=E1+E2+En+Es+Ens+Exc

c         E=E+Exc-Ex
c         write(*,*)
c         write(*,450) E
c       else
c         E=E-Ex
c       endif
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
      do j=1,M
        do i=j,M
          kk=kk+1
          RMM(kk)=0.D0
c
          do k=1,NCOa
            k0=M18+M*(k-1)-1
            ki=k0+i
            kj=k0+j
            RMM(kk)=RMM(kk)-RMM(M13+k-1)*RMM(ki)*RMM(kj)
          enddo
c
          do k=1,NCOb
            k0=M18b+M*(k-1)-1
            ki=k0+i
            kj=k0+j
            RMM(kk)=RMM(kk)-RMM(M22+k-1)*RMM(ki)*RMM(kj)
          enddo
c
          if (i.ne.j) then
            RMM(kk)=2.0D0*RMM(kk)
          endif
c
        enddo
      enddo

c
c-----------------------------------------------------------------
c PROPERTIES CALCULATION
c calculates dipole moment
c
c      if (idip.eq.1) then
c        call dip(ux,uy,uz)
c       u=sqrt(ux**2+uy**2+uz**2)
c
c      write(*,*)
c      write(*,*) 'DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
c       write(*,900) ux,uy,uz,u
c      write(*,*)
c u in Debyes
c      endif
c
c calculates Mulliken poputations
c       if (ipop.eq.1) then
c        call int1(En)
c
c
c        do n=1,natom
c          q(n)=Iz(n)
c        enddo
c
c        do i=1,M
c
c          do j=1,i-1
c            kk=i+(M2-j)*(j-1)/2
c            t0=RMM(kk)*RMM(+kk-1)/2.D0
c            q(Nuc(i))=q(Nuc(i))-t0
c          enddo
c
c          kk=i+(M2-i)*(i-1)/2
c          t0=RMM(kk)*RMM(M5+kk-1)
c          q(Nuc(i))=q(Nuc(i))-t0
c
c          do j=i+1,M
c            kk=j+(M2-i)*(i-1)/2
c            t0=RMM(kk)*RMM(M5+kk-1)/2.D0
c            q(Nuc(i))=q(Nuc(i))-t0
c          enddo
c        enddo
c
c        write(*,*) 'MULLIKEN POPULATION ANALYSIS'
c        write(*,770)

c        do n=1,natom
c         write(*,760) n,Iz(n),q(n)
c        enddo
c
c        write(*,*)
c
c UNPAIRED SPIN POPULATION
c M7 spin alpha density matrix, M11 spin beta density matrix
c        kk=M7-1
c        kk1=M11-1
c        do j=1,M
c          do i=j,M
c            kk=kk+1
c            kk1=kk1+1
c            RMM(kk)=0.0D0
c            RMM(kk1)=0.0D0
c
c            do k=1,NCOa
c              k0=M18+M*(k-1)-1
c              ki=k0+i
c              kj=k0+j
c              RMM(kk)=RMM(kk) + RMM(ki)*RMM(kj)
c            enddo
c
c            do k=1,NCOb
c              k0=M18b+M*(k-1)-1
c              ki=k0+i
c              kj=k0+j
c              RMM(kk1)=RMM(kk1)  + RMM(ki)*RMM(kj)
c            enddo
c
c            if (i.ne.j) then
c              RMM(kk)=2.0D0*RMM(kk)
c              RMM(kk1)=2.0D0*RMM(kk1)
c            endif
c
c          enddo
c        enddo
       
c        do n=1,natom
c          q(n)=0.0D0
c        enddo
c
c        do i=1,M
c
c          do j=1,i-1
c            kk=i+(M2-j)*(j-1)/2
c            kka=kk+M7-1
c            kkb=kk+M11-1
c         
c            t0=(RMM(kka)-RMM(kkb))*RMM(M5+kk-1)/2.D0
c            q(Nuc(i))=q(Nuc(i))-t0
c          enddo
c
c          kk=i+(M2-i)*(i-1)/2
c          kka=kk+M7-1
c          kkb=kk+M11-1
c          t0=(RMM(kka)-RMM(kkb))*RMM(M5+kk-1)
c          q(Nuc(i))=q(Nuc(i))-t0
c
c          do j=i+1,M
c            kk=j+(M2-i)*(i-1)/2
c            kka=kk+M7-1
c            kkb=kk+M11-1
c
c            t0=(RMM(kka)-RMM(kkb))*RMM(M5+kk-1)/2.D0
c            q(Nuc(i))=q(Nuc(i))-t0
c          enddo
c        enddo
c
c        write(*,*) 'UNPAIRED SPIN MULLIKEN POPULATION ANALYSIS'
c        write(*,770)
c
c        do n=1,natom
c          write(*,760) n,Iz(n),q(n)
c        enddo
c
c      endif
c
c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
c
c        if (icharge.eq.1) then
c          Q1=-(2*NCO+Nunp)
c         do n=1,natom
c          Q1=Q1+Iz(n)
c         enddo
c         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c     >            c,a,RMM,map,Q1)
c        endif
c-----------------------------------------------------
c      do l=1,M
c        do n=1,NCOa
c          kk=M18+(l-1)+M*(n-1)
c          X(index(l),M+n)=RMM(kk)
c        enddo
c      enddo
c
c      do l=1,M
c        write(2,400) (X(l,M+n),n=1,NCOa)
c      enddo
c
c-------------------------------------------------
c       do l=1,M
c         do n=1,NCOb
c           kk=M18b+(l-1)+M*(n-1)
c           X(index(l),M+n)=RMM(kk)
c         enddo
c       enddo
c
c      do l=1,M
c        write(2,400) (X(l,M+n),n=1,NCOb)
c      enddo
c-------------------------------------------------
c

         deallocate (kkind,kkinds)
         deallocate(cool,cools)

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
      call g2g_timer_stop('SCF')
       return
       end
C  -------------------------                                            
