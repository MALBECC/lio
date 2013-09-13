      subroutine init_lio(natomin,Izin,nclatom,charge)

      use garcha_mod
c      use qmmm_module, only : qmmm_struct,qmmm_nml
      implicit real*8 (a-h,o-z)
c
c PARAMETERS - DYNAMICAL VECTOR ONLY -------------------
c
c ngDyn : number of atoms * number basis functions
c ngdDyn: number of atoms * number auxiliar functions
c Ngrid : number of grid points (LS-SCF part)
c norbit : number of MO
c
c Ngrid may be set to 0 , in the case of using Num. Integ.
c
c      parameter (ngDyn=700)
c      parameter (ngdDyn=850)


c      include 'param'
        integer , intent(in) :: charge, nclatom
       integer , intent(in)  :: natomin
       integer , intent(in)  :: Izin(natom)
c      parameter (norbit=800,Ngrid=0)

c        write(*,*) 'charge=',charge
         natom=natomin

         
c       integer, intent(in) :: Iiiz(natom)
         ntatom=natom+nclatom
         ntatom=ntatom ! the number of clasical atoms can change
         ngnu=natom*ng0
         ngdnu=natom*ngd0
         ngDyn=ngnu
         ngdDyn=ngdnu
c
        ng3=4*ngDyn
c para version en memoria
      ng2=5*ngDyn*(ngDyn+1)/2+3*ngdDyn*(ngdDyn+1)/2+
     >           ngDyn+ngDyn*norbit+Ngrid

c      write(*,*) 'ng2 en init',ng2,ngdyn,ngddyn

      allocate(X(ngDyn,ng3),XX(ngdDyn,ngdDyn))
      allocate(RMM(ng2),RMM1(ng2),RMM2(ng2), RMM3(ng2))
      
       allocate (c(ngnu,nl),a(ngnu,nl),Nuc(ngnu),ncont(ngnu)
     >  ,cx(ngdnu,nl),ax(ngdnu,nl),Nucx(ngdnu),ncontx(ngdnu)
     > ,cd(ngdnu,nl),ad(ngdnu,nl),Nucd(ngdnu),ncontd(ngdnu)
     > ,indexii(ngnu),indexiid(ngdnu))

      allocate (r(ntatom,3),v(ntatom,3),rqm(natom,3),Em(ntatom)
     >,Rm(ntatom),pc(ntatom),Iz(natom),nnat(ntatom),af(natom*ngd0),
     >  B(natom*ngd0,3))
      allocate(d(natom,natom))
         Iz=Izin

c        write(92,*) izin
c      write(*,*) (ngDyn*ng3+ngdDyn**2+ng2+ngDyn*
c     > (NgDyn+1)/2*NgdDyn)*8*1.0D-06, '  Memoria en MB'
c      write(*,*)   ng2
!#ifdef G2G
      call g2g_init()
!#endif
        nqnuc=0
       do i=1,natom
        nqnuc=nqnuc+Iz(i)
        enddo

        nco=(nqnuc - charge)/2

        write(*,*) 'NCO=',NCO
       write(*,*) natom,ntatom,ngDyn,ngdDyn,ng0,ngd0
c--------------------------------------------------------
      call drive(ng2,ngDyn,ngdDyn)

      end
c---------------------------------------------------------------------
