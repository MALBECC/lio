c MAIN SUBROUTINE ----------------------------------------------------
C DFT calculation with gaussian basis sets
c---------------------------------------------------------------------

      use garcha_mod
      implicit real*8 (a-h,o-z)

      character(len=20)::argument,inpfile,inpbasis,inpcoords
      integer::charge
      logical::filexist
      REAL*8, dimension (:,:), ALLOCATABLE   :: dxyzqm
      namelist /lio/ natom,nsol,charge,OPEN,NMAX,Nunp,VCINP,frestartin, 
     > GOLD,told,rmax,rmaxs,predcoef,  
     > idip,writexyz,intsoldouble,DIIS,ndiis,dgtrig, 
     > Iexch,integ,dens,igrid,igrid2,timedep, tdstep, ntdstep, 
     > propagator,NBCH, 
     > field,a0,epsilon,exter,Fx,Fy,Fz, tdrestart, writedens
 
      integer :: ifind, ierr

     !defaults
      basis='input'  ! name of the base file
      output='output'
      fcoord='qm.xyz'
      fmulliken='mulliken'
      frestart='restart.out'
      frestartin='restart.in'
      verbose=.false.
      OPEN=.false.
      NMAX= 100
      NUNP= 0
      VCINP= .false.
      GOLD= 10.
      told=1.0D-6
      rmax=16
      rmaxs=5
      predcoef=.false.
      idip=1
      writexyz=.true.
      intsoldouble=.true.
      DIIS=.true.
      ndiis=30
      dgtrig=100.
      Iexch=9
      integ=.true.
      DENS = .true.
      IGRID = 2
      IGRID2 = 2
      timedep = 0
      tdstep = 2.D-3
      field=.false.
      a0=1000.0
      epsilon=1.D0
      Fx=0.05
      Fy=0.05
      Fz=0.05
      NBCH=10
      propagator=1
      tdrestart=.false.
      writedens=.true.

       narg=command_argument_count()

       do i=1, narg
         call get_command_argument(i,argument)
          select case(adjustl(argument))
           case("-i")
            call get_command_argument(i+1,inpfile)
           case("-b")
            call get_command_argument(i+1,basis)
           case("-c")
            call get_command_argument(i+1,inpcoords)
            case("-v")
            verbose=.true.
           case default
          end select
       enddo

       inquire(file=inpfile,exist=filexist)
      if(filexist) then
      open(unit=100,file=inpfile,iostat=ios)
      else
      write(*,*) 'input file ',adjustl(inpfile),' not found'
      stop
      endif
      read(100,nml=lio,iostat=ierr)
      if(ierr.gt.0) stop 'input error in lio namelist'
      inquire(file=inpcoords,exist=filexist)
      if(filexist) then
      open(unit=101,file=inpcoords,iostat=ios)
      else
      write(*,*) 'input file ',adjustl(inpcoords),' not found'
      stop
      endif
c        write(*,*) natom,nsol
        write(*,nml=lio)
       ntatom=natom+nsol

         ngnu=natom*ng0
         ngdnu=natom*ngd0
         ngDyn=ngnu
         ngdDyn=ngdnu
        ng3=4*ngDyn
      ng2=5*ngDyn*(ngDyn+1)/2+3*ngdDyn*(ngdDyn+1)/2+
     >           ngDyn+ngDyn*norbit+Ngrid
c      write(*,*)ng2,ngDyn,ngdDyn,norbit,Ngrid

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

       do i=1,natom
       read(101,*) iz(i),r(i,1:3)
                  rqm(i,1:3)=r(i,1:3)
c       write(*,*) iz(i),r(i,1:3)
       enddo
       do i=natom+1,ntatom
       read(101,*) pc(i),r(i,1:3)   ! o es pc(i-natom)???
c       write(*,*) pc(i),r(i,1:3)
       enddo
       r=r/0.529177D0
       rqm=rqm/0.529177D0
     
       call g2g_init()   !initialize g2g

        nqnuc=0
       do i=1,natom
         nqnuc=nqnuc+Iz(i)
       enddo

       nco=((nqnuc - charge)-Nunp)/2

c       write(*,*) 'NCO=',NCO
c       write(*,*) natom,ntatom,ngDyn,ngdDyn,ng0,ngd0
c       write(*,*) ng2,ngDyn,ngdDyn
c--------------------------------------------------------
       call drive(ng2,ngDyn,ngdDyn)
c--------------------------------------------------------
       if(OPEN) then
         call SCFOP(escf)
       else
         call SCF(escf,dipxyz)
       endif
c-------------------------------------------------------- 

       write(*,*) 'SCF ENRGY=',escf 
        
       allocate (dxyzqm(3,natom))
       dxyzqm=0.0
c       call dft_get_qm_forces(dxyzqm)
       call g2g_solve_groups(3, Exc, dxyzqm)
c       write(*,*) dxyzqm

c       do k=1,natom
c         write(*,'("fuerza",I,D,D,D)') 
c     >     k,dxyzqm(k,1),dxyzqm(k,2),dxyzqm(k,3)
c       enddo
       
       call lio_finalize()     
       end
