      program liosolo
c MAIN SUBROUTINE ----------------------------------------------------
C DFT calculation with gaussian basis sets
c---------------------------------------------------------------------
      use garcha_mod
      use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP
     > ,cutECP,local_nonlocal, ecp_debug,ecp_full_range_int
     > ,verbose_ECP,Cnorm,FOCK_ECP_read, FOCK_ECP_write,Fulltimer_ECP
     > ,cut2_0,cut3_0
#ifdef CUBLAS
      use cublasmath
#endif
      implicit real*8 (a-h,o-z)

      character(len=20)::argument,inpfile,inpbasis,inpcoords
      integer::charge
      logical::filexist
      REAL*8, dimension (:,:), ALLOCATABLE   :: dxyzqm,dxyzcl
      namelist /lio/ natom,nsol,charge,OPEN,NMAX,Nunp,VCINP,frestartin,
     > GOLD,told,Etold,rmax,rmaxs,predcoef,
     > idip,writexyz,intsoldouble,DIIS,ndiis,dgtrig,
     > Iexch,integ,dens,igrid,igrid2,timedep, tdstep, ntdstep,
     > propagator,NBCH,
     > field,a0,epsilon,exter,Fx,Fy,Fz, tdrestart, writedens,
     > writeforces,basis_set,fitting_set,int_basis,

!%% Effective Core Potential Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     > ecpmode,ecptypes,tipeECP,ZlistECP,
     > cutECP,local_nonlocal, ecp_debug,ecp_full_range_int,verbose_ECP,
     > verbose,FOCK_ECP_read, FOCK_ECP_write,
     > Fulltimer_ECP,cut2_0,cut3_0,
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%% Hybrid damping-diis Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     > hybrid_converg,good_cut,
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%% Output format Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     > style, allnml,
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     > cubegen_only,cube_res,
     > cube_dens,cube_dens_file,
     > cube_orb,cube_sel,cube_orb_file,cube_elec,cube_elec_file

      integer :: ifind, ierr

      ! Defaults
      narg=command_argument_count()
      call lio_defaults()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	call LIO_LOGO2()

      do i=1, narg
        call get_command_argument(i,argument)
        select case(adjustl(argument))
          case("-i")
            call get_command_argument(i+1,inpfile)
          case("-b")
            call get_command_argument(i+1,basis)
          case("-bs")
            omit_bas=.true.
            call get_command_argument(i+1,basis_set)
          case("-fs")
            call get_command_argument(i+1,fitting_set)
          case("-ib")
            int_basis=.true.
            !call get_command_argument(i+1,int_basis)
          case("-c")
            call get_command_argument(i+1,inpcoords)
          case("-v")
            verbose=.true.
          case default
        end select
      enddo


      call g2g_timer_sum_start("Total")



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


      ! Coordinate file needs to be read twice; first time to get Iz.
      allocate (iz(natom), r(1,3))
      do i=1, natom
        read(101,*) iz(i),r(1,1:3)
      enddo
      deallocate (r)

      call lio_init(natom, Iz, sol, charge)

      ! Coordinates file is read again to get all available data.
      do i=1,natom
        read(101,*) iz(i),r(i,1:3)
        rqm(i,1:3)=r(i,1:3)
      enddo
      do i=natom+1,ntatom
        read(101,*) pc(i),r(i,1:3)  
      enddo
      r=r/0.529177D0
      rqm=rqm/0.529177D0


c--------------------------------------------------------
 
       call liomain()       !no hace nada!!!!!!
       if (.not.allocated(Smat))    allocate(Smat(M,M))
       if (.not.allocated(RealRho)) allocate(RealRho(M,M))
c--------------------------------------------------------

       if(OPEN) then
         if (ecpmode) stop "Lio havent got ECP for open shell"         
         call SCFOP(escf,dipxyz)
       else
         call SCF(escf,dipxyz)
       endif
c--------------------------------------------------------

c       write(*,*) 'SCF ENRGY=',escf

      if(writeforces) then
	if (ecpmode) stop "Lio havent got ECP forces"
       open(unit=123,file='fuerzas')
       allocate (dxyzqm(3,natom))
       dxyzqm=0.0

       if(nsol.gt.0) then
          allocate (dxyzcl(3,natom+nsol))
          dxyzcl=0.
       endif

       call dft_get_qm_forces(dxyzqm)
       if (nsol.gt.0) then
         call dft_get_mm_forces(dxyzcl,dxyzqm)
       endif
c       call g2g_solve_groups(3, Exc, dxyzqm)
c       write(*,*) dxyzqm

       do k=1,natom
       write(123,100)
     >     k,dxyzqm(k,1),dxyzqm(k,2),dxyzqm(k,3)
       enddo
         if(nsol.gt.0) then
          do k=natom,natom+nsol
!         write(123,'("fuerza",I,D,D,D)')
           write(123,100)
     >     k,dxyzcl(k,1),dxyzcl(k,2),dxyzcl(k,3)
          enddo

         endif
       deallocate (dxyzqm)
       if(nsol.gt.0) deallocate(dxyzcl)
       endif
       call lio_finalize()
100    format (I5,2x,f10.6,2x,f10.6,2x,f10.6)
       end program

