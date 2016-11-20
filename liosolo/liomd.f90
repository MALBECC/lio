! Make an NVE ensamble MD.
      program liomd

      use garcha_mod
      use liosubs
#ifdef CUBLAS
      use cublasmath
#endif

      implicit none
      real*8              :: escf, dipxyz(3), time_step, Kenergy, Uenergy
      integer             :: ii, kk, narg, nqnuc, ios, nn, md_steps, charge
      logical             :: filexist

      character(len=20)   :: argument,inpfile,inpbasis,inpcoords
      real*8, allocatable :: oldpos(:,:), nowpos(:,:), newpos(:,:),  &
                             nowvel(:,:), nowfrc(:,:), atom_mass(:), &
                             dxyzqm(:,:), dxyzcl(:,:)

      namelist /lio/ natom, nsol, charge, OPEN, NMAX, Nunp, VCINP, frestartin, &
                     GOLD, told, rmax, rmaxs, predcoef, idip, writexyz, DIIS,  &
                     intsoldouble, ndiis, dgtrig, Iexch, integ, dens, igrid,   &
                     igrid2, timedep, tdstep, ntdstep, propagator, NBCH, field,&
                     a0, epsilon, exter, Fx, Fy, Fz, tdrestart, writedens,     &
                     writeforces, basis_set, fitting_set, int_basis, cube_res, &
                     cubegen_only, cube_dens, cube_dens_file, cube_orb,        &
                     cube_sel,cube_orb_file, cube_elec, cube_elec_file, frestart

      ! Loads default variable values.
      call lio_defaults()

      ! Reads command line arguments for LIO.
      narg=command_argument_count()
      do ii=1, narg
          call get_command_argument(ii,argument)
          select case(adjustl(argument))
              case("-i")
                  call get_command_argument(ii+1,inpfile)
              case("-b")
                  call get_command_argument(ii+1,basis)
              case("-bs")
                 omit_bas=.true.
                 call get_command_argument(ii+1,basis_set)
              case("-fs")
                 call get_command_argument(ii+1,fitting_set)
              case("-ib")
                 int_basis=.true.
              case("-c")
                 call get_command_argument(ii+1,inpcoords)
              case("-v")
                 verbose=.true.
              case default
          endselect
      enddo

      ! Checks input files' existence.
      call g2g_timer_sum_start("Total")

      inquire(file=inpfile,exist=filexist)
      if(filexist) then
          open(unit=100,file=inpfile,iostat=ios)
      else
          write(*,*) 'input file ',adjustl(inpfile),' not found'
          stop
      endif
      read(100,nml=lio,iostat=ios)
      if(ios.gt.0) stop 'input error in lio namelist'

      inquire(file=inpcoords,exist=filexist)
      if(filexist) then
          open(unit=101,file=inpcoords,iostat=ios)
      else
          write(*,*) 'input file ',adjustl(inpcoords),' not found'
          stop
      endif

      ! Reads coordinates file.
      ntatom = natom + nsol
      allocate (iz(natom), r(ntatom,3), rqm(natom,3), pc(ntatom))
      do ii=1,natom
        read(101,*) iz(ii),r(ii,1:3)
        rqm(ii,1:3)=r(ii,1:3)
      enddo
      do ii=natom+1,ntatom
        read(101,*) pc(ii),r(ii,1:3)
      enddo
      r=r/0.529177D0
      rqm=rqm/0.529177D0

      call lio_init(natom, Iz, nsol, charge, 0)

      call liomain()

! ALLOCATION AND INITIALIZATION
!--------------------------------------------------------------------!
      if (.not.allocated(Smat))    allocate(Smat(M,M))
      if (.not.allocated(RealRho)) allocate(RealRho(M,M))
      if (allocated(dxyzqm)) deallocate (dxyzqm)
      allocate(dxyzqm(3,natom), oldpos(natom,3), nowpos(natom,3), &
               newpos(natom,3), nowvel(natom,3), nowfrc(natom,3))

      if (allocated(atom_mass)) deallocate(atom_mass)
      allocate(atom_mass(natom))
      call set_masses(natom,Iz,atom_mass)

! SETUP
!--------------------------------------------------------------------!
      do kk=1,3
          do nn=1,natom
              oldpos(nn,kk)=r(nn,kk)
              nowpos(nn,kk)=r(nn,kk)
              newpos(nn,kk)=r(nn,kk)
              nowvel(nn,kk)=0.0d0
              nowfrc(nn,kk)=0.0d0
          enddo
      enddo

! RUN DYNAMIC
!--------------------------------------------------------------------!
      md_steps  = 0
      time_step = 0.0001                ! time in ps
      time_step = time_step * 41341.0d0 ! time from ps to au
      open(unit=501, file='liomd-trays.xyz' )
      open(unit=502, file='liomd-energy.dat')
      call write_energy(-1.0d0, escf, Kenergy, escf+Kenergy, 502)

      do nn=1,md_steps+1
        call SCF(escf,dipxyz)
        Uenergy=escf
        call dft_get_qm_forces(dxyzqm)
        do kk=1,3
        do ii=1,natom
          nowfrc(ii,kk)=-dxyzqm(kk,ii)
        enddo
        enddo
        call nuclear_verlet(natom,time_step,atom_mass,nowfrc, &
                            oldpos,nowpos,newpos,nowvel,Kenergy)
        call write_geom(natom,Iz,nowpos,501)
        call write_energy(nn*tdstep,escf,Kenergy,escf+Kenergy,502)

        oldpos=nowpos
        nowpos=newpos
        do kk=1,3
        do ii=1,natom
          r(ii,kk)   = nowpos(ii,kk)
          rqm(ii,kk) = nowpos(ii,kk)
        enddo
        enddo
      enddo

      write(*,*) 'SCF ENRGY=',escf

      if(writeforces) then
          open(unit=123,file='forces')

          if (allocated(dxyzqm)) deallocate(dxyzqm)
          allocate ( dxyzqm(3, natom) )
          dxyzqm = 0.0

          if(nsol.gt.0) then
              allocate ( dxyzcl(3, natom+nsol) )
              dxyzcl = 0.0
          endif

          call dft_get_qm_forces(dxyzqm)
          if (nsol.gt.0) then
              call dft_get_mm_forces(dxyzcl, dxyzqm)
          endif

          do kk=1,natom
              write(123,100) kk, dxyzqm(kk,1), dxyzqm(kk,2), dxyzqm(kk,3)
          enddo

          if(nsol.gt.0) then
              do kk=natom,natom+nsol
                  write(123,100) kk, dxyzcl(kk,1), dxyzcl(kk,2), dxyzcl(kk,3)
              enddo
          endif

          deallocate (dxyzqm)
          if(nsol.gt.0) deallocate (dxyzcl)
       endif

       call lio_finalize()

100    format (I5,2x,f10.6,2x,f10.6,2x,f10.6)
       end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
