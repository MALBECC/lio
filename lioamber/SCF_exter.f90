#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% SCF_EXTER.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! File SCF_exter.f90 contains SCF setup and calls when used in tantem with     !
! AMBER or GROMACS. Routines included:                                         !
! * SCF_in  (called from AMBER)                                                !
! * SCF_gro (called from GROMACS)                                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% SCF_in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs SCF setup and routine calls from AMBER.                             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      subroutine SCF_in(E,qmcoords,clcoords,nsolin,dipxyz)
      use garcha_mod, only : r, v, Em, Rm, pc, natom, ntatom, nsol, rqm, &
                             writexyz, Iz

          implicit none
          integer, intent(in) :: nsolin
          LIODBLE , intent(in) :: qmcoords(3,natom), clcoords(4,nsolin)

          LIODBLE :: E, dipxyz(3)
          integer :: i, j, n

          nsol = nsolin ; ntatom = nsol + natom ;

          if (allocated(r))  deallocate(r)
          if (allocated(v))  deallocate(v)
          if (allocated(Em)) deallocate(Em)
          if (allocated(Rm)) deallocate(Rm)
          if (allocated(pc)) deallocate(pc)

          allocate ( r(ntatom,3), v(ntatom,3), Em(ntatom), Rm(ntatom), &
                     pc(ntatom) )

          ! This section converts the coordinates array and partial charges    !
          ! array received from Gromacs into the r (all positions), rqm (QM    !
          ! region positions) and pc (MM partial charges) arrays. Also the xyz !
          ! file containing the QM region is written.                          ! 
          if(writexyz) write(18,*) natom
          !if(writexyz) write(18,*) natom
          if(writexyz) write(18,*)

          do i=1,natom
              do j=1,3
                 r(i,j)  = qmcoords(j,i) /0.529177D0
                 rqm(i,j)= qmcoords(j,i) /0.529177D0
              enddo
              if(writexyz) write(18,345) Iz(i), rqm(i,:)*0.52917725D0
          enddo

          do i=1,nsol
              n = natom + i
              pc(n) = clcoords(4,i)
              do j=1,3
                  r(n,j) = clcoords(j,i) / 0.529177D0
              enddo
              !if(writexyz) write(18,346) pc(n), r(n,:)*0.529
          enddo
          call recenter_coords(rqm, r, natom, nsol)

           ! Calls liomain, which performs common procedures and SCF.
           call liomain(E, dipxyz)

 345  format(2x, I2,    2x, 3(f10.6,2x))
 !346  format(2x, f10.6, 2x, 3(f10.6,2x))

      return
      end subroutine SCF_in
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% ehren_in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs ehrenfest setup and routine calls from AMBER.                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehren_in( qmcoords, qmvels, clcoords, nsolin, dipxyz, E)

   use garcha_mod,    only: natom, nucpos, nucvel, atom_mass, Iz
   use td_data,       only: tdstep
   use ehrensubs,     only: ehrenaux_masses
   use debug_tools,   only: Check_posvel
   use constants_mod, only: bohr
   implicit none
   integer, intent(in)    :: nsolin
   LIODBLE,  intent(in)    :: qmcoords(3,natom)
   LIODBLE,  intent(in)    :: qmvels(3,natom)
   LIODBLE,  intent(in)    :: clcoords(4,nsolin)
   LIODBLE,  intent(inout) :: dipxyz(3)
   LIODBLE,  intent(inout) :: E

   integer :: ii, kk

   if (allocated(nucpos)) deallocate(nucpos)
   if (allocated(nucvel)) deallocate(nucvel)
   allocate(nucpos(3,natom))
   allocate(nucvel(3,natom))

   do ii=1,natom
   do kk=1,3
      ! velocity units in angstrom per 1/20.455 pico-second, and must go 
      ! to atomic units
      nucpos(kk,ii) = qmcoords(kk,ii) / bohr
      nucvel(kk,ii) = qmvels(kk,ii)
      nucvel(kk,ii) = nucvel(kk,ii)*(20.455d0)
      nucvel(kk,ii) = nucvel(kk,ii)*(2.418884326505E-5)*(1.889725989d0)
   enddo
   enddo

   call Check_posvel( tdstep, natom, nucpos, nucvel, 'Check_posvel.log' )

   if (allocated(atom_mass)) deallocate(atom_mass)
   allocate(atom_mass(natom))
   call ehrenaux_masses( natom, Iz, atom_mass )
   call SCF_in(E,qmcoords,clcoords,nsolin,dipxyz)

end subroutine ehren_in

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% SCF_GRO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs SCF setup and routine calls from GROMACS.                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      subroutine SCF_gro(E, qmcoords, clcoords, clcharge, nsolin)
      use garcha_mod, only : nsol, ntatom, natom, r, v, Em, Rm, pc, rqm, Iz, &
                             writexyz

          implicit none
          integer, intent(in) :: nsolin
          LIODBLE , intent(in) :: qmcoords(3*natom), clcoords(3*nsolin), &
                                 clcharge(nsolin )

          LIODBLE  :: dipxyz(3), E
          integer :: i, j, n


          dipxyz = 0d0 ; nsol = nsolin ; ntatom = nsol + natom ;

          call g2g_timer_sum_start("Total")

          deallocate (r, v, Em, Rm, pc)
          allocate (r(ntatom,3), v(ntatom,3), Em(ntatom), Rm(ntatom), &
                    pc(ntatom))

          ! This section converts the coordinates array and partial charges    !
          ! array received from Gromacs into the r (all positions), rqm (QM    !
          ! region positions) and pc (MM partial charges) arrays. Also the xyz !
          ! file containing the QM region is written.                          !
          if(writexyz) write(18,*) natom
          !if(writexyz) write(18,*) ntatom
          if(writexyz) write(18,*)

          do i = 1, natom
              do j = 1, 3
                  r(i,j)  = qmcoords((i-1)*3+j)
                  rqm(i,j)= qmcoords((i-1)*3+j)
              enddo
              if(writexyz) write(18,345) Iz(i), rqm(i,:)*0.52917725D0
          enddo

          do i = 1, nsol
              n     = natom+i
              pc(n) = clcharge(i)
              do j = 1, 3
                  r(n, j) = clcoords((i-1)*3+j)
              enddo
              !if(writexyz) write(18,346) pc(n), r(n,:)*0.52917725D0
          enddo
          call recenter_coords(rqm, r, natom, nsol)

          ! Calls liomain, which performs common procedures and SCF.
          call liomain(E, dipxyz)

 345  format(2x, I2,    2x, 3(f10.6,2x))
! 346  format(2x, f10.6, 2x, 3(f10.6,2x))

      return
      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% SCF_hyb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs SCF & forces calculation calls from hybrid                          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine SCF_hyb(hyb_natom, mm_natom, hyb_r, E, fdummy, Iz_cl,do_SCF, do_QM_forces, do_properties, &
                   vel, do_HOPP, do_ElecInterp)
    use garcha_mod   , only : r,rqm,pc, Iz, natom, nsol, ntatom, calc_propM
    use fstsh_data   , only : FSTSH, call_number
    use ehrensubs    , only : ehrenaux_masses
    use fstsh_data   , only : FSTSH, call_number
    use fstshsubs    , only : do_electronic_interpolation
    use constants_mod, only : H_to_eV
    implicit none
    integer, intent(in) :: hyb_natom, mm_natom !number of QM and MM atoms
    LIODBLE, intent(in) :: hyb_r(3,hyb_natom+mm_natom), Iz_cl(mm_natom) !positions and charge of MM atoms
    LIODBLE, intent(out) :: E !total LIO energy
    LIODBLE, intent(out) :: fdummy(3,hyb_natom+mm_natom) !forces
    LIODBLE, intent(inout) :: vel(3,hyb_natom)
    logical, intent(in)    :: do_ElecInterp
    logical, intent(out)   :: do_HOPP
    LIODBLE, dimension(:,:), allocatable :: fa, fmm !QM and MM forces
    LIODBLE :: dipxyz(3) !dipole
    integer :: i,j !auxiliar
    logical, intent(in) :: do_SCF, do_QM_forces !SCF & forces control variable
    logical, intent(in) :: do_properties !properties control

    ! This is the second call from HYB in order to perform Electronic Interpolation
    ! in Surface Hopping. This decides whether the HOPP occurs or not.
    if ( do_ElecInterp ) then
       if ( .not. FSTSH ) stop "HYBRID wants to perform a TSH calculation, but &
                          & in LIO input FSTSH=FALSE, please turn on this variable"

       write(*,*) "Electronic Interpolation in TSH-LIO"
       do_HOPP = .false.
       call do_electronic_interpolation( vel, do_HOPP )
       if ( do_HOPP ) then
          write(*,*) "Forces in the New Surface"
          call_number = 2
       else
          ! In the second call, if there is not hopp, is no neccessary 
          ! to calculated forces
          return
       endif
    else
       ! This variable is only valid in TSH
       call_number = 1
    endif

    allocate(fa(3,hyb_natom), fmm(3,mm_natom))
    calc_propM = do_properties

    nsol = mm_natom
    ntatom = nsol + natom

    write(*,*) "atoms QM, MM, totals", natom, nsol, ntatom
    write(*,*) "doing SCF?", do_SCF
    write(*,*) "doing forces?", do_QM_forces

    E=0.d0
    fa=0.d0
    dipxyz = 0d0

    if (allocated(pc)) deallocate(pc)
    if (allocated(r)) deallocate(r)
    allocate ( pc(ntatom), r(ntatom,3) )

    do i=1, ntatom
      do j=1, 3
        r(i,j)=hyb_r(j,i)
        if (i .le. hyb_natom) rqm(i,j)=hyb_r(j,i) !positions on QM subsystem
      enddo
        if (i .le. hyb_natom) pc(i)= Iz(i) !nuclear charge
        if (i .gt. hyb_natom) pc(i) = Iz_cl(i-hyb_natom) ! MM force-field charge
    end do

    ! Variables to TSH
    if ( .not. FSTSH ) call recenter_coords(rqm, r, natom, nsol)

! Calls main procedures.
    if (do_SCF)  call liomain(E, dipxyz)
    write(*,*) "Lio  E(H)", E
    write(*,*) "Lio  E(eV)", E*H_to_eV
    fa=0.d0
    fmm=0.d0

    if (do_QM_forces) then
      call  dft_get_qm_forces(fa)
      call  dft_get_mm_forces(fmm,fa)
    end if

    fa=-1.d0*fa  ! - change gradient to forces 
    fmm=-1.d0*fmm

    do j=1, 3
      do i=1, hyb_natom
        fdummy(j,i)=fa(j,i)
      end do
      do i=1, mm_natom
        fdummy(j,hyb_natom+i)=fmm(j,i)
      end do
    end do

    deallocate(fa, fmm)
    return
end  subroutine SCF_hyb
