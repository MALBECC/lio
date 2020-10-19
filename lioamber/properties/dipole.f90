!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DIPOLE.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine dipole(uDip, Pmat_v, nElec, at_pos, at_dists, atom_z, mm_charges, &
                  print_dip)
   use faint_cpu , only: intdip

   implicit none
   integer, intent(in)           :: nElec
   integer, intent(in)           :: atom_z(:)
   logical, intent(in), optional :: print_dip
   LIODBLE, intent(in)           :: Pmat_v(:)
   LIODBLE, intent(in)           :: at_pos(:,:)
   LIODBLE, intent(in)           :: at_dists(:,:)
   LIODBLE, intent(in)           :: mm_charges(:)
   LIODBLE, intent(inout)        :: uDip(3)

   LIODBLE, allocatable :: dip_mat(:,:), uDipAt(:)
   LIODBLE :: Qc, factor
   integer :: iCount
   
   allocate(dip_mat(3, size(Pmat_v,1)), uDipAt(3))
   dip_mat = 0.0D0

   call intdip(dip_mat, at_pos, at_dists, Pmat_v)

   uDip = 0.0D0
   do iCount = 1, size(Pmat_v,1)
      uDip(1) = uDip(1) + dip_mat(1, iCount) 
      uDip(2) = uDip(2) + dip_mat(2, iCount) 
      uDip(3) = uDip(3) + dip_mat(3, iCount) 
   enddo

   uDipAt = 0.0D0
   Qc     = 0.0D0
   do iCount = 1, size(atom_z,1)
      Qc  = Qc  + atom_z(iCount)
      uDipAt(1) = uDipAt(1) + atom_z(iCount) * at_pos(iCount,1)
      uDipAt(2) = uDipAt(2) + atom_z(iCount) * at_pos(iCount,2)
      uDipAt(3) = uDipAt(3) + atom_z(iCount) * at_pos(iCount,3)
   enddo

   if (size(mm_charges,1) > size(atom_z,1)) then
      do iCount = size(atom_z,1)+1, size(mm_charges,1)
         Qc = Qc + mm_charges(iCount)
      enddo
   endif

! Factor : For charged species dipole moment depends on the definition of the  !
! coordinates' origin. Using this factor, it is defined with respect to the    !
! center of charge (important in Reaction Field calculations). For neutral     !
! systems this is not necessary.                                               !

   factor = Qc / dble(nElec)
   uDip   = (uDipAt - uDip * factor) * 2.54D0

   if (present(print_dip)) then
      if (print_dip) call print_dipole(uDip)
   endif
   deallocate(dip_mat, uDipAt)
end subroutine dipole

!%% WRITE_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole moment to output, where dipxyz is the dipole moment vector.!
subroutine print_dipole(dipxyz)
   use fileio_data    , only: get_style
   use properties_data, only: UIDs
   implicit none
   LIODBLE, intent(in) :: dipxyz(3)
   
   integer :: UID
   LIODBLE :: u_abs
   logical :: style

   u_abs = dipxyz(1) * dipxyz(1) + dipxyz(2) * dipxyz(2) + &
           dipxyz(3) * dipxyz(3)
   u_abs = sqrt(u_abs)

   UID = UIDs%dip
   call get_style(style)
   if (style) then
      write(UID,8698)
      write(UID,8699)
      write(UID,8700)
      write(UID,8701)
      write(UID,8702)
      write(UID,8704) dipxyz(1), dipxyz(2), dipxyz(3), u_abs
      write(UID,8703)
   else
      write(UID,*)
      write(UID,'(A)') '#DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
      write(UID,*)
      write(UID,'(4(2x,F13.9))') dipxyz(1), dipxyz(2), dipxyz(3), u_abs
   endif

 8698 FORMAT(4x,"╔════════════════",&
      "═════════════════════",&
      "═════════════════════","═════╗")
 8699 FORMAT(4x,"║                         Dipole Moment            ", &
      "             ║")
 8700 FORMAT(4x,"╠═══════════════╦", &
      "═══════════════╦═════",       &
      "══════════╦══════════",       &
      "═════╣")
 8701 FORMAT(4x,"║       ux      ║       uy      ║       uz     ",&
      " ║       u       ║")
 8702 FORMAT(4x,"╠═══════════════╬", &
      "═══════════════╬═════",       &
      "══════════╬══════════",       &
      "═════╣")
 8704 FORMAT(4x,4("║",F13.9,2x),"║")
 8703 FORMAT(4x,"╚═══════════════╩", &
      "═══════════════╩═════",       &
      "══════════╩══════════",       &
      "═════╝")
end subroutine print_dipole

!%% WRITE_DIPOLE_TD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole momment vector in TD calculations, where dipxyz is the     !
! dipole moment vector, time is the current time (in fs) in TD.                !
subroutine write_dipole_td(dipxyz, time)
   use properties_data, only: UIDs
   implicit none
   LIODBLE, intent(in) :: dipxyz(3), time

   write(UIDs%diptd,100) time, dipxyz(1), dipxyz(2), dipxyz(3)

100 format (e15.8,1x, e15.8,1x,e15.8,1x,e15.8)
end subroutine write_dipole_td

!%% WRITE_DIPOLE_TD_HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole moment file's header in TD calculations, where time_step is!
! the one used in TD, fx-fy-fz are the perturbation field coordinates.         !
subroutine write_dipole_td_header(time_step, fx, fy, fz)
   use properties_data, only: UIDs
   implicit none
   LIODBLE, intent(in) :: time_step, fx, fy, fz

   write(UIDs%diptd, 100) time_step, fx, fy, fz

100 format ('# ',e12.5,1x, e12.5,1x,e12.5,1x,e12.5)
end subroutine write_dipole_td_header