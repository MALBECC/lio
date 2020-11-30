!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DIPOLE.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains all dipole-calculation and dipole-printing related        ! 
! subroutines.                                                                 !
! Called after a single-point or similar calculation:                          !
!  * dipole():              calculates the XYZ dipole moment and prints it if  !
!                           needed. This includes a "per-region" dipole moment !
!                           calculation.                                       !
!  * print_dipole():        Prints dipole moment to output. Called by dipole().!
!  * region_print_dipole(): Same as above, but as per-region dipole.           !
! Called during TD, the last two called by dipole():                                                            !
!  * write_dipole_td_header(): Prints TD dipole file header.                   !
!  * write_dipole_td():        Prints dipole for TD calculations.              !
!  * region_write_dipole_td(): Prints per-region dipole for TD calculations.   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

! Calculates dipole moment from dipole matrix.
! Inputs:
!   Pmat_v:     Atomic-basis density matrix in vector form.
!   nElec:      Number of electrons.
!   at_pos:     Array containing all atomic positions (size ntatom, 3).
!   at_dists:   Array containing distances between atoms (size natom, natom).
!   atom_z:     Array containing atomic numbers for QM region. (size natom).
!   mm_charges: Array containg MM partial charges for QM/MM (size ntatom).
!   print_dip:  (Optional) If present and > 0, prints dipole to output.
!   td_time:    (Optional) TD time, used if print_dip=2 (i.e. print TD dipole).
! Output:
!   uDip:       Total dipole moment vector.
subroutine dipole(uDip, Pmat_v, nElec, at_pos, at_dists, atom_z, mm_charges, &
                  dipNUC, print_dip, td_time)
   use faint_cpu      , only: intdip
   use properties_data, only: prop_regions

   implicit none
   integer, intent(in)           :: nElec
   integer, intent(in)           :: atom_z(:)
   LIODBLE, intent(in)           :: Pmat_v(:)
   LIODBLE, intent(in)           :: at_pos(:,:)
   LIODBLE, intent(in)           :: at_dists(:,:)
   LIODBLE, intent(in)           :: mm_charges(:)
   integer, intent(in), optional :: print_dip
   LIODBLE, intent(in), optional :: td_time
   LIODBLE, intent(inout)        :: uDip(3)
   logical, intent(in)           :: dipNUC

   LIODBLE, allocatable :: dip_mat(:,:), uDipAt(:), uDip_reg(:,:), &
                           dip_mat_unpacked(:,:)
   LIODBLE :: Qc, factor
   integer :: iCount, iFunct, jFunct, iReg, iatom, tIndex, Msize
   
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
   if ( dipNUC ) then
      uDip   = (uDipAt - uDip * factor) * 2.54D0
   endif

   if (present(print_dip)) then
      if (print_dip == 1) call print_dipole(uDip)
      if (print_dip == 2) call write_dipole_td(uDip, td_time)
   endif

   ! Now separates dipole contributions by region.
   if (prop_regions%n_regions < 1) then
      deallocate(dip_mat, uDipAt)
      return
   endif

   Msize = int(sqrt(2.0D0 * size(Pmat_v,1)))
   allocate(uDip_reg(3,prop_regions%n_regions))

   allocate(dip_mat_unpacked(Msize,Msize))
   uDip_reg = 0.0D0

   ! We use a Mulliken-like projection for the dipole matrix,
   ! but we accumulate the atoms into the regions themselves.
   do iReg = 1, prop_regions%n_regions
      uDip = 0.0D0

      ! We work over each component of the dipole moment, X-Y-Z
      do iCount = 1, 3
         call spunpack_rho('L', Msize, dip_mat(iCount,:), dip_mat_unpacked)

         do iFunct = 1, prop_regions%nfuncs(iReg)
            tIndex = prop_regions%funcs(iReg,iFunct)
            
            do jFunct = 1, Msize
               uDip(iCount) = uDip(iCount) + dip_mat_unpacked(tIndex,jFunct)
            enddo
         enddo
      enddo

      uDipAt = 0.0D0
      do iCount = 1, prop_regions%natoms(iReg)
         iatom = prop_regions%atoms(iReg,iCount)
         uDipAt(1) = uDipAt(1) + atom_z(iatom) * at_pos(iatom,1)
         uDipAt(2) = uDipAt(2) + atom_z(iatom) * at_pos(iatom,2)
         uDipAt(3) = uDipAt(3) + atom_z(iatom) * at_pos(iatom,3)
      enddo

      uDip_reg(:,iReg) = (uDipAt - uDip * factor) * 2.54D0
   enddo

   if (present(print_dip)) then
      if (print_dip == 1) call region_print_dipole(uDip_reg)
      if (print_dip == 2) call region_write_dipole_td(uDip_reg, td_time)
   endif
  
   deallocate(dip_mat, uDipAt, uDip_reg, dip_mat_unpacked)
end subroutine dipole

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

! Prints the dipole momment vector in TD calculations, where dipxyz is the     !
! dipole moment vector, time is the current time (in fs) in TD.                !
subroutine write_dipole_td(dipxyz, time)
   use properties_data, only: UIDs
   implicit none
   LIODBLE, intent(in) :: dipxyz(3), time

   write(UIDs%diptd,100) time, dipxyz(1), dipxyz(2), dipxyz(3)

100 format (e15.8,1x, e15.8,1x,e15.8,1x,e15.8)
end subroutine write_dipole_td

! Prints the dipole moment file's header in TD calculations, where time_step is!
! the one used in TD, fx-fy-fz are the perturbation field coordinates.         !
subroutine write_dipole_td_header(time_step, fx, fy, fz)
   use properties_data, only: UIDs
   implicit none
   LIODBLE, intent(in) :: time_step, fx, fy, fz

   write(UIDs%diptd, 100) time_step, fx, fy, fz

100 format ('# ',e12.5,1x, e12.5,1x,e12.5,1x,e12.5)
end subroutine write_dipole_td_header

! Prints the dipole moment vector per region.                                  !
subroutine region_print_dipole(dipxyz)
   use properties_data, only: UIDs, prop_regions
   implicit none
   LIODBLE, intent(in) :: dipxyz(:,:)
   
   integer :: UID, iReg
   LIODBLE :: u_abs

   UID = UIDs%dip+50
   write(UID,*)
   write(UID,'(A6,5x,A5,10x,A5,10x,A5,10x,A13)') 'REGION', 'DIP_X', 'DIP_Y', &
                                              'DIP_Z', 'NORM (DEBYES)'
   do iReg = 1, prop_regions%n_regions
      u_abs = dipxyz(1,iReg) * dipxyz(1,iReg) + &
              dipxyz(2,iReg) * dipxyz(2,iReg) + &
              dipxyz(3,iReg) * dipxyz(3,iReg)
      u_abs = sqrt(u_abs)
            
      write(UID,'(1x,I3,4(2x,F13.9))') iReg, dipxyz(1,iReg), dipxyz(2,iReg), &
                                    dipxyz(3,iReg), u_abs
   enddo

end subroutine region_print_dipole

! Prints the dipole momment vector in TD calculations, separated by regions.   !
subroutine region_write_dipole_td(dipxyz, time)
   use properties_data, only: UIDs, prop_regions

   implicit none
   LIODBLE, intent(in) :: dipxyz(:,:), time
   integer :: iReg

   do iReg = 1, prop_regions%n_regions
      write(UIDs%diptd+50,100) time, iReg, dipxyz(1,iReg), dipxyz(2,iReg), &
                               dipxyz(3,iReg)
   enddo

100 format (e15.8,1x,I3,1x,e15.8,1x,e15.8,1x,e15.8)
end subroutine region_write_dipole_td