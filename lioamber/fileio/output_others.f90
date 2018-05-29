!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_OUTPUT.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains several output-file printing routines. Currently includes:!
! * write_dipole     (handles dipole moment printing)                          !
! * write_dipole_td  (handles dipole moment printing in TD)                    !
! * write_forces     (handles grandient printing to output)                    !
! * write_fukui      (handles Fukui function printing to output)               !
! * write_orbitals   (prints orbitals and energies to output)                  !
! * write_orbitals_op(prints orbitals and energies to output, open shell)      !
! * write_population (handles population/charge printing to output)            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole moment to output, where dipxyz is the dipole moment vector,!
! u is its norm, uid is the output UID, and header decides whether to print a  !
! header or not.                                                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole(dipxyz, u, uid, header)
   use fileio_data, only : style
   implicit none
   double precision, intent(in) :: dipxyz(3), u
   integer         , intent(in) :: uid
   logical         , intent(in) :: header

   open(unit = uid, file = "dipole_moment")
   if (style) then
      if (header) then
         write(UID,8698)
         write(UID,8699)
         write(UID,8700)
         write(UID,8701)
         write(UID,8702)
      else
         write(UID,8704) dipxyz(1), dipxyz(2), dipxyz(3), u
      endif
   else
      if (header) then
         write(UID,*)
         write(UID,*) '#DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
         write(UID,*)
      else
         write(UID,*) dipxyz(1), dipxyz(2), dipxyz(3), u
      endif
   endif

   return
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
 8704 FORMAT(4x,4("║"F13.9,2x),"║")
end subroutine write_dipole
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_DIPOLE_TD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole momment vector in TD calculations, where dipxyz is the     !
! dipole moment vector, time is the current time (in fs) in TD, and uid is the !
! output file UID.                                                             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole_td(dipxyz, time, uid)
   implicit none
   double precision, intent(in) :: dipxyz(3), time
   integer         , intent(in) :: uid

   write(UID,100) time, dipxyz(1), dipxyz(2), dipxyz(3)

   return
100 format (e15.8,' ', e15.8,' ',e15.8,' ',e15.8)
end subroutine write_dipole_td
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_DIPOLE_TD_HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole moment file's header in TD calculations, where time_step is!
! the one used in TD, fx-fy-fz are the perturbation field coordinates, and uid !
! is the output UID.                                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole_td_header(time_step, fx, fy, fz, uid)
   implicit none
   double precision, intent(in) :: time_step, fx, fy, fz
   integer         , intent(in) :: uid

   write(UID, 100) time_step, fx, fy, fz

   return
100 format ('# ',e12.5,' ', e12.5,' ',e12.5,' ',e12.5)
end subroutine write_dipole_td_header
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the final forces for chosen atoms to output. dxyz is the gradient     !
! vector, natom is the total amount of atoms to print, offset is the starting  !
! atom number, and uid is the output file UID.                                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_forces(dxyz, natom, offset, uid)
   implicit none
   integer         , intent(in) :: uid, natom, offset
   double precision, intent(in) :: dxyz(3, natom+offset)
   integer :: icount

   do icount = offset+1, offset+natom
      write(uid,100) icount, dxyz(1, icount), dxyz(2, icount), dxyz(3, icount)
   enddo

   return
100 format (I5,2x,f10.6,2x,f10.6,2x,f10.6)
end subroutine write_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes Fukui function and local softness to output, where fukuiZZZ are the   !
! atomic Fukui functions for the nuclei (positive, negative or radical), N is  !
! the total amount of atoms, Iz their atomic number, and soft is the global    !
! softness for the molecule.                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_fukui(fukuiNeg, fukuiPos, fukuiRad, N, Iz, soft)
   implicit none
   integer         , intent(in) :: N, Iz(N)
   double precision, intent(in) :: fukuiNeg(N), fukuiPos(N), fukuiRad(N), soft
   integer :: i

   write(*,*) "Global Softness (A.U.):  ", soft
   write(*,*) "N", "Fukui-", "Fukui+", "Fukui0", "Local Softness (A.U.)"
   do i=1, N
      write(*,*) Iz(i), fukuiNeg(i), fukuiPos(i), fukuiRad(i), &
                 abs(soft*fukuiRad(i))
   enddo

   return
end subroutine write_fukui
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_ORBITALS / WRITE_ORBITALS_OPEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints orbital energies and coefficients. M is the total number of basis,    !
! NCO is the number of occupied orbitals, Eorbs is the energy of each MO,      !
! MO_coeff is the MO coefficient matrix, and uid is the output file UID.       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_orbitals(M, NCO, E_orbs, MO_coeff, uid)
   implicit none
   integer, intent(in)          :: M, NCO, uid
   double precision, intent(in) :: E_orbs(M), MO_coeff(M*M)
   integer                      :: icount, jcount

   write(uid,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
   do icount = 1, NCO
      write(UID, 850) icount, E_orbs(icount)
      do jcount = 1, M
         write(UID, 400) MO_coeff(jcount + icount)
      enddo
   enddo
   do icount = NCO+1, M
      write(uid, 851) icount, E_orbs(icount)
      do jcount = 1, M
         write(uid, 400) MO_coeff(jcount + icount)
      enddo
   enddo

   return
400 format(4(E14.7E2, 2x))
850 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
851 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7, '(NON OCC.)')
end subroutine write_orbitals

subroutine write_orbitals_op(M, NCO, NUnp, E_orbs, E_orbs_b, MO_coeff, &
                               MO_coeff_b, uid)
   implicit none
   integer, intent(in)          :: M, NCO, NUnp, uid
   double precision, intent(in) :: E_orbs(M), E_orbs_b(M), MO_coeff(M*M), &
                                   MO_coeff_b(M*M)
   integer                      :: icount, jcount

   write(uid,*) 'ORBITAL COEFFICIENTS AND ENERGIES, OPEN SHELL ALPHA'
   do icount = 1, NCO
      write(UID, 850) icount, E_orbs(icount)
      do jcount = 1, M
         write(UID, 400) MO_coeff(jcount + icount)
      enddo
   enddo
   do icount = NCO+1, M
      write(uid, 851) icount, E_orbs(icount)
      do jcount = 1, M
         write(uid, 400) MO_coeff(jcount + icount)
      enddo
   enddo

   write(uid,*) 'ORBITAL COEFFICIENTS AND ENERGIES, OPEN SHELL BETA'
   do icount = 1, NCO+NUnp
      write(UID, 850) icount, E_orbs_b(icount)
      do jcount = 1, M
         write(UID, 400) MO_coeff_b(jcount + icount)
      enddo
   enddo
   do icount = NCO+NUnp+1, M
      write(uid, 851) icount, E_orbs(icount)
      do jcount = 1, M
         write(uid, 400) MO_coeff_b(jcount + icount)
      enddo
   enddo

   return
400 format(4(E14.7E2, 2x))
850 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
851 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7, '(NON OCC.)')
end subroutine write_orbitals_op
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_POPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes Mulliken/Löwdin charges to output. N is the number of atoms, q0 is    !
! their atomic number, q is their charge, pop is the type of population        !
! analysis performed (0=Mulliken, 1=Löwdin) and UID is the output file UID.    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_population(N, q0, q, pop, UID)
   use fileio_data, only : style
   implicit none
   integer         , intent(in) :: UID, N, q0(N), pop
   double precision, intent(in) :: q(N)
   double precision :: qtotal
   integer          :: icount

   qtotal = 0.0D0
   write(UID,*)
   if (style) then
      write(UID,300)
      if (pop.eq.0) write(UID,301)
      if (pop.eq.1) write(UID,309)
      write(UID,302); write(UID,303); write(UID,304)
      do icount = 1, N
         qtotal = qtotal + q(icount)
         write(UID,305) icount, q0(icount), q(icount)
      enddo
      write(UID,306)
      write(UID,307) qtotal
      write(UID,308)
   else
      if (pop.eq.0) write(UID,402) "# Mulliken Population Analysis"
      if (pop.eq.1) write(UID,402) "# Löwdin Population Analysis"
      write(UID,402) "# Atom   Type   Population"
      do icount = 1, N
         qtotal = qtotal + q(icount)
         write(UID,400) icount, q0(icount), q(icount)
      enddo
      write(UID,401) qtotal
   endif
   write(UID,*)

   return
300 FORMAT(8x,"╔════════════════", &
    "═════════════════╗")
301 FORMAT(8x,"║   MULLIKEN POPULATION ANALYSIS  ║")
302 FORMAT(8x,"╠════════╦═══════════╦════════════╣")
303 FORMAT(8x,"║ ATOM # ║ ATOM TYPE ║ POPULATION ║")
304 FORMAT(8x,"╠════════╬═══════════╬════════════╣")
305 FORMAT(8x,"║",2x,i3,3x,"║"3x,i3,5x,"║",1x,F10.7,1x,"║")
306 FORMAT(8x,"╚════════╬═══════════╬════════════╣")
307 FORMAT(8x,"         ║   TOTAL   ║",1x,F10.7,1x,"║")
308 FORMAT(8x,"         ╚═══════════╩════════════╝")
309 FORMAT(8x,"║    LÖWDIN POPULATION ANALYSIS   ║")
400 FORMAT(2x,i3,4x,i3,5x,F10.7)
401 FORMAT(2x,"Total Charge = ", F10.7)
402 FORMAT(A)
end subroutine write_population
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
