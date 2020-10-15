!%% WRITE_POPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes Mulliken/Löwdin charges to output. q0 is the array of atomic numbers, !
! q is their charge, pop is the type of population, analysis performed         !
! (0 = Mulliken, 1 = Löwdin, 2 = Becke) and UID is the output file unit ID.    !
subroutine write_population(q0, q, pop, UID, filename)
   use fileio_data, only : style
   implicit none
   integer         , intent(in) :: UID
   integer         , intent(in) :: pop
   integer         , intent(in) :: q0(:)
   LIODBLE         , intent(in) :: q(:)
   character(len=*), intent(in) :: filename

   LIODBLE :: qtotal
   integer :: icount
   logical :: file_open

   ! Checks if file is open.
   inquire(unit = UID, opened = file_open)
   if (.not. file_open) open(unit = UID, file = filename)

   qtotal = 0.0D0
   write(UID,*)
   if (style) then
      write(UID,300)
      if (pop == 0) write(UID,301)
      if (pop == 1) write(UID,309)
      if (pop == 2) write(UID,310)
      if (pop == 3) write(UID,311)
      if (pop == 4) write(UID,312)
      if (pop == 5) write(UID,313)
      write(UID,302); write(UID,303); write(UID,304)
      do icount = 1, N
         qtotal = qtotal + q(icount)
         write(UID,305) icount, q0(icount), q(icount)
      enddo
      write(UID,306)
      write(UID,307) qtotal
      write(UID,308)
   else
      if (pop == 0) write(UID,402) "# Mulliken Population Analysis"
      if (pop == 1) write(UID,402) "# Mulliken Spin Population Analysis"
      if (pop == 2) write(UID,402) "# Löwdin Population Analysis"
      if (pop == 3) write(UID,402) "# Löwdin Spin Population Analysis"
      if (pop == 4) write(UID,402) "# Becke Population Analysis"
      if (pop == 5) write(UID,402) "# Becke Spin Population Analysis"
      write(UID,402) "# Atom   Type   Population"
      do icount = 1, size(q0,1)
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
305 FORMAT(8x,"║",2x,i3,3x,"║",3x,i3,5x,"║",1x,F10.7,1x,"║")
306 FORMAT(8x,"╚════════╬═══════════╬════════════╣")
307 FORMAT(8x,"         ║   TOTAL   ║",1x,F10.7,1x,"║")
308 FORMAT(8x,"         ╚═══════════╩════════════╝")
309 FORMAT(8x,"║MULLIKEN SPIN POPULATION ANALYSIS║")
310 FORMAT(8x,"║    LÖWDIN POPULATION ANALYSIS   ║")
311 FORMAT(8x,"║ LÖWDIN SPIN POPULATION ANALYSIS ║")
312 FORMAT(8x,"║     BECKE POPULATION ANALYSIS   ║")
313 FORMAT(8x,"║  BECKE SPIN POPULATION ANALYSIS ║")
400 FORMAT(2x,i3,4x,i3,5x,F10.7)
401 FORMAT(2x,"Total Charge = ", F10.7)
402 FORMAT(A)
end subroutine write_population

!%% WRITE_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes Fukui function and local softness to output, where fukuiZZZ are the   !
! atomic Fukui functions for the nuclei (positive, negative or radical), N is  !
! the total amount of atoms, Iz their atomic number, and soft is the global    !
! softness for the molecule.                                                   !
subroutine write_fukui_core(fukuiNeg, fukuiPos, fukuiRad, Iz, soft)
   implicit none
   integer, intent(in) :: atom_z(:)
   LIODBLE, intent(in) :: fukuiNeg(:)
   LIODBLE, intent(in) :: fukuiPos(:)
   LIODBLE, intent(in) :: fukuiRad(:)
   LIODBLE, intent(in) :: soft
   
   integer :: icount
   logical :: is_open

   inquire(unit = 1984, opened = is_open)
   if (.not. is_open) open(file = 'fukui', unit = 1984)

   write(1984,'(A)') "Condensed to Atoms Fukui Function"
   write(1984,'(A26,F14.7)') "Global Softness (A.U.): ", soft
   write(1984,'(A)') "  N     Fukui-       Fukui+       Fukui0  &
                     &  Local Softness (A.U.)"
   do icount = 1, size(atom_z,1)
      write(1984,'(I3,2x,F12.9,2x,F12.9,2x,F12.9,2x,F14.7)') atom_z(icount),  &
                        fukuiNeg(icount), fukuiPos(icount), fukuiRad(icount), &
                        abs(soft*fukuiRad(icount))
   enddo

end subroutine write_fukui_core
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!