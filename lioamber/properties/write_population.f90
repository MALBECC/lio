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
      do icount = 1, size(q0,1)
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
   call write_population_per_region(q, pop, UID+50, filename)

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

subroutine write_population_per_region(q, pop, UID, filename)
   use properties_data, only : prop_regions
   implicit none
   integer         , intent(in) :: UID
   integer         , intent(in) :: pop
   LIODBLE         , intent(in) :: q(:)
   character(len=*), intent(in) :: filename

   LIODBLE :: qtotal
   integer :: ireg, iatom
   logical :: file_open
   character(len=60) :: new_fname

   if (prop_regions%n_regions == 0) return
   new_fname = "regions_"//trim(filename)

   ! Checks if file is open.
   inquire(unit = UID, opened = file_open)
   if (.not. file_open) open(unit = UID, file = new_fname)

   write(UID,*)
   if (pop == 0) write(UID,402) "# Mulliken Population Analysis"
   if (pop == 1) write(UID,402) "# Mulliken Spin Population Analysis"
   if (pop == 2) write(UID,402) "# Löwdin Population Analysis"
   if (pop == 3) write(UID,402) "# Löwdin Spin Population Analysis"
   if (pop == 4) write(UID,402) "# Becke Population Analysis"
   if (pop == 5) write(UID,402) "# Becke Spin Population Analysis"
   write(UID,402) "# Region  N° atoms  Population"
   do ireg = 1, prop_regions%n_regions

      qtotal = 0.0D0
      do iatom = 1, prop_regions%natoms(ireg)
         qtotal = qtotal + q(prop_regions%atoms(ireg,iatom))
      enddo
      write(UID,400) ireg, prop_regions%natoms(ireg), qtotal
   enddo
   write(UID,*)

400 FORMAT(2x,i3,4x,i3,5x,F10.7)
402 FORMAT(A)
end subroutine write_population_per_region