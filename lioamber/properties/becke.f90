! Prints Becke population analysis, intended
! to be called from liomain. Most of the
! calculations are done on g2g, so everything
! should already be initialised there.
! Inputs:
!   * atom_z:     Array containing atomic numbers.
!   * open_shell: Logical indicating open shell calculation.
subroutine print_becke(atom_z, open_shell)
   use properties_data, only: UIDs, fbecke
   implicit none
   integer, intent(in)  :: atom_z(:)
   logical, intent(in)  :: open_shell

   LIODBLE, allocatable :: q(:)
   character(len=100)   :: spinfile

   call g2g_timer_sum_start("Becke")
   allocate(q(size(atom_z,1)))
   
   call g2g_get_becke_dens(q)
   call write_population(atom_z, q, 4, UIDs%bec, fbecke)

   if (open_shell) then
      spinfile = trim(fbecke) // "_spin"
      call g2g_get_becke_spin(q)
      call write_population(atom_z, q, 5, UIDs%becs, spinfile)
   endif

   deallocate(q)
   call g2g_timer_sum_pause("Becke")
end subroutine print_becke