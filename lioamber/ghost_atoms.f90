!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% GHOST_ATOMS.F90 %%%%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This files contains variables ghost_atoms and n_ghosts, as well as subroutine!
! summon_ghosts which sets the desired atomic nuclei (atom_Z) to 0.            !
! n_ghost is the number of ghost_atoms in the system, while ghost_atoms()      !
! contains the index of each ghost.                                            !
! WARNING: The following module assumes the maximum amount of QM atoms is 300. !
!          While this is not optimal, a similar cap is set in module G2G.      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ghost_atoms_data
   implicit none
   integer      :: ghost_atoms(300)
   integer      :: n_ghosts = 0
contains
end module ghost_atoms_data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ghost_atoms_subs
contains
   subroutine summon_ghosts(Iz, natom)
      use ghost_atoms_data, only: ghost_atoms, n_ghosts

      implicit none
      integer, intent(in)    :: natom
      integer, intent(inout) :: Iz(natom)
      integer :: icount

      if (n_ghosts .lt. 1) return;
      do icount=1, n_ghosts
         if (ghost_atoms(icount) .lt. 1) then
            write(*,'(A)') " ERROR - SUMMON_GHOSTS: Ghost_atoms(x) cannot be"&
                 &" less than 1."
         endif
         Iz(ghost_atoms(icount)) = 0
      enddo

      return
   end subroutine summon_ghosts
end module ghost_atoms_subs
