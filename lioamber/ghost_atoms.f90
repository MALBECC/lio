!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% GHOST_ATOMS.F90 %%%%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This files contains variables ghost_atoms and n_ghosts, as well as subroutine!
! summon_ghosts which sets the desired atomic nuclei (atom_Z) to 0.            !
! n_ghost is the number of ghost_atoms in the system, while ghost_atoms()      !
! contains the index of each ghost.                                            !
! To use this module, set n_ghosts > 1 and ghost_atoms=x,y,z,..,n where x,y,z  !
! and n indicate the index for each ghost atom (in the order that are placed in!
! the coordinates file).                                                       !
!                                                                              !
! WARNING: The following module assumes the maximum amount of QM atoms is 300. !
!          While this is not optimal, a similar cap is set in module G2G.      !
! F. Pedron - July 2018                                                        !
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
   subroutine adjust_ghost_charge(atom_Z, n_atoms, nuc_charge)
      use ghost_atoms_data, only: ghost_atoms, n_ghosts
      implicit none
      integer, intent(in)    :: n_atoms, atom_Z(n_atoms)
      integer, intent(inout) :: nuc_charge
      integer :: icount

      if (n_ghosts .lt. 1) return;

      do icount = 1, n_ghosts
         nuc_charge = nuc_charge - atom_Z(ghost_atoms(icount))
      enddo

      return
   end subroutine adjust_ghost_charge

   subroutine summon_ghosts(atom_Z, n_atoms, verbose)
      use ghost_atoms_data, only: ghost_atoms, n_ghosts

      implicit none
      integer, intent(in)    :: n_atoms, verbose
      integer, intent(inout) :: atom_Z(n_atoms)
      integer :: icount

      if (n_ghosts .lt. 1) return;
      if (verbose.gt.3) then
         write(*,'(A)', advance = 'no') "  Ghost atom indexes = "
         write(*,*) ghost_atoms
      endif

      do icount=1, n_ghosts
         if (ghost_atoms(icount) .lt. 1) then
            write(*,'(A)') " ERROR - SUMMON_GHOSTS: Ghost_atoms(x) cannot be"&
                 &" less than 1."
         endif
         atom_Z(ghost_atoms(icount)) = 0
      enddo

      return
   end subroutine summon_ghosts
end module ghost_atoms_subs
