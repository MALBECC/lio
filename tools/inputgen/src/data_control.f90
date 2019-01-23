!##############################################################################!
module data_control
!##############################################################################!
   implicit none

   integer :: nof_atoms
   integer, allocatable :: numid_of_atom(:)
   real*8 , allocatable :: coord_of_atom(:,:)

   
   integer :: nof_basis
   integer, allocatable :: atom_of_basis(:)
   integer, allocatable :: angm_of_basis(:)
   real*8 , allocatable :: densmat_valij(:,:)

contains
!
!
!
!
!##############################################################################!

subroutine reads_data()
   use keywords       , only: output_type, inpname_fchk

   use parser_gaussian, only: gaufchk_atoms, gaufchk_atoms_info &
                           &, gaufchk_basis, gaufchk_basis_info &
                           &, gaufchk_densmat

   implicit none
   integer :: density_state


!  SETUP PROCEDURES
   select case( output_type )
      case ( "full_ground" )
         density_state = 0

      case ( "full_excited" )
         density_state = 1

      case default
         print*, "ERROR: reads_data"
         print*, "Unidentified kind of request..."
         print*, "output_type = ", output_type
         stop

   end select

   call gaufchk_atoms( inpname_fchk, nof_atoms )
   allocate( numid_of_atom(nof_atoms) )
   allocate( coord_of_atom(3, nof_atoms) )
   call gaufchk_atoms_info( inpname_fchk, nof_atoms, numid_of_atom, coord_of_atom )

   call gaufchk_basis( inpname_fchk, nof_basis )
   allocate( atom_of_basis(nof_basis) )
   allocate( angm_of_basis(nof_basis) )
   allocate( densmat_valij(nof_basis,nof_basis) )
   call gaufchk_basis_info( inpname_fchk, nof_basis, atom_of_basis, angm_of_basis )
   call gaufchk_densmat( inpname_fchk, nof_basis, density_state, densmat_valij )


end subroutine reads_data

!------------------------------------------------------------------------------!

subroutine shape_data()
   use auxmod_subs, only: get_ordered_vec, reorder_int1vec, reorder_dbl3vec, &
                        & reorder_dblmat, update_idref

   implicit none
   integer              :: iatom
   integer              :: atomic_number
   integer              :: new_index
   integer, allocatable :: ordered_atoms(:)
   integer, allocatable :: ordered_basis(:)

   allocate( ordered_atoms(nof_atoms) )
   allocate( ordered_basis(nof_basis) )

   call get_ordered_vec( nof_atoms, numid_of_atom, ordered_atoms )
   call reorder_int1vec( nof_atoms, ordered_atoms, numid_of_atom )
   call reorder_dbl3vec( nof_atoms, ordered_atoms, coord_of_atom )

   call update_idref( nof_atoms, nof_basis, ordered_atoms, atom_of_basis )

   call get_ordered_vec( nof_basis, atom_of_basis, ordered_basis )
   call reorder_int1vec( nof_basis, ordered_basis, atom_of_basis )
   call reorder_int1vec( nof_basis, ordered_basis, angm_of_basis )
   call reorder_dblmat(  nof_basis, ordered_basis, densmat_valij )

   call get_ordered_vec( nof_basis, angm_of_basis, ordered_basis )
   call reorder_int1vec( nof_basis, ordered_basis, atom_of_basis )
   call reorder_int1vec( nof_basis, ordered_basis, angm_of_basis )
   call reorder_dblmat(  nof_basis, ordered_basis, densmat_valij )


end subroutine shape_data

!------------------------------------------------------------------------------!

subroutine write_data()
   use keywords   , only: output_xyz, output_rst
   use auxmod_subs, only: safe_open, safe_close
   use parser_lio , only: writelio_zyx, writelio_rst
   implicit none

   call safe_open( 201, output_xyz )
   call writelio_zyx( 201, nof_atoms, numid_of_atom, coord_of_atom )
   call safe_close( 201 )

   call safe_open( 201, output_rst )
   call writelio_rst( 201, nof_basis, densmat_valij )
   call safe_close( 201 )


end subroutine write_data

!##############################################################################!
end module data_control
!##############################################################################!
