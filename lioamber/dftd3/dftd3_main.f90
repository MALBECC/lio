! General setup
subroutine dftd3_setup(n_atoms, atom_z)
   use dftd3_data, only: dftd3, c6_ab, r0_ab, c8_ab, c6_cn, r_cov, c8_coef
   implicit none
   integer, intent(in)  :: n_atoms, atom_z(:)

   if (.not. dftd3) return
   if (.not. allocated(c6_ab  )) allocate(c6_ab(n_atoms,n_atoms))
   if (.not. allocated(r0_ab  )) allocate(r0_ab(n_atoms,n_atoms))
   if (.not. allocated(c8_ab  )) allocate(c8_ab(n_atoms,n_atoms))
   if (.not. allocated(c6_cn  )) allocate(c6_cn(n_atoms,n_atoms,5,5,3))
   if (.not. allocated(r_cov  )) allocate(r_cov(n_atoms))
   if (.not. allocated(c8_coef)) allocate(c8_coef(n_atoms))


   call dftd3_read_c6(c6_cn  , n_atoms, atom_z)
   call dftd3_read_r0(r0_ab  , n_atoms, atom_z)
   call dftd3_read_rc(r_cov  , n_atoms, atom_z)
   call dftd3_read_c8(c8_coef, n_atoms, atom_z)
end subroutine dftd3_setup

subroutine dftd3_finalise()
   use dftd3_data, only: dftd3, c6_ab, r0_ab, c8_ab, c6_cn, r_cov, c8_coef
   implicit none

   if (.not. dftd3) return
   if (allocated(c6_ab)) deallocate(c6_ab)
   if (allocated(r0_ab)) deallocate(r0_ab)
   if (allocated(c8_ab)) deallocate(c8_ab)

   if (allocated(c6_cn))   deallocate(c6_cn)
   if (allocated(r_cov))   deallocate(r_cov)
   if (allocated(c8_coef)) deallocate(c8_coef)
end subroutine dftd3_finalise

! Energy calculations
subroutine dftd3_energy(e_disp, dists, n_atoms)
   use dftd3_data, only: dftd3
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:)
   real(kind=8), intent(inout) :: e_disp
   
   real(kind=8) :: e_disp2, e_disp3

   if (.not. dftd3) return
   e_disp2 = 0.0D0
   e_disp3 = 0.0D0
   
   call dftd3_set_c6c8(dists, n_atoms)
   call dftd3_2bodies_e(e_disp2, dists, n_atoms)
   call dftd3_3bodies_e(e_disp3, dists, n_atoms)

   ! The 3-body term should be negative (as C9 is negative),
   ! so E3 is added and not substracted.
   e_disp = e_disp - e_disp2 + e_disp3
end subroutine dftd3_energy

! Gradient calculations. This needs to run the energy routine
! for the same set of positions, in order to keep the C6 and
! C8 coefficients appropriate.
subroutine dftd3_gradients(grad, dists, pos, n_atoms)
   use dftd3_data, only: dftd3
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:), pos(:,:)
   real(kind=8), intent(inout) :: grad(:,:)
   
   if (.not. dftd3) return
   call dftd3_2bodies_g(grad, dists, pos, n_atoms)
   call dftd3_3bodies_g(grad, dists, pos, n_atoms)
end subroutine dftd3_gradients
