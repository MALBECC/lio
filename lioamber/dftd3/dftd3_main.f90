! General setup
subroutine dftd3_setup()
   implicit none
end subroutine dftd3_setup

subroutine dftd3_finalise()
   implicit none
end subroutine dftd3_finalise

! Energy calculations
subroutine dftd3_energy(e_disp, dists, n_atoms)
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:)
   real(kind=8), intent(inout) :: e_disp
   
   real(kind=8) :: e_disp2, e_disp3

   e_disp2 = 0.0D0
   e_disp3 = 0.0D0
   
   call dftd3_2bodies_e(e_disp2, dists, n_atoms)
   call dftd3_3bodies_e(e_disp3, dists, n_atoms)

   ! The 3-body term should be negative (as C9 is negative),
   ! so E3 is added and not substracted.
   e_disp = e_disp - e_disp2 + e_disp3
end subroutine
