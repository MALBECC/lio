! 3-body terms for DFTD3 corrections to energy and gradients.
subroutine dftd3_2bodies_e(e_2, dists, n_atoms)
   use dftd3_data, only: c6_ab, c8_ab, r0_ab
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:)
   real(kind=8), intent(inout) :: e_2

   real(kind=8) :: term
   integer      :: iatom, jatom

   do iatom = 1      , n_atoms
   do jatom = iatom+1, n_atoms
      ! r^6 terms.
      term = 1.0D0 + 6.0D0 * ( (1.217D0 * r0_ab(iatom,jatom)) / &
                                dists(iatom,jatom) ) ** 14
      e_2 = e_2 + c6_ab(iatom,jatom) / (term * dists(iatom,jatom) ** 6)

      ! r^8 term.
      term = 1.0D0 + 6.0D0 * (r0_ab(iatom,jatom) / dists(iatom,jatom)) ** 16
      e_2 = e_2 + c8_ab(iatom,jatom) / (term * dists(iatom,jatom) ** 8)
   enddo
   enddo   
end subroutine dftd3_2bodies_e

subroutine dftd3_2bodies_g(e_2, dists, pos, n_atoms)
   use dftd3_data, only: c6_ab, c8_ab, r0_ab
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:), pos(:)
   real(kind=8), intent(inout) :: e_2

   real(kind=8) :: term
   integer      :: iatom, jatom

   do iatom = 1      , n_atoms
   do jatom = iatom+1, n_atoms
      ! r^6 terms.
      term = 1.0D0 + 6.0D0 * ( (1.217D0 * r0_ab(iatom,jatom)) / &
                                dists(iatom,jatom) ) ** 14
      e_2 = e_2 + c6_ab(iatom,jatom) / (term * dists(iatom,jatom) ** 6)

      ! r^8 term.
      term = 1.0D0 + 6.0D0 * (r0_ab(iatom,jatom) / dists(iatom,jatom)) ** 16
      e_2 = e_2 + c8_ab(iatom,jatom) / (term * dists(iatom,jatom) ** 8)
   enddo
   enddo   
end subroutine dftd3_2bodies_g
