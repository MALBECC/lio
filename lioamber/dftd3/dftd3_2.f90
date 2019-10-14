! 3-body terms for DFTD3 corrections to energy and gradients.
subroutine dftd3_2bodies_e(e_2, dists, n_atoms)
   use dftd3_data, only: c6_ab, c8_ab, r0_ab, dftd3_s6, dftd3_s8, dftd3_sr6
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:)
   real(kind=8), intent(inout) :: e_2

   real(kind=8) :: term, e_c6, e_c8
   integer      :: iatom, jatom

   e_c6 = 0.0D0
   e_c8 = 0.0D0
   do iatom = 1      , n_atoms
   do jatom = iatom+1, n_atoms
      ! r^6 terms.
      term = dftd3_sr6 * r0_ab(iatom,jatom) / dists(iatom,jatom)
      term = 1.0D0 / (1.0D0 + 6.0D0 * term ** 14)
      e_c6 = e_c6 + &
             dftd3_s6 * term * c6_ab(iatom,jatom) / dists(iatom,jatom) ** 6

      ! r^8 term.
      term = 1.0D0 + 6.0D0 * (r0_ab(iatom,jatom) / dists(iatom,jatom)) ** 16
      e_c8 = e_c8 + dftd3_s8 * c8_ab(iatom,jatom) / &
                             (term * dists(iatom,jatom) ** 8)
   enddo
   enddo

   e_2 = e_c6 + e_c8
end subroutine dftd3_2bodies_e

subroutine dftd3_2bodies_g(grad, dists, pos, n_atoms)
   use dftd3_data, only: c6_ab, c8_ab, r0_ab, dftd3_s6, dftd3_s8, dftd3_sr6
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:), pos(:,:)
   real(kind=8), intent(inout) :: grad(:,:)

   real(kind=8) :: rtemp, btemp, grad6, grad8
   integer      :: iatom, jatom

   do iatom = 1      , n_atoms
   do jatom = iatom+1, n_atoms
      ! r^6 terms.
      rtemp = dists(iatom,jatom) ** 7
      btemp = dftd3_sr6 * r0_ab(iatom,jatom)
      btemp = btemp ** 7

      ! dE/dR = (- 6 * C6 * f / R^7) + (C6 / R^6 * df/dR)
      grad6 = 1.0D0 + 6.0D0 * (btemp / rtemp) ** 2
      grad6 = -6.0D0 * c6_ab(iatom,jatom) * grad6 / rtemp + &
              (c6_ab(iatom,jatom) / rtemp) * 84.0D0 /     &
              (6.0D0 * btemp / rtemp + rtemp / btemp ) ** 2
      grad6 = grad6 * dftd3_s6

      ! r^8 terms.
      rtemp = dists(iatom,jatom) ** 8
      btemp = r0_ab(iatom,jatom) ** 8

      ! dE/dR = (- 8 * C8 * f / R^9) + (C8 / R^8 * df/dR)
      grad8 = 1.0D0 + 6.0D0 * ( btemp / rtemp ) ** 2
      grad8 = -8.0D0 * c8_ab(iatom,jatom) * grad8 / rtemp + &
              (c8_ab(iatom,jatom) / rtemp) * 96.0D0 /     &
              (6.0D0 * btemp / rtemp + rtemp / btemp ) ** 2
      grad8 = dftd3_s8 * grad8 / dists(iatom,jatom)

      ! Finally, gradients are expressed as dE/dR * dR/dx
      ! and dR/dx = (xi - xj) / R
      btemp = (grad8 + grad6) / dists(iatom,jatom)

      grad(iatom,1) = btemp * (pos(iatom,1) - pos(jatom,1))
      grad(iatom,2) = btemp * (pos(iatom,2) - pos(jatom,2))
      grad(iatom,3) = btemp * (pos(iatom,3) - pos(jatom,3))
      grad(jatom,1) = btemp * (pos(jatom,1) - pos(iatom,1))
      grad(jatom,2) = btemp * (pos(jatom,2) - pos(iatom,2))
      grad(jatom,3) = btemp * (pos(jatom,3) - pos(iatom,3))
   enddo
   enddo   
end subroutine dftd3_2bodies_g