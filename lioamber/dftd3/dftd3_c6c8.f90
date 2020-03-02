! These routines set the values for C6 and C8 coefficients.
subroutine dftd3_set_c6c8(dists, n_atoms)
   use dftd3_data, only: c6_ab, c8_ab, c6_cn, c8_coef, r_cov
   implicit none
   integer     , intent(in) :: n_atoms
   LIODBLE, intent(in) :: dists(:,:)

   LIODBLE :: Lij, Wsum, Zsum, cna, cnb, c6_tmp, r_min
   integer      :: iatom, jatom, cni, cnj
   ! This is the atomic coordination number.
   LIODBLE, allocatable :: atom_cn(:)

   allocate(atom_cn(n_atoms))
   call dftd3_calc_cn(atom_cn, dists, n_atoms, r_cov)

   c6_ab  = 0.0D0
   c8_ab  = 0.0D0
   do iatom = 1    , n_atoms
   do jatom = iatom, n_atoms
      Wsum   = 0.0D0
      Zsum   = 0.0D0
      r_min  = 1.0D99
      c6_tmp = 0.0D0
      
      do cni = 1, 5
      do cnj = 1, 5
         if (c6_cn(iatom, jatom, cni, cnj, 1) > 0.0D0) then
            cna = atom_cn(iatom) - c6_cn(iatom, jatom, cni, cnj, 2)
            cnb = atom_cn(jatom) - c6_cn(iatom, jatom, cni, cnj, 3)

            Lij = cna * cna + cnb * cnb
            if (Lij < r_min) then
               r_min  = Lij
               c6_tmp = c6_cn(iatom, jatom, cni, cnj, 1)
            endif
            Lij = exp(-4.0D0 * Lij)

            Wsum = Wsum + Lij
            Zsum = Zsum + c6_cn(iatom, jatom, cni, cnj, 1) * Lij
         endif
      enddo
      enddo
      
      if (Wsum > 1.0D-99) then
         c6_ab(iatom, jatom) = Zsum / Wsum
      else
         c6_ab(iatom, jatom) = c6_tmp
      endif
      c6_ab(jatom, iatom) = c6_ab(iatom, jatom)

      c8_ab(iatom, jatom) = 3.0D0 * c6_ab(iatom, jatom) * c8_coef(iatom) *&
                                    c8_coef(jatom)
      c8_ab(jatom, iatom) = c8_ab(iatom, jatom)
   enddo
   enddo

   deallocate(atom_cn)
end subroutine dftd3_set_c6c8

subroutine dftd3_calc_cn(atom_cn, dists, n_atoms, r_cov)
   implicit none
   integer     , intent(in)    :: n_atoms
   LIODBLE, intent(in)    :: dists(:,:), r_cov(:)
   LIODBLE, intent(inout) :: atom_cn(:)

   LIODBLE :: term
   integer      :: iatom, jatom

   atom_cn = 0.0D0
   do iatom = 1, n_atoms
   do jatom = 1, n_atoms
      if (iatom /= jatom) then
         term = (r_cov(iatom) + r_cov(jatom)) / dists(iatom,jatom)
         term = 1.0D0 + exp(-16.0D0 * (term - 1.0D0) )
         atom_cn(iatom) = atom_cn(iatom) + 1.0D0 / term
      endif
   enddo
   enddo
end subroutine dftd3_calc_cn