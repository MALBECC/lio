! 3-body terms for DFTD3 corrections to energy and gradients.
subroutine dftd3_3bodies_e(e_3, dists, n_atoms)
   use dftd3_data, only: c6_ab, r0_ab
   implicit none
   integer     , intent(in)    :: n_atoms
   real(kind=8), intent(in)    :: dists(:,:)
   real(kind=8), intent(inout) :: e_3

   real(kind=8) :: coef, term, dsq1, dsq2, dsq3, td1, td2, td3, dist3
   integer      :: iatom, jatom, katom

   do iatom = 1      , n_atoms
   do jatom = iatom+1, n_atoms
   do katom = jatom+1, n_atoms
      ! Computes front coefficients.
      dist3 = dists(iatom,jatom) * dists(iatom,katom) * dists(jatom,katom)
      coef  = r0_ab(iatom,jatom) * r0_ab(iatom,katom) * r0_ab(jatom,katom) / &
              dist3
      coef  = 1.0D0 + 6.0D0 * ( (4.0D0 / 3.0D0) * coef ** (1.0D0/3.0D0) ) ** 16
      coef  = sqrt(c6_ab(iatom,jatom)*c6_ab(iatom,katom)* c6_ab(jatom,katom)) /&
              coef
      
      ! Computes cosine angular terms.
      dsq1 = dists(iatom,jatom) * dists(iatom,jatom)
      dsq2 = dists(iatom,katom) * dists(iatom,katom)
      dsq3 = dists(jatom,katom) * dists(jatom,katom)
      
      td1  = dsq1 + dsq2 - dsq3
      td2  = dsq1 + dsq3 - dsq2
      td3  = dsq2 + dsq3 - dsq1

      term = 1.0D0 + 0.375D0 * td1 * td2 * td3 / (dist3 * dist3)
      term = term / (dist3 * dist3 * dist3)

      e_3 = e_3 + coef * term
   enddo
   enddo
   enddo   
end subroutine dftd3_3bodies_e