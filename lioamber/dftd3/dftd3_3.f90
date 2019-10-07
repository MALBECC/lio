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

subroutine dftd3_3bodies_g(grad, dists, pos, n_atoms)
        use dftd3_data, only: c6_ab, r0_ab
        implicit none
        integer     , intent(in)    :: n_atoms
        real(kind=8), intent(in)    :: dists(:,:), pos(:,:)
        real(kind=8), intent(inout) :: grad(:,:)
     
        real(kind=8) :: coef, dsq1, dsq2, dsq3, td1, td2, td3, dist3, &
                        r0_abc, fdamp, fang, c9_abc, term_ab, term_ac

        integer      :: iatom, jatom, katom
     
        do iatom = 1      , n_atoms
        do jatom = iatom+1, n_atoms
        do katom = jatom+1, n_atoms

           ! Computes front coefficients.
           dist3  = dists(iatom,jatom) * dists(iatom,katom) * dists(jatom,katom)
           dist3  = dist3 ** (1.0D0/3.0D0)
           r0_abc = r0_ab(iatom,jatom) * r0_ab(iatom,katom) * r0_ab(jatom,katom)
           r0_abc = (4.0D0 / 3.0D0) * r0_abc ** (1.0D0/3.0D0)
           coef   = (r0_abc / dist3) ** 8
           c9_abc = - sqrt(c6_ab(iatom,jatom) * c6_ab(iatom,katom) *&
                           c6_ab(jatom,katom))

           ! Since fdamp will be used in explicit dE/dRab derivatives,
           ! which only involve the angular terms, the 3/4 factor is
           ! already included here.
           fdamp  = 0.75D0 * (c9_abc / dist3 ** 15) / &
                             (1.0D0 + 6.0D0 * (coef * coef))

           ! Computes cosine angular terms.
           dsq1 = dists(iatom,jatom) * dists(iatom,jatom)
           dsq2 = dists(iatom,katom) * dists(iatom,katom)
           dsq3 = dists(jatom,katom) * dists(jatom,katom)
           
           td1  = dsq1 + dsq2 - dsq3
           td2  = dsq1 + dsq3 - dsq2
           td3  = dsq2 + dsq3 - dsq1
     
           fang = 0.375D0 * td1 * td2 * td3 / (dsq1 * dsq2 * dsq3)
        
           ! Gradient due to explicit dE/dRab and dE/dRac derivatives.
           term_ab = ((dsq2 - dsq3) ** 2 + 2.0D0 * dsq1 * (dsq2 + dsq3) - &
                      3.0D0 * dsq1 ** 2) * dists(iatom,jatom) * fdamp
           term_ac = ((dsq1 - dsq3) ** 2 + 2.0D0 * dsq2 * (dsq1 + dsq3) - &
                      3.0D0 * dsq2 ** 2) * dists(iatom,katom) * fdamp
           
           ! Computes dE/dRabc. fdamp is now used just as temporary storage, 
           ! and includes the Rabc/3 term due to dRabc/dRab = Rabc/3Rab
           fdamp = (c9_abc / (dist3 * (r0_abc ** 8))) / &
                   (18.0D0 * coef + 3.0D0 / coef)
           fdamp = fdamp * (fang * (2.0D0 * coef - 5.0D0 / coef) + &
                            14.0D0 * coef - 3.0D0 / coef)
           term_ab = term_ab + fdamp / dists(iatom,jatom)
           term_ab = term_ab / dists(iatom,jatom)
           term_ac = term_ac + fdamp / dists(iatom,katom)
           term_ac = term_ac / dists(iatom,katom)
           
           grad(iatom,1) = term_ab * (pos(iatom,1) - pos(jatom,1))
           grad(iatom,2) = term_ab * (pos(iatom,2) - pos(jatom,2))
           grad(iatom,3) = term_ab * (pos(iatom,3) - pos(jatom,3))
           grad(jatom,1) = term_ab * (pos(jatom,1) - pos(iatom,1))
           grad(jatom,2) = term_ab * (pos(jatom,2) - pos(iatom,2))
           grad(jatom,3) = term_ab * (pos(jatom,3) - pos(iatom,3))

           grad(iatom,1) = term_ac * (pos(iatom,1) - pos(katom,1))
           grad(iatom,2) = term_ac * (pos(iatom,2) - pos(katom,2))
           grad(iatom,3) = term_ac * (pos(iatom,3) - pos(katom,3))
           grad(katom,1) = term_ac * (pos(katom,1) - pos(iatom,1))
           grad(katom,2) = term_ac * (pos(katom,2) - pos(iatom,2))
           grad(katom,3) = term_ac * (pos(katom,3) - pos(iatom,3))
        enddo
        enddo
        enddo   
end subroutine dftd3_3bodies_g