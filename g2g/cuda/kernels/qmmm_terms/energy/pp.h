//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (S-S)-------------------------------------------
        scalar_type F_mU[3];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
          lio_gamma<scalar_type,2>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          scalar_type term;

          uint curr_ind = 0;
          for (uint p1_l = 0; p1_l < 3; p1_l++) {
            scalar_type p_s0 = PmA[p1_l] * F_mU[0] - PmC[p1_l] * F_mU[1];
            scalar_type p_s1 = PmA[p1_l] * F_mU[1] - PmC[p1_l] * F_mU[2];
            for (uint p2_l = 0; p2_l < 3; p2_l++) {
              bool skip = same_func && (p2_l > p1_l);

              term  = PmB[p2_l] * p_s0 - PmC[p2_l] * p_s1;
              term += (p1_l == p2_l) * inv_two_zeta * (F_mU[0] - F_mU[1]);

              my_fock[curr_ind] += !skip * clatom_charge_sh[j] * term;
              curr_ind += !skip;
            }
          }
  
        }
        // END individual force terms 
//------------------------------------------END TERM-TYPE DEPENDENT PART (S-S)----------------------------------------------
