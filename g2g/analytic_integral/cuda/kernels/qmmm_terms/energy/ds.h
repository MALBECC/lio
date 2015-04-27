//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (S-S)-------------------------------------------
        scalar_type F_mU[3];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
          lio_gamma<scalar_type,2>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          scalar_type preterm,term;

          uint curr_ind = 0;
          for (uint d_l1 = 0; d_l1 < 3; d_l1++) {
            scalar_type p_s0 = PmA[d_l1] * F_mU[0] - PmC[d_l1] * F_mU[1];
            scalar_type p_s1 = PmA[d_l1] * F_mU[1] - PmC[d_l1] * F_mU[2];
            for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {

              preterm  = (d_l1 != d_l2) * 1.0f + (d_l1 == d_l2) * G2G::gpu_normalization_factor;
              term     = PmA[d_l2] * p_s0 - PmC[d_l2] * p_s1;
              term    += (d_l1 == d_l2) * inv_two_zeta * (F_mU[0] - F_mU[1]);

              my_fock[curr_ind] += preterm * clatom_charge_sh[j] * term;
              curr_ind++;
            }
          }
  
        }
        // END individual force terms 
//------------------------------------------END TERM-TYPE DEPENDENT PART (S-S)----------------------------------------------
