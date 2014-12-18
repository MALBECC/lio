//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (S-S)-------------------------------------------
        scalar_type F_mU[4];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
          lio_gamma<scalar_type,3>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          scalar_type preterm,term;

          uint curr_ind = 0;
          for (uint d_l1 = 0; d_l1 < 3; d_l1++) {
            scalar_type p_s0_1 = PmA[d_l1] * F_mU[0] - PmC[d_l1] * F_mU[1];
            scalar_type p_s1_1 = PmA[d_l1] * F_mU[1] - PmC[d_l1] * F_mU[2];
            scalar_type p_s2_1 = PmA[d_l1] * F_mU[2] - PmC[d_l1] * F_mU[3];
            for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {
              scalar_type d_s0 = PmA[d_l2] * p_s0_1 - PmC[d_l2] * p_s1_1;
              d_s0            += (d_l1 == d_l2) * inv_two_zeta * (F_mU[0] - F_mU[1]);
              scalar_type d_s1 = PmA[d_l2] * p_s1_1 - PmC[d_l2] * p_s2_1;
              d_s1            += (d_l1 == d_l2) * inv_two_zeta * (F_mU[1] - F_mU[2]);

              scalar_type p_s0_2 = PmA[d_l2] * F_mU[0] - PmC[d_l2] * F_mU[1];
              scalar_type p_s1_2 = PmA[d_l2] * F_mU[1] - PmC[d_l2] * F_mU[2];
              preterm  = (d_l1 != d_l2) * 1.0f + (d_l1 == d_l2) * gpu_normalization_factor;
              for (uint p_l = 0; p_l < 3; p_l++) {

                term  = PmB[p_l] * d_s0 - PmC[p_l] * d_s1;
                term += (d_l1 == p_l) * inv_two_zeta * (p_s0_2 - p_s1_2);
                term += (d_l2 == p_l) * inv_two_zeta * (p_s0_1 - p_s1_1);

                my_fock[curr_ind] += preterm * clatom_charge_sh[j] * term;
                curr_ind++;
              }
            }
          }
  
        }
        // END individual force terms 
//------------------------------------------END TERM-TYPE DEPENDENT PART (S-S)----------------------------------------------
