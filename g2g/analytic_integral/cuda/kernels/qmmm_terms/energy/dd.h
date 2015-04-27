//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (S-S)-------------------------------------------
        scalar_type F_mU[5];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
          lio_gamma<scalar_type,4>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          scalar_type preterm1,preterm2,term;

          uint curr_ind = 0;
          for (uint d1_l1 = 0; d1_l1 < 3; d1_l1++) {

            scalar_type p_s0_1 = PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1];
            scalar_type p_s1_1 = PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2];
            scalar_type p_s2_1 = PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3];
            scalar_type p_s3_1 = PmA[d1_l1] * F_mU[3] - PmC[d1_l1] * F_mU[4];

            for (uint d1_l2 = 0; d1_l2 <= d1_l1; d1_l2++) {

              scalar_type d_s0 = PmA[d1_l2] * p_s0_1 - PmC[d1_l2] * p_s1_1;
              d_s0            += (d1_l1 == d1_l2) * inv_two_zeta * (F_mU[0] - F_mU[1]);
              scalar_type d_s1 = PmA[d1_l2] * p_s1_1 - PmC[d1_l2] * p_s2_1;
              d_s1            += (d1_l1 == d1_l2) * inv_two_zeta * (F_mU[1] - F_mU[2]);
              scalar_type d_s2 = PmA[d1_l2] * p_s2_1 - PmC[d1_l2] * p_s3_1;
              d_s2            += (d1_l1 == d1_l2) * inv_two_zeta * (F_mU[2] - F_mU[3]);

              scalar_type p_s0_2 = PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1];
              scalar_type p_s1_2 = PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2];
              scalar_type p_s2_2 = PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3];

              preterm1  = (d1_l1 != d1_l2) * 1.0f + (d1_l1 == d1_l2) * G2G::gpu_normalization_factor;

              for (uint d2_l1 = 0; d2_l1 < 3; d2_l1++) {

                scalar_type d_p0 = (d1_l1 == d2_l1) * (p_s0_2 - p_s1_2);
                d_p0            += (d1_l2 == d2_l1) * (p_s0_1 - p_s1_1);
                d_p0            *= inv_two_zeta;
                d_p0            += PmB[d2_l1] * d_s0 - PmC[d2_l1] * d_s1;
                scalar_type d_p1 = (d1_l1 == d2_l1) * (p_s1_2 - p_s2_2);
                d_p1            += (d1_l2 == d2_l1) * (p_s1_1 - p_s2_1);
                d_p1            *= inv_two_zeta;
                d_p1            += PmB[d2_l1] * d_s1 - PmC[d2_l1] * d_s2;

                scalar_type p_p0_1 = PmB[d2_l1] * p_s0_1 - PmC[d2_l1] * p_s1_1;
                p_p0_1            += (d1_l1 == d2_l1) * inv_two_zeta * (F_mU[0] - F_mU[1]);
                scalar_type p_p1_1 = PmB[d2_l1] * p_s1_1 - PmC[d2_l1] * p_s2_1;
                p_p1_1            += (d1_l1 == d2_l1) * inv_two_zeta * (F_mU[1] - F_mU[2]);

                scalar_type p_p0_2 = PmB[d2_l1] * p_s0_2 - PmC[d2_l1] * p_s1_2;
                p_p0_2            += (d1_l2 == d2_l1) * inv_two_zeta * (F_mU[0] - F_mU[1]);
                scalar_type p_p1_2 = PmB[d2_l1] * p_s1_2 - PmC[d2_l1] * p_s2_2;
                p_p1_2            += (d1_l2 == d2_l1) * inv_two_zeta * (F_mU[1] - F_mU[2]);

                for (uint d2_l2 = 0; d2_l2 <= d2_l1; d2_l2++) {

                  bool skip = same_func && (d2_l1 > d1_l1 || (d2_l1 == d1_l1 && d2_l2 > d1_l2));
                  preterm2  = (d2_l1 != d2_l2) * 1.0f + (d2_l1 == d2_l2) * G2G::gpu_normalization_factor;

                  term  = (d1_l1 == d2_l2) * (p_p0_2 - p_p1_2);
                  term += (d1_l2 == d2_l2) * (p_p0_1 - p_p1_1);
                  term += (d2_l1 == d2_l2) * (d_s0 - d_s1);
                  term *= inv_two_zeta;
                  term += PmB[d2_l2] * d_p0 - PmC[d2_l2] * d_p1;

                  my_fock[curr_ind] += !skip * preterm1 * preterm2 * clatom_charge_sh[j] * term;
                  curr_ind += !skip;
                }
              }
            }
          }
  
        }
        // END individual force terms 
//------------------------------------------END TERM-TYPE DEPENDENT PART (S-S)----------------------------------------------
