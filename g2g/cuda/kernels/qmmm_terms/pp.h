//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (P-P)-------------------------------------------
        scalar_type F_mU[4];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * zeta;
          //F_mU[0] = (SQRT_PI / (2*sqrtU)) * erff(sqrtU);
          for (int m = 0; m <= 3; m++) 
          {
            // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
            F_mU[m] = lio_gamma<scalar_type>(m,U);
            //F_mU[m] = fetch(qmmm_F_values_tex,(float)(U/gamma_inc-0.5f),(float)(m+0.5f));
          }
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common, p_p0, p_p1;
          scalar_type mm_charge = clatom_charge_sh[j];
          uint dens_ind = 0;

          for (int p1_l = 0; p1_l < 3; p1_l++)
          {
            //scalar_type p_s0 = PmA[p1_l] * F_mU[0] - PmC[p1_l] * F_mU[1];
            //scalar_type p_s1 = PmA[p1_l] * F_mU[1] - PmC[p1_l] * F_mU[2];
            for (int p2_l = 0; p2_l < 3; p2_l++)
            {
              bool del12 = p1_l == p2_l;
              scalar_type pre_term;
              {
                bool skip = same_func && (p2_l > p1_l);
                pre_term  = !skip * mm_charge * dens[dens_ind];
                dens_ind += !skip;
              }

              p_p0  = PmB[p2_l] * (PmA[p1_l] * F_mU[0] - PmC[p1_l] * F_mU[1]); // p_s0
              p_p0 -= PmC[p2_l] * (PmA[p1_l] * F_mU[1] - PmC[p1_l] * F_mU[2]); // p_s1
              p_p0 += del12 * inv_two_zeta * (F_mU[0] - F_mU[1]);
              p_p1  = PmB[p2_l] * (PmA[p1_l] * F_mU[1] - PmC[p1_l] * F_mU[2]); // p_s1
              p_p1 -= PmC[p2_l] * (PmA[p1_l] * F_mU[2] - PmC[p1_l] * F_mU[3]); // p_s2
              p_p1 += del12 * inv_two_zeta * (F_mU[1] - F_mU[2]);
              for (int grad_l = 0; grad_l < 3; grad_l++)
              {
                bool del1l = p1_l == grad_l; bool del2l = p2_l == grad_l;
                C_force_term  = del1l * (PmB[p2_l] * F_mU[1] - PmC[p2_l] * F_mU[2]); // p_s1 (B)
                C_force_term += del2l * (PmA[p1_l] * F_mU[1] - PmC[p1_l] * F_mU[2]); // p_s1
                C_force_term  = PmC[grad_l] * p_p1 + inv_two_zeta * C_force_term;

                AB_common     = del1l * (PmB[p2_l] * F_mU[0] - PmC[p2_l] * F_mU[1]); // p_s0 (B)
                AB_common    += del2l * (PmA[p1_l] * F_mU[0] - PmC[p1_l] * F_mU[1]); // p_s0

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[grad_l] * p_p0 + A_force_term;
                A_force_term += PmA[grad_l] * p_p0;

                A_force[grad_l]     += pre_term * (2.0f * ai * A_force_term - del1l * (PmB[p2_l] * F_mU[0] - PmC[p2_l] * F_mU[1])); // p_s0 (B)
                B_force[grad_l]     += pre_term * (2.0f * aj * B_force_term - del2l * (PmA[p1_l] * F_mU[0] - PmC[p1_l] * F_mU[1])); // p_s0
                C_force[grad_l][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (P-P)-------------------------------------------
