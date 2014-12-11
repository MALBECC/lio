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

          for (int orb1 = 0; orb1 < 3; orb1++)
          {
            //scalar_type p_s0 = PmA[orb1] * F_mU[0] - PmC[orb1] * F_mU[1];
            //scalar_type p_s1 = PmA[orb1] * F_mU[1] - PmC[orb1] * F_mU[2];
            for (int orb2 = 0; orb2 < 3; orb2++)
            {
              bool del12 = orb1 == orb2;
              scalar_type pre_term;
              {
                bool skip = same_func && (orb2 > orb1);
                pre_term = (1-skip) * mm_charge * dens[(orb1*3+orb2)];
              }

              p_p0  = PmB[orb2] * (PmA[orb1] * F_mU[0] - PmC[orb1] * F_mU[1]); // p_s0
              p_p0 -= PmC[orb2] * (PmA[orb1] * F_mU[1] - PmC[orb1] * F_mU[2]); // p_s1
              p_p0 += del12 * inv_two_zeta * (F_mU[0] - F_mU[1]);
              p_p1  = PmB[orb2] * (PmA[orb1] * F_mU[1] - PmC[orb1] * F_mU[2]); // p_s1
              p_p1 -= PmC[orb2] * (PmA[orb1] * F_mU[2] - PmC[orb1] * F_mU[3]); // p_s2
              p_p1 += del12 * inv_two_zeta * (F_mU[1] - F_mU[2]);
              for (int grad_l = 0; grad_l < 3; grad_l++)
              {
                bool del1l = orb1 == grad_l; bool del2l = orb2 == grad_l;
                C_force_term  = del1l * (PmB[orb2] * F_mU[1] - PmC[orb2] * F_mU[2]); // p_s1 (B)
                C_force_term += del2l * (PmA[orb1] * F_mU[1] - PmC[orb1] * F_mU[2]); // p_s1
                C_force_term  = PmC[grad_l] * p_p1 + inv_two_zeta * C_force_term;

                AB_common     = del1l * (PmB[orb2] * F_mU[0] - PmC[orb2] * F_mU[1]); // p_s0 (B)
                AB_common    += del2l * (PmA[orb1] * F_mU[0] - PmC[orb1] * F_mU[1]); // p_s0

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[grad_l] * p_p0 + A_force_term;
                A_force_term += PmA[grad_l] * p_p0;

                A_force[grad_l]     += pre_term * (2.0f * ai * A_force_term - del1l * (PmB[orb2] * F_mU[0] - PmC[orb2] * F_mU[1])); // p_s0 (B)
                B_force[grad_l]     += pre_term * (2.0f * aj * B_force_term - del2l * (PmA[orb1] * F_mU[0] - PmC[orb1] * F_mU[1])); // p_s0
                C_force[grad_l][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (P-P)-------------------------------------------
