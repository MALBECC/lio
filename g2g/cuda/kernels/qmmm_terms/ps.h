//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (P-S)-------------------------------------------
        scalar_type F_mU[3];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * zeta;
          //F_mU[0] = (SQRT_PI / (2*sqrtU)) * erff(sqrtU);
          for (int m = 0; m <= 2; m++) 
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
          //scalar_type mm_charge = clatom_charge_sh[j];

          for (int p_l = 0; p_l < 3; p_l++)
          {
            scalar_type pre_term = /*mm_charge*/clatom_charge_sh[j] * dens[p_l];
            scalar_type p_s0 = PmA[p_l] * F_mU[0] - PmC[p_l] * F_mU[1];
            scalar_type p_s1 = PmA[p_l] * F_mU[1] - PmC[p_l] * F_mU[2];
            for (int grad_l = 0; grad_l < 3; grad_l++)
            {
              bool del = p_l == grad_l;
              C_force_term  = PmC[grad_l] * p_s1;//(PmA[p_l] * F_mU[1] - PmC[p_l] * F_mU[2]) + del * inv_two_zeta * F_mU[1];
              C_force_term += del * inv_two_zeta * F_mU[1];

              A_force_term  = del * inv_two_zeta * F_mU[0] - C_force_term;
              B_force_term  = PmB[grad_l] * p_s0 + A_force_term;
              A_force_term += PmA[grad_l] * p_s0;

              A_force[grad_l]     += pre_term * (2.0f * ai * A_force_term - del * F_mU[0]);
              B_force[grad_l]     += pre_term * (2.0f * aj * B_force_term);
              C_force[grad_l][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (P-S)-------------------------------------------
