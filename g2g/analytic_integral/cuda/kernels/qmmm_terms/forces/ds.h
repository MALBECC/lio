//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-S)-------------------------------------------
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
          scalar_type AB_common, d_s0, d_s1;
          //scalar_type mm_charge = clatom_charge_sh[j];
          uint dens_ind = 0;

          for (int d_l1 = 0; d_l1 < 3; d_l1++)
          {
            for (int d_l2 = 0; d_l2 <= d_l1; d_l2++)
            {
              bool del12 = d_l1 == d_l2;
              scalar_type pre_term = /*mm_charge*/clatom_charge_sh[j] * dens[dens_ind];
              pre_term *= !del12*1.0f + del12*G2G::gpu_normalization_factor;
              //if (valid_thread) printf("%.4e\n", dens[dens_ind]);
              dens_ind++;

              d_s0  = PmA[d_l1] * (PmA[d_l2] * F_mU[0] - PmC[d_l2] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[d_l1] * (PmA[d_l2] * F_mU[1] - PmC[d_l2] * F_mU[2]); // p_s1 (d_l2)
              d_s0 += del12 * inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1  = PmA[d_l1] * (PmA[d_l2] * F_mU[1] - PmC[d_l2] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[d_l1] * (PmA[d_l2] * F_mU[2] - PmC[d_l2] * F_mU[3]); // p_s2 (d_l2)
              d_s1 += del12 * inv_two_zeta * (F_mU[1] - F_mU[2]);
              for (int grad_l = 0; grad_l < 3; grad_l++)
              {
                bool del1l = d_l1 == grad_l; bool del2l = d_l2 == grad_l;
                C_force_term  = del1l * (PmA[d_l2] * F_mU[1] - PmC[d_l2] * F_mU[2]); // p_s1 (d_l2)
                C_force_term += del2l * (PmA[d_l1] * F_mU[1] - PmC[d_l1] * F_mU[2]); // p_s1 (d_l1)
                C_force_term  = PmC[grad_l] * d_s1 + inv_two_zeta * C_force_term;

                AB_common     = del1l * (PmA[d_l2] * F_mU[0] - PmC[d_l2] * F_mU[1]); // p_s0 (d_l2)
                AB_common    += del2l * (PmA[d_l1] * F_mU[0] - PmC[d_l1] * F_mU[1]); // p_s0 (d_l1)

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[grad_l] * d_s0 + A_force_term;
                A_force_term += PmA[grad_l] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= del1l * (PmA[d_l2] * F_mU[0] - PmC[d_l2] * F_mU[1]); // p_s0 (d_l2)
                A_force_term -= del2l * (PmA[d_l1] * F_mU[0] - PmC[d_l1] * F_mU[1]); // p_s0 (d_l1)
                A_force[grad_l]     += pre_term * A_force_term;
                B_force[grad_l]     += pre_term * 2.0f * aj * B_force_term;
                C_force[grad_l][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-S)-------------------------------------------
