// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (P-S)-------------------------------------------
        scalar_type F_mU[3];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * zeta;
          for (int m = 0; m <= 2; m++) 
          {
            F_mU[m] = lio_gamma<scalar_type>(m,U);
          }
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          // p_l == 0
          {
            scalar_type pre_term = clatom_charge_sh[j] * dens[0];
            scalar_type p_s0 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
            scalar_type p_s1 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
            // grad_l == 0
            {
              C_force_term  = PmC[0] * p_s1;
              C_force_term += inv_two_zeta * F_mU[1];

              A_force_term  = inv_two_zeta * F_mU[0] - C_force_term;
              B_force_term  = PmB[0] * p_s0 + A_force_term;
              A_force_term += PmA[0] * p_s0;

              A_force[0]     += pre_term * (2.0f * ai * A_force_term - F_mU[0]);
              B_force[0]     += pre_term * (2.0f * aj * B_force_term);
              C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
            // grad_l == 1
            {
              C_force_term  = PmC[1] * p_s1;

              A_force_term  = -C_force_term;
              B_force_term  = PmB[1] * p_s0 + A_force_term;
              A_force_term += PmA[1] * p_s0;

              A_force[1]     += pre_term * (2.0f * ai * A_force_term);
              B_force[1]     += pre_term * (2.0f * aj * B_force_term);
              C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
            // grad_l == 2
            {
              C_force_term  = PmC[2] * p_s1;

              A_force_term  = -C_force_term;
              B_force_term  = PmB[2] * p_s0 + A_force_term;
              A_force_term += PmA[2] * p_s0;

              A_force[2]     += pre_term * (2.0f * ai * A_force_term);
              B_force[2]     += pre_term * (2.0f * aj * B_force_term);
              C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
          }
          // p_l == 1
          {
            scalar_type pre_term = clatom_charge_sh[j] * dens[1];
            scalar_type p_s0 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
            scalar_type p_s1 = PmA[1] * F_mU[1] - PmC[1] * F_mU[2];
            // grad_l == 0
            {
              C_force_term  = PmC[0] * p_s1;

              A_force_term  = -C_force_term;
              B_force_term  = PmB[0] * p_s0 + A_force_term;
              A_force_term += PmA[0] * p_s0;

              A_force[0]     += pre_term * (2.0f * ai * A_force_term);
              B_force[0]     += pre_term * (2.0f * aj * B_force_term);
              C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
            // grad_l == 1
            {
              C_force_term  = PmC[1] * p_s1;
              C_force_term += inv_two_zeta * F_mU[1];

              A_force_term  = inv_two_zeta * F_mU[0] - C_force_term;
              B_force_term  = PmB[1] * p_s0 + A_force_term;
              A_force_term += PmA[1] * p_s0;

              A_force[1]     += pre_term * (2.0f * ai * A_force_term - F_mU[0]);
              B_force[1]     += pre_term * (2.0f * aj * B_force_term);
              C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
            // grad_l == 2
            {
              C_force_term  = PmC[2] * p_s1;

              A_force_term  = -C_force_term;
              B_force_term  = PmB[2] * p_s0 + A_force_term;
              A_force_term += PmA[2] * p_s0;

              A_force[2]     += pre_term * (2.0f * ai * A_force_term);
              B_force[2]     += pre_term * (2.0f * aj * B_force_term);
              C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
          }
          // p_l == 2
          {
            scalar_type pre_term = clatom_charge_sh[j] * dens[2];
            scalar_type p_s0 = PmA[2] * F_mU[0] - PmC[2] * F_mU[1];
            scalar_type p_s1 = PmA[2] * F_mU[1] - PmC[2] * F_mU[2];
            // grad_l == 0
            {
              C_force_term  = PmC[0] * p_s1;

              A_force_term  = -C_force_term;
              B_force_term  = PmB[0] * p_s0 + A_force_term;
              A_force_term += PmA[0] * p_s0;

              A_force[0]     += pre_term * (2.0f * ai * A_force_term);
              B_force[0]     += pre_term * (2.0f * aj * B_force_term);
              C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
            // grad_l == 1
            {
              C_force_term  = PmC[1] * p_s1;

              A_force_term  = -C_force_term;
              B_force_term  = PmB[1] * p_s0 + A_force_term;
              A_force_term += PmA[1] * p_s0;

              A_force[1]     += pre_term * (2.0f * ai * A_force_term);
              B_force[1]     += pre_term * (2.0f * aj * B_force_term);
              C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
            // grad_l == 2
            {
              C_force_term  = PmC[2] * p_s1;
              C_force_term += inv_two_zeta * F_mU[1];

              A_force_term  = inv_two_zeta * F_mU[0] - C_force_term;
              B_force_term  = PmB[2] * p_s0 + A_force_term;
              A_force_term += PmA[2] * p_s0;

              A_force[2]     += pre_term * (2.0f * ai * A_force_term - F_mU[0]);
              B_force[2]     += pre_term * (2.0f * aj * B_force_term);
              C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (P-S)-------------------------------------------
