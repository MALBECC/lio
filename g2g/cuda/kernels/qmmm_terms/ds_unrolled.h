// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-S)-------------------------------------------
        scalar_type F_mU[4];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          //for (int m = 0; m <= 3; m++) 
          //{
          //  F_mU[m] = lio_gamma<scalar_type>(m,U);
          lio_gamma<scalar_type,3>(F_mU,U);
          //}
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common, d_s0, d_s1;
          uint dens_ind = 0;

          // d_l1 == 0
          {
            // d_l2 == 0
            {
              scalar_type pre_term = gpu_normalization_factor * clatom_charge_sh[j] * dens[dens_ind];
              dens_ind++;

              d_s0  = PmA[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1  = PmA[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              d_s0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // grad_l == 0
              {
                C_force_term  = PmA[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1 (d_l2)
                AB_common     = PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l2)
                C_force_term += PmA[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1 (d_l1)
                AB_common    += PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l1)
                C_force_term  = PmC[0] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * d_s0 + A_force_term;
                A_force_term += PmA[0] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l2)
                A_force_term -= PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l1)
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * 2.0f * aj * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmC[1] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[1] * d_s0 + A_force_term;
                A_force_term += PmA[1] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * 2.0f * aj * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmC[2] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[2] * d_s0 + A_force_term;
                A_force_term += PmA[2] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * 2.0f * aj * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
          // d_l1 == 1
          {
            // d_l2 == 0
            {
              scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
              dens_ind++;

              d_s0  = PmA[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1  = PmA[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              // grad_l == 0
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1 (d_l1)
                AB_common    += PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l1)
                C_force_term  = PmC[0] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * d_s0 + A_force_term;
                A_force_term += PmA[0] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l1)
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * 2.0f * aj * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmA[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1 (d_l2)
                AB_common     = PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l2)
                C_force_term  = PmC[1] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * d_s0 + A_force_term;
                A_force_term += PmA[1] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l2)
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * 2.0f * aj * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmC[2] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[2] * d_s0 + A_force_term;
                A_force_term += PmA[2] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * 2.0f * aj * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // d_l2 == 1
            {
              scalar_type pre_term = gpu_normalization_factor * clatom_charge_sh[j] * dens[dens_ind];
              dens_ind++;

              d_s0  = PmA[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              d_s1  = PmA[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
              d_s0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // grad_l == 0
              {
                C_force_term  = PmC[0] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[0] * d_s0 + A_force_term;
                A_force_term += PmA[0] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * 2.0f * aj * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmA[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1 (d_l2)
                AB_common     = PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l2)
                C_force_term += PmA[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1 (d_l1)
                AB_common    += PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l1)
                C_force_term  = PmC[1] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * d_s0 + A_force_term;
                A_force_term += PmA[1] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l2)
                A_force_term -= PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l1)
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * 2.0f * aj * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmC[2] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[2] * d_s0 + A_force_term;
                A_force_term += PmA[2] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * 2.0f * aj * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
          // d_l1 == 2
          {
            // d_l2 == 0
            {
              scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
              dens_ind++;

              d_s0  = PmA[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1  = PmA[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              // grad_l == 0
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1 (d_l1)
                AB_common    += PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l1)
                C_force_term  = PmC[0] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * d_s0 + A_force_term;
                A_force_term += PmA[0] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l1)
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * 2.0f * aj * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmC[1] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[1] * d_s0 + A_force_term;
                A_force_term += PmA[1] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * 2.0f * aj * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmA[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1 (d_l2)
                AB_common     = PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l2)
                C_force_term  = PmC[2] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * d_s0 + A_force_term;
                A_force_term += PmA[2] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (d_l2)
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * 2.0f * aj * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // d_l2 == 1
            {
              scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
              dens_ind++;

              d_s0  = PmA[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              d_s1  = PmA[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
              // grad_l == 0
              {
                C_force_term  = PmC[0] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[0] * d_s0 + A_force_term;
                A_force_term += PmA[0] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * 2.0f * aj * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1 (d_l1)
                AB_common    += PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l1)
                C_force_term  = PmC[1] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * d_s0 + A_force_term;
                A_force_term += PmA[1] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l1)
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * 2.0f * aj * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmA[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1 (d_l2)
                AB_common     = PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l2)
                C_force_term  = PmC[2] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * d_s0 + A_force_term;
                A_force_term += PmA[2] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (d_l2)
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * 2.0f * aj * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // d_l2 == 2
            {
              scalar_type pre_term = gpu_normalization_factor * clatom_charge_sh[j] * dens[dens_ind];
              dens_ind++;

              d_s0  = PmA[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
              d_s1  = PmA[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l2)
              d_s0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // grad_l == 0
              {
                C_force_term  = PmC[0] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[0] * d_s0 + A_force_term;
                A_force_term += PmA[0] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * 2.0f * aj * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmC[1] * d_s1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[1] * d_s0 + A_force_term;
                A_force_term += PmA[1] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * 2.0f * aj * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmA[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1 (d_l2)
                AB_common     = PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l2)
                C_force_term += PmA[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1 (d_l1)
                AB_common    += PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l1)
                C_force_term  = PmC[2] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * d_s0 + A_force_term;
                A_force_term += PmA[2] * d_s0;

                A_force_term *= 2.0f * ai;
                A_force_term -= PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l2)
                A_force_term -= PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (d_l1)
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * 2.0f * aj * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-S)-------------------------------------------
