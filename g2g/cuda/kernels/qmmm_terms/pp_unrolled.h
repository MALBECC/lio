// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (P-P)-------------------------------------------
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
          scalar_type AB_common, p_p0, p_p1;
          scalar_type mm_charge = clatom_charge_sh[j];
          uint dens_ind = 0;

          // p1_l == 0
          {
            // p2_l == 0
            {
              scalar_type pre_term;
              {
                pre_term  = mm_charge * dens[dens_ind];
                dens_ind++;
              }

              p_p0  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0
              p_p0 -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1
              p_p1  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1
              p_p1 -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2
              p_p0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              p_p1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // grad_l == 0
              {
                C_force_term  = PmB[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (B)
                C_force_term += PmA[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1
                AB_common    += PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0
                C_force_term  = PmC[0] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (B)
                B_force_term -= PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmC[1] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmC[2] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // p2_l == 1
            {
              scalar_type pre_term;
              {
                pre_term  = !same_func * mm_charge * dens[dens_ind];
                dens_ind += !same_func;
              }

              p_p0  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0
              p_p0 -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1
              p_p1  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1
              p_p1 -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2
              // grad_l == 0
              {
                C_force_term  = PmB[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (B)
                C_force_term  = PmC[0] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (B)
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1
                AB_common    += PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0
                C_force_term  = PmC[1] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                B_force_term -= PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmC[2] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // p2_l == 2
            {
              scalar_type pre_term;
              {
                pre_term  = !same_func * mm_charge * dens[dens_ind];
                dens_ind += !same_func;
              }

              p_p0  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0
              p_p0 -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1
              p_p1  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1
              p_p1 -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2
              // grad_l == 0
              {
                C_force_term  = PmB[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (B)
                C_force_term  = PmC[0] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (B)
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmC[1] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1
                AB_common    += PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0
                C_force_term  = PmC[2] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                B_force_term -= PmA[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
          // p1_l == 1
          {
            // p2_l == 0
            {
              scalar_type pre_term;
              {
                pre_term  = mm_charge * dens[dens_ind];
                dens_ind++;
              }

              p_p0  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0
              p_p0 -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1
              p_p1  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1
              p_p1 -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2
              // grad_l == 0
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1
                AB_common    += PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0
                C_force_term  = PmC[0] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                B_force_term -= PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmB[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (B)
                C_force_term  = PmC[1] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (B)
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmC[2] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // p2_l == 1
            {
              scalar_type pre_term;
              {
                pre_term  = mm_charge * dens[dens_ind];
                dens_ind++;
              }

              p_p0  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0
              p_p0 -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1
              p_p1  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1
              p_p1 -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2
              p_p0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              p_p1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // grad_l == 0
              {
                C_force_term  = PmC[0] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmB[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (B)
                C_force_term += PmA[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1
                AB_common    += PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0
                C_force_term  = PmC[1] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (B)
                B_force_term -= PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmC[2] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // p2_l == 2
            {
              scalar_type pre_term;
              {
                pre_term  = !same_func * mm_charge * dens[dens_ind];
                dens_ind += !same_func;
              }

              p_p0  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0
              p_p0 -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1
              p_p1  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1
              p_p1 -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2
              // grad_l == 0
              {
                C_force_term  = PmC[0] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmB[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (B)
                C_force_term  = PmC[1] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (B)
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1
                AB_common    += PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0
                C_force_term  = PmC[2] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                B_force_term -= PmA[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
          // p1_l == 2
          {
            // p2_l == 0
            {
              scalar_type pre_term;
              {
                pre_term  = mm_charge * dens[dens_ind];
                dens_ind++;
              }

              p_p0  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0
              p_p0 -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1
              p_p1  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1
              p_p1 -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2
              // grad_l == 0
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1
                AB_common    += PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0
                C_force_term  = PmC[0] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                B_force_term -= PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmC[1] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmB[0] * F_mU[1] - PmC[0] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (B)
                C_force_term  = PmC[2] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[0] * F_mU[0] - PmC[0] * F_mU[1]; // p_s0 (B)
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // p2_l == 1
            {
              scalar_type pre_term;
              {
                pre_term  = mm_charge * dens[dens_ind];
                dens_ind++;
              }

              p_p0  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0
              p_p0 -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1
              p_p1  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1
              p_p1 -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2
              // grad_l == 0
              {
                C_force_term  = PmC[0] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
                C_force_term += PmA[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1
                AB_common    += PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0
                C_force_term  = PmC[1] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                B_force_term -= PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmB[1] * F_mU[1] - PmC[1] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (B)
                C_force_term  = PmC[2] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[1] * F_mU[0] - PmC[1] * F_mU[1]; // p_s0 (B)
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
            // p2_l == 2
            {
              scalar_type pre_term;
              {
                pre_term  = mm_charge * dens[dens_ind];
                dens_ind++;
              }

              p_p0  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0
              p_p0 -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1
              p_p1  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1
              p_p1 -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2
              p_p0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              p_p1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // grad_l == 0
              {
                C_force_term  = PmC[0] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[0] * p_p0 + A_force_term;
                A_force_term += PmA[0] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[0]     += pre_term * A_force_term;
                B_force[0]     += pre_term * B_force_term;
                C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 1
              {
                C_force_term  = PmC[1] * p_p1;

                A_force_term  = -C_force_term;
                B_force_term  = PmB[1] * p_p0 + A_force_term;
                A_force_term += PmA[1] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force[1]     += pre_term * A_force_term;
                B_force[1]     += pre_term * B_force_term;
                C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
              // grad_l == 2
              {
                C_force_term  = PmB[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (B)
                C_force_term += PmA[2] * F_mU[1] - PmC[2] * F_mU[2]; // p_s1
                AB_common    += PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0
                C_force_term  = PmC[2] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[2] * p_p0 + A_force_term;
                A_force_term += PmA[2] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
                A_force_term -= PmB[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0 (B)
                B_force_term -= PmA[2] * F_mU[0] - PmC[2] * F_mU[1]; // p_s0
                A_force[2]     += pre_term * A_force_term;
                B_force[2]     += pre_term * B_force_term;
                C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (P-P)-------------------------------------------
