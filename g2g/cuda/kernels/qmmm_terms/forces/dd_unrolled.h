// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-D)-------------------------------------------
        scalar_type F_mU[6];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          //F_mU[0] = (SQRT_PI / (2*sqrtU)) * erff(sqrtU);
          //for (int m = 0; m <= 5; m++) 
          //{
            // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
          //  F_mU[m] = lio_gamma<scalar_type>(m,U);
          lio_gamma<scalar_type,5>(F_mU,U);
            //F_mU[m] = fetch(qmmm_F_values_tex,(float)(U/gamma_inc-0.5f),(float)(m+0.5f));
          //}
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common;
          //scalar_type mm_charge = clatom_charge_sh[j];
          uint dens_ind = 0;

          {
            {

              scalar_type d1_s0  = PmA[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
              d1_s0             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
              scalar_type d1_s1  = PmA[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
              d1_s1             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d1_l2)
              scalar_type d1_s2  = PmA[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d1_l2)
              d1_s2             -= PmC[0] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s3 (d1_l2)
              scalar_type d1_s3  = PmA[0] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s2 (d1_l2)
              d1_s3             -= PmC[0] * (PmA[0] * F_mU[4] - PmC[0] * F_mU[5]); // p_s3 (d1_l2)
              d1_s0             += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d1_s1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d1_s2             += inv_two_zeta * (F_mU[2] - F_mU[3]);
              d1_s3             += inv_two_zeta * (F_mU[3] - F_mU[4]);
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[0] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l1)
                p_p0_d1l1_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l1_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l1_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1  = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p1_d2l1  = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p2_d2l1  = (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]) - (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p0_d2l1 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p1_d2l1 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p2_d2l1 += (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]) - (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                d1_p1_d2l1 += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                d1_p2_d2l1 += PmB[0] * d1_s2 - PmC[0] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1  = (PmB[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmB[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmB[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1  = p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[1] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l1)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                d1_p1_d2l1 += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                d1_p2_d2l1 += PmB[1] * d1_s2 - PmC[1] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2   = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p1_d2l2   = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1  = (PmB[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmB[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1  = p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[2] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l1)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[2] * d1_s0 - PmC[2] * d1_s1;
                d1_p1_d2l1 += PmB[2] * d1_s1 - PmC[2] * d1_s2;
                d1_p2_d2l1 += PmB[2] * d1_s2 - PmC[2] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2   = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p1_d2l2   = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1  = p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                  d1_p1_d2l2  += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[2] * p_p0_d1l1_d2l1 - PmC[2] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[2] * p_p1_d1l1_d2l1 - PmC[2] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[2] * d1_p0_d2l1 - PmC[2] * d1_p1_d2l1;
                  d_d1 += PmB[2] * d1_p1_d2l1 - PmC[2] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
            }
          }
          {
            {

              scalar_type d1_s0  = PmA[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
              d1_s0             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
              scalar_type d1_s1  = PmA[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
              d1_s1             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d1_l2)
              scalar_type d1_s2  = PmA[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d1_l2)
              d1_s2             -= PmC[1] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s3 (d1_l2)
              scalar_type d1_s3  = PmA[1] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s2 (d1_l2)
              d1_s3             -= PmC[1] * (PmA[0] * F_mU[4] - PmC[0] * F_mU[5]); // p_s3 (d1_l2)
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[0] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l1)
                scalar_type p_p0_d1l2_d2l1  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[0] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l2)
                p_p0_d1l2_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l2_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l2_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p1_d2l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p2_d2l1 += (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]) - (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                d1_p1_d2l1 += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                d1_p2_d2l1 += PmB[0] * d1_s2 - PmC[0] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l2  = (PmB[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmB[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmB[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[1] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l1)
                p_p0_d1l1_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l1_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l1_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);
                scalar_type p_p0_d1l2_d2l1  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[1] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l2)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1  = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p1_d2l1  = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p2_d2l1  = (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]) - (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                d1_p1_d2l1 += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                d1_p2_d2l1 += PmB[1] * d1_s2 - PmC[1] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2  = (PmB[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmB[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1  = (PmB[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmB[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[1] * p_p0_d1l2_d2l1 - PmC[1] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[1] * p_p1_d1l2_d2l1 - PmC[1] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l2_d2l1 - p_p1_d1l2_d2l1;
                  d_d1  = p_p1_d1l2_d2l1 - p_p2_d1l2_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[2] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l1)
                scalar_type p_p0_d1l2_d2l1  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[2] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l2)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[2] * d1_s0 - PmC[2] * d1_s1;
                d1_p1_d2l1 += PmB[2] * d1_s1 - PmC[2] * d1_s2;
                d1_p2_d2l1 += PmB[2] * d1_s2 - PmC[2] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2   = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p1_d2l2   = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                  d1_p1_d2l2  += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[1] * p_p0_d1l2_d2l1 - PmC[1] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[1] * p_p1_d1l2_d2l1 - PmC[1] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l2_d2l1 - p_p1_d1l2_d2l1;
                  d_d1  = p_p1_d1l2_d2l1 - p_p2_d1l2_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[2] * p_p0_d1l1_d2l1 - PmC[2] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[2] * p_p1_d1l1_d2l1 - PmC[2] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[2] * p_p0_d1l2_d2l1 - PmC[2] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[2] * p_p1_d1l2_d2l1 - PmC[2] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[2] * d1_p0_d2l1 - PmC[2] * d1_p1_d2l1;
                  d_d1 += PmB[2] * d1_p1_d2l1 - PmC[2] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
            }
            {

              scalar_type d1_s0  = PmA[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l2)
              d1_s0             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l2)
              scalar_type d1_s1  = PmA[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l2)
              d1_s1             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d1_l2)
              scalar_type d1_s2  = PmA[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d1_l2)
              d1_s2             -= PmC[1] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s3 (d1_l2)
              scalar_type d1_s3  = PmA[1] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s2 (d1_l2)
              d1_s3             -= PmC[1] * (PmA[1] * F_mU[4] - PmC[1] * F_mU[5]); // p_s3 (d1_l2)
              d1_s0             += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d1_s1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d1_s2             += inv_two_zeta * (F_mU[2] - F_mU[3]);
              d1_s3             += inv_two_zeta * (F_mU[3] - F_mU[4]);
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[0] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l1)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                d1_p1_d2l1 += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                d1_p2_d2l1 += PmB[0] * d1_s2 - PmC[0] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[1] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l1)
                p_p0_d1l1_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l1_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l1_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1  = (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p1_d2l1  = (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p2_d2l1  = (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]) - (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p0_d2l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p1_d2l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p2_d2l1 += (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]) - (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                d1_p1_d2l1 += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                d1_p2_d2l1 += PmB[1] * d1_s2 - PmC[1] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1  = (PmB[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmB[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1  = p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[2] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[2] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l1)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[2] * d1_s0 - PmC[2] * d1_s1;
                d1_p1_d2l1 += PmB[2] * d1_s1 - PmC[2] * d1_s2;
                d1_p2_d2l1 += PmB[2] * d1_s2 - PmC[2] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2   = (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p1_d2l2   = (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                  d1_p1_d2l2  += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1  = p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[2] * p_p0_d1l1_d2l1 - PmC[2] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[2] * p_p1_d1l1_d2l1 - PmC[2] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[2] * d1_p0_d2l1 - PmC[2] * d1_p1_d2l1;
                  d_d1 += PmB[2] * d1_p1_d2l1 - PmC[2] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
            }
          }
          {
            {

              scalar_type d1_s0  = PmA[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
              d1_s0             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
              scalar_type d1_s1  = PmA[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
              d1_s1             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d1_l2)
              scalar_type d1_s2  = PmA[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d1_l2)
              d1_s2             -= PmC[2] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s3 (d1_l2)
              scalar_type d1_s3  = PmA[2] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s2 (d1_l2)
              d1_s3             -= PmC[2] * (PmA[0] * F_mU[4] - PmC[0] * F_mU[5]); // p_s3 (d1_l2)
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)
                scalar_type p_p0_d1l2_d2l1  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[0] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l2)
                p_p0_d1l2_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l2_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l2_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p1_d2l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p2_d2l1 += (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]) - (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                d1_p1_d2l1 += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                d1_p2_d2l1 += PmB[0] * d1_s2 - PmC[0] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l2  = (PmB[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmB[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmB[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)
                scalar_type p_p0_d1l2_d2l1  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[1] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l2)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                d1_p1_d2l1 += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                d1_p2_d2l1 += PmB[1] * d1_s2 - PmC[1] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2  = (PmB[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmB[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[1] * p_p0_d1l2_d2l1 - PmC[1] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[1] * p_p1_d1l2_d2l1 - PmC[1] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)
                p_p0_d1l1_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l1_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l1_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);
                scalar_type p_p0_d1l2_d2l1  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[2] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s1 (d1_l2)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1  = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p1_d2l1  = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p2_d2l1  = (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]) - (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[2] * d1_s0 - PmC[2] * d1_s1;
                d1_p1_d2l1 += PmB[2] * d1_s1 - PmC[2] * d1_s2;
                d1_p2_d2l1 += PmB[2] * d1_s2 - PmC[2] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                  d1_p1_d2l2  += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[1] * p_p0_d1l2_d2l1 - PmC[1] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[1] * p_p1_d1l2_d2l1 - PmC[1] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[2] * p_p0_d1l1_d2l1 - PmC[2] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[2] * p_p1_d1l1_d2l1 - PmC[2] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[2] * p_p0_d1l2_d2l1 - PmC[2] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[2] * p_p1_d1l2_d2l1 - PmC[2] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l2_d2l1 - p_p1_d1l2_d2l1;
                  d_d1  = p_p1_d1l2_d2l1 - p_p2_d1l2_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[2] * d1_p0_d2l1 - PmC[2] * d1_p1_d2l1;
                  d_d1 += PmB[2] * d1_p1_d2l1 - PmC[2] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
            }
            {

              scalar_type d1_s0  = PmA[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l2)
              d1_s0             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l2)
              scalar_type d1_s1  = PmA[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l2)
              d1_s1             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d1_l2)
              scalar_type d1_s2  = PmA[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d1_l2)
              d1_s2             -= PmC[2] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s3 (d1_l2)
              scalar_type d1_s3  = PmA[2] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s2 (d1_l2)
              d1_s3             -= PmC[2] * (PmA[1] * F_mU[4] - PmC[1] * F_mU[5]); // p_s3 (d1_l2)
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)
                scalar_type p_p0_d1l2_d2l1  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[0] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l2)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                d1_p1_d2l1 += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                d1_p2_d2l1 += PmB[0] * d1_s2 - PmC[0] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)
                scalar_type p_p0_d1l2_d2l1  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[1] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l2)
                p_p0_d1l2_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l2_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l2_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p1_d2l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p2_d2l1 += (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]) - (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                d1_p1_d2l1 += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                d1_p2_d2l1 += PmB[1] * d1_s2 - PmC[1] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l2  = (PmB[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmB[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[1] * p_p0_d1l2_d2l1 - PmC[1] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[1] * p_p1_d1l2_d2l1 - PmC[1] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)
                p_p0_d1l1_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l1_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l1_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);
                scalar_type p_p0_d1l2_d2l1  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d1_l2)
                scalar_type p_p1_d1l2_d2l1  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s1 (d1_l2)
                scalar_type p_p2_d1l2_d2l1  = PmB[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[2] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s1 (d1_l2)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1  = (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p1_d2l1  = (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p2_d2l1  = (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]) - (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[2] * d1_s0 - PmC[2] * d1_s1;
                d1_p1_d2l1 += PmB[2] * d1_s1 - PmC[2] * d1_s2;
                d1_p2_d2l1 += PmB[2] * d1_s2 - PmC[2] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[0] * p_p0_d1l2_d2l1 - PmC[0] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[0] * p_p1_d1l2_d2l1 - PmC[0] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = 1.0f;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2  += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2  *= inv_two_zeta;
                  d1_p1_d2l2  *= inv_two_zeta;
                  d1_p0_d2l2  += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                  d1_p1_d2l2  += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[1] * p_p0_d1l2_d2l1 - PmC[1] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[1] * p_p1_d1l2_d2l1 - PmC[1] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= !same_func * clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind += !same_func;
                  }

                  scalar_type p_d2_0_d1l2 = 0.0f, p_d2_1_d1l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l2 *= inv_two_zeta;
                  p_d2_1_d1l2 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[2] * p_p0_d1l1_d2l1 - PmC[2] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[2] * p_p1_d1l1_d2l1 - PmC[2] * p_p2_d1l1_d2l1;
                  p_d2_0_d1l2 += PmB[2] * p_p0_d1l2_d2l1 - PmC[2] * p_p1_d1l2_d2l1;
                  p_d2_1_d1l2 += PmB[2] * p_p1_d1l2_d2l1 - PmC[2] * p_p2_d1l2_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l2_d2l1 - p_p1_d1l2_d2l1;
                  d_d1  = p_p1_d1l2_d2l1 - p_p2_d1l2_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[2] * d1_p0_d2l1 - PmC[2] * d1_p1_d2l1;
                  d_d1 += PmB[2] * d1_p1_d2l1 - PmC[2] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l2;
                    AB_common     = p_d2_0_d1l2;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l2;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
            }
            {

              scalar_type d1_s0  = PmA[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l2)
              d1_s0             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l2)
              scalar_type d1_s1  = PmA[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l2)
              d1_s1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d1_l2)
              scalar_type d1_s2  = PmA[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d1_l2)
              d1_s2             -= PmC[2] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s3 (d1_l2)
              scalar_type d1_s3  = PmA[2] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s2 (d1_l2)
              d1_s3             -= PmC[2] * (PmA[2] * F_mU[4] - PmC[2] * F_mU[5]); // p_s3 (d1_l2)
              d1_s0             += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d1_s1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d1_s2             += inv_two_zeta * (F_mU[2] - F_mU[3]);
              d1_s3             += inv_two_zeta * (F_mU[3] - F_mU[4]);
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[0] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                d1_p1_d2l1 += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                d1_p2_d2l1 += PmB[0] * d1_s2 - PmC[0] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[1] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1 += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                d1_p1_d2l1 += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                d1_p2_d2l1 += PmB[1] * d1_s2 - PmC[1] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
              {

                scalar_type p_p0_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d1_l1)
                scalar_type p_p1_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s1 (d1_l1)
                scalar_type p_p2_d1l1_d2l1  = PmB[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[2] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s1 (d1_l1)
                p_p0_d1l1_d2l1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1l1_d2l1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                p_p2_d1l1_d2l1             += inv_two_zeta * (F_mU[2] - F_mU[3]);

                scalar_type d1_p0_d2l1 = 0.0f, d1_p1_d2l1 = 0.0f, d1_p2_d2l1 = 0.0f;
                d1_p0_d2l1  = (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p1_d2l1  = (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p2_d2l1  = (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]) - (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p0_d2l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p1_d2l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p2_d2l1 += (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]) - (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p0_d2l1 *= inv_two_zeta;
                d1_p1_d2l1 *= inv_two_zeta;
                d1_p2_d2l1 *= inv_two_zeta;
                d1_p0_d2l1 += PmB[2] * d1_s0 - PmC[2] * d1_s1;
                d1_p1_d2l1 += PmB[2] * d1_s1 - PmC[2] * d1_s2;
                d1_p2_d2l1 += PmB[2] * d1_s2 - PmC[2] * d1_s3;
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[0] * d1_s0 - PmC[0] * d1_s1;
                  d1_p1_d2l2  += PmB[0] * d1_s1 - PmC[0] * d1_s2;
                  p_d2_0_d1l1 += PmB[0] * p_p0_d1l1_d2l1 - PmC[0] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[0] * p_p1_d1l1_d2l1 - PmC[0] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[0] * d1_p0_d2l1 - PmC[0] * d1_p1_d2l1;
                  d_d1 += PmB[0] * d1_p1_d2l1 - PmC[0] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[0] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type d1_p0_d2l2 = 0.0f, d1_p1_d2l2 = 0.0f;
                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  d1_p0_d2l2  += PmB[1] * d1_s0 - PmC[1] * d1_s1;
                  d1_p1_d2l2  += PmB[1] * d1_s1 - PmC[1] * d1_s2;
                  p_d2_0_d1l1 += PmB[1] * p_p0_d1l1_d2l1 - PmC[1] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[1] * p_p1_d1l1_d2l1 - PmC[1] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0 += PmB[1] * d1_p0_d2l1 - PmC[1] * d1_p1_d2l1;
                  d_d1 += PmB[1] * d1_p1_d2l1 - PmC[1] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[1] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l2;
                    AB_common    += d1_p0_d2l2;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l2;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
                {
                  scalar_type pre_term;
                  {
                    pre_term  = gpu_normalization_factor * gpu_normalization_factor;
                    pre_term *= clatom_charge_sh[j] * dens[dens_ind];
                    dens_ind++;
                  }

                  scalar_type p_d2_0_d1l1 = 0.0f, p_d2_1_d1l1 = 0.0f;
                  p_d2_0_d1l1  = (PmB[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1  = (PmB[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmB[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1 *= inv_two_zeta;
                  p_d2_1_d1l1 *= inv_two_zeta;
                  p_d2_0_d1l1 += PmB[2] * p_p0_d1l1_d2l1 - PmC[2] * p_p1_d1l1_d2l1;
                  p_d2_1_d1l1 += PmB[2] * p_p1_d1l1_d2l1 - PmC[2] * p_p2_d1l1_d2l1;

                  scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                  d_d0  = p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1  = p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1;
                  d_d1 += p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1;
                  d_d0 += d1_s0 - d1_s1;
                  d_d1 += d1_s1 - d1_s2;
                  d_d0 *= inv_two_zeta;
                  d_d1 *= inv_two_zeta;
                  d_d0 += PmB[2] * d1_p0_d2l1 - PmC[2] * d1_p1_d2l1;
                  d_d1 += PmB[2] * d1_p1_d2l1 - PmC[2] * d1_p2_d2l1;
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[0] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[0] * d_d0 + A_force_term;
                    A_force_term += PmA[0] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[0]     += pre_term * A_force_term;
                    B_force[0]     += pre_term * B_force_term;
                    C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = PmC[1] * d_d1;

                    A_force_term  = -C_force_term;
                    B_force_term  = PmB[1] * d_d0 + A_force_term;
                    A_force_term += PmA[1] * d_d0;
                    A_force_term *= 2.0f * ai;
                    B_force_term *= 2.0f * aj;
                    A_force[1]     += pre_term * A_force_term;
                    B_force[1]     += pre_term * B_force_term;
                    C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                  {
                    C_force_term  = 0.0f;
                    AB_common     = 0.0f;
                    C_force_term  = p_d2_1_d1l1;
                    AB_common     = p_d2_0_d1l1;
                    C_force_term += p_d2_1_d1l1;
                    AB_common    += p_d2_0_d1l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term += d1_p1_d2l1;
                    AB_common    += d1_p0_d2l1;
                    C_force_term  = PmC[2] * d_d1 + inv_two_zeta * C_force_term; 

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[2] * d_d0 + A_force_term;
                    A_force_term += PmA[2] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= p_d2_0_d1l1;
                    A_force_term -= p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= d1_p0_d2l1;
                    B_force_term -= d1_p0_d2l1;
                    A_force[2]     += pre_term * A_force_term;
                    B_force[2]     += pre_term * B_force_term;
                    C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-D)-------------------------------------------
