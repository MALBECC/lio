// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-P)-------------------------------------------
        scalar_type F_mU[5];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          lio_gamma<scalar_type,4>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common;
          uint dens_ind = 0;

          // d_l1 = 0
          {
            // d_l2 = 0
            {

              scalar_type d_s0  = PmA[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
              d_s0             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              scalar_type d_s1  = PmA[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              scalar_type d_s2  = PmA[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              d_s2             -= PmC[0] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s3 (d_l2)
              d_s0             += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d_s2             += inv_two_zeta * (F_mU[2] - F_mU[3]);
              // p_l = 0
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l1)
                p_p0_d1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p0_d2  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                p_p0_d2             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d2             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0  = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d_l2) - p_s1 (d_l2)
                d_p0 += (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d_l1) - p_s1 (d_l1)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[0] * d_s0 - PmC[0] * d_s1;
                d_p1  = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s1 (d_l2) - p_s2 (d_l2)
                d_p1 += (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s1 (d_l1) - p_s2 (d_l1)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[0] * d_s1 - PmC[0] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 1
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[1] * d_s0 - PmC[1] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[1] * d_s1 - PmC[1] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 2
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[2] * d_s0 - PmC[2] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[2] * d_s1 - PmC[2] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
            }
          }
          // d_l1 = 1
          {
            // d_l2 = 0
            {

              scalar_type d_s0  = PmA[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
              d_s0             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              scalar_type d_s1  = PmA[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              scalar_type d_s2  = PmA[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              d_s2             -= PmC[1] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s3 (d_l2)
              // p_l = 0
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                p_p0_d2             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d2             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d_l1) - p_s1 (d_l1)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[0] * d_s0 - PmC[0] * d_s1;
                d_p1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s1 (d_l1) - p_s2 (d_l1)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[0] * d_s1 - PmC[0] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 1
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l1)
                p_p0_d1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p0_d2  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0  = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d_l2) - p_s1 (d_l2)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[1] * d_s0 - PmC[1] * d_s1;
                d_p1  = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s1 (d_l2) - p_s2 (d_l2)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[1] * d_s1 - PmC[1] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 2
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[2] * d_s0 - PmC[2] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[2] * d_s1 - PmC[2] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
            }
            // d_l2 = 1
            {

              scalar_type d_s0  = PmA[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
              d_s0             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              scalar_type d_s1  = PmA[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              d_s1             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
              scalar_type d_s2  = PmA[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
              d_s2             -= PmC[1] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s3 (d_l2)
              d_s0             += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d_s2             += inv_two_zeta * (F_mU[2] - F_mU[3]);
              // p_l = 0
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[0] * d_s0 - PmC[0] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[0] * d_s1 - PmC[0] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 1
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l1)
                p_p0_d1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p0_d2  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
                p_p0_d2             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d2             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0  = (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d_l2) - p_s1 (d_l2)
                d_p0 += (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d_l1) - p_s1 (d_l1)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[1] * d_s0 - PmC[1] * d_s1;
                d_p1  = (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s1 (d_l2) - p_s2 (d_l2)
                d_p1 += (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s1 (d_l1) - p_s2 (d_l1)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[1] * d_s1 - PmC[1] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 2
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[2] * d_s0 - PmC[2] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[2] * d_s1 - PmC[2] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
            }
          }
          // d_l1 = 2
          {
            // d_l2 = 0
            {

              scalar_type d_s0  = PmA[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
              d_s0             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              scalar_type d_s1  = PmA[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
              d_s1             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              scalar_type d_s2  = PmA[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
              d_s2             -= PmC[2] * (PmA[0] * F_mU[3] - PmC[0] * F_mU[4]); // p_s3 (d_l2)
              // p_l = 0
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[0] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[0] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[0] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                p_p0_d2             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d2             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d_l1) - p_s1 (d_l1)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[0] * d_s0 - PmC[0] * d_s1;
                d_p1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s1 (d_l1) - p_s2 (d_l1)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[0] * d_s1 - PmC[0] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 1
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[1] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[1] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[1] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[1] * d_s0 - PmC[1] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[1] * d_s1 - PmC[1] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 2
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                p_p0_d1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p0_d2  = PmB[2] * (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[2] * (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[2] * (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0  = (PmA[0] * F_mU[0] - PmC[0] * F_mU[1]) - (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]);  // p_s0 (d_l2) - p_s1 (d_l2)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[2] * d_s0 - PmC[2] * d_s1;
                d_p1  = (PmA[0] * F_mU[1] - PmC[0] * F_mU[2]) - (PmA[0] * F_mU[2] - PmC[0] * F_mU[3]);  // p_s1 (d_l2) - p_s2 (d_l2)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[2] * d_s1 - PmC[2] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
            }
            // d_l2 = 1
            {

              scalar_type d_s0  = PmA[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
              d_s0             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              scalar_type d_s1  = PmA[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
              d_s1             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
              scalar_type d_s2  = PmA[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
              d_s2             -= PmC[2] * (PmA[1] * F_mU[3] - PmC[1] * F_mU[4]); // p_s3 (d_l2)
              // p_l = 0
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[0] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[0] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[0] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[0] * d_s0 - PmC[0] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[0] * d_s1 - PmC[0] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 1
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[1] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[1] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[1] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
                p_p0_d2             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d2             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d_l1) - p_s1 (d_l1)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[1] * d_s0 - PmC[1] * d_s1;
                d_p1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s1 (d_l1) - p_s2 (d_l1)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[1] * d_s1 - PmC[1] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 2
              {
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                p_p0_d1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p0_d2  = PmB[2] * (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[2] * (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[2] * (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0  = (PmA[1] * F_mU[0] - PmC[1] * F_mU[1]) - (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]);  // p_s0 (d_l2) - p_s1 (d_l2)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[2] * d_s0 - PmC[2] * d_s1;
                d_p1  = (PmA[1] * F_mU[1] - PmC[1] * F_mU[2]) - (PmA[1] * F_mU[2] - PmC[1] * F_mU[3]);  // p_s1 (d_l2) - p_s2 (d_l2)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[2] * d_s1 - PmC[2] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
            }
            // d_l2 = 2
            {

              scalar_type d_s0  = PmA[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l2)
              d_s0             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
              scalar_type d_s1  = PmA[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
              d_s1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l2)
              scalar_type d_s2  = PmA[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l2)
              d_s2             -= PmC[2] * (PmA[2] * F_mU[3] - PmC[2] * F_mU[4]); // p_s3 (d_l2)
              d_s0             += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d_s2             += inv_two_zeta * (F_mU[2] - F_mU[3]);
              // p_l = 0
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[0] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[0] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[0] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[0] * d_s0 - PmC[0] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[0] * d_s1 - PmC[0] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 1
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                scalar_type p_p0_d2  = PmB[1] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[1] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[1] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l2)
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[1] * d_s0 - PmC[1] * d_s1;
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[1] * d_s1 - PmC[1] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
              // p_l = 2
              {
                scalar_type pre_term = G2G::gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;

                scalar_type p_p0_d1  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l1)
                p_p0_d1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p0_d2  = PmB[2] * (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[2] * (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[2] * (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]); // p_s2 (d_l2)
                p_p0_d2             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d2             += inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
                d_p0  = (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d_l2) - p_s1 (d_l2)
                d_p0 += (PmA[2] * F_mU[0] - PmC[2] * F_mU[1]) - (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]);  // p_s0 (d_l1) - p_s1 (d_l1)
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[2] * d_s0 - PmC[2] * d_s1;
                d_p1  = (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s1 (d_l2) - p_s2 (d_l2)
                d_p1 += (PmA[2] * F_mU[1] - PmC[2] * F_mU[2]) - (PmA[2] * F_mU[2] - PmC[2] * F_mU[3]);  // p_s1 (d_l1) - p_s2 (d_l1)
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[2] * d_s1 - PmC[2] * d_s2;
                // grad_l = 0
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[0] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[0] * d_p0 + A_force_term;
                  A_force_term += PmA[0] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[0]     += pre_term * A_force_term;
                  B_force[0]     += pre_term * B_force_term;
                  C_force[0][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 1
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term  = PmC[1] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[1] * d_p0 + A_force_term;
                  A_force_term += PmA[1] * d_p0;
                  A_force_term *= 2.0f * ai;
                  B_force_term *= 2.0f * aj;
                  A_force[1]     += pre_term * A_force_term;
                  B_force[1]     += pre_term * B_force_term;
                  C_force[1][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
                // grad_l = 2
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
                  C_force_term += d_s1;
                  AB_common    += d_s0;
                  C_force_term  = PmC[2] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[2] * d_p0 + A_force_term;
                  A_force_term += PmA[2] * d_p0;
                  A_force_term *= 2.0f * ai;
                  A_force_term -= p_p0_d2;
                  A_force_term -= p_p0_d1;
                  B_force_term *= 2.0f * aj;
                  B_force_term -= d_s0;
                  A_force[2]     += pre_term * A_force_term;
                  B_force[2]     += pre_term * B_force_term;
                  C_force[2][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                }
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-P)-------------------------------------------
