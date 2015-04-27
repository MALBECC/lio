//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-D)-------------------------------------------
        scalar_type F_mU[6];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * zeta;
          //F_mU[0] = (SQRT_PI / (2*sqrtU)) * erff(sqrtU);
          for (int m = 0; m <= 5; m++) 
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
          scalar_type AB_common;
          //scalar_type mm_charge = clatom_charge_sh[j];
          uint dens_ind = 0;

          for (int d1_l1 = 0; d1_l1 < 3; d1_l1++)
          {
            for (int d1_l2 = 0; d1_l2 <= d1_l1; d1_l2++)
            {
              bool del_d1 = d1_l1 == d1_l2;

              scalar_type d1_s0  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]); // p_s0 (d1_l2)
              d1_s0             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s1 (d1_l2)
              d1_s0             += del_d1 * inv_two_zeta * (F_mU[0] - F_mU[1]);
              scalar_type d1_s1  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s1 (d1_l2)
              d1_s1             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s2 (d1_l2)
              d1_s1             += del_d1 * inv_two_zeta * (F_mU[1] - F_mU[2]);
              scalar_type d1_s2  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s2 (d1_l2)
              d1_s2             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]); // p_s3 (d1_l2)
              d1_s2             += del_d1 * inv_two_zeta * (F_mU[2] - F_mU[3]);
              scalar_type d1_s3  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]); // p_s2 (d1_l2)
              d1_s3             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[4] - PmC[d1_l2] * F_mU[5]); // p_s3 (d1_l2)
              d1_s3             += del_d1 * inv_two_zeta * (F_mU[3] - F_mU[4]);
              for (int d2_l1 = 0; d2_l1 < 3; d2_l1++)
              {
                bool del_d1l1_d2l1 = d1_l1 == d2_l1, del_d1l2_d2l1 = d1_l2 == d2_l1;

                scalar_type p_p0_d1l1_d2l1  = PmB[d2_l1] * (PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]); // p_s0 (d1_l1)
                p_p0_d1l1_d2l1             -= PmC[d2_l1] * (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]); // p_s1 (d1_l1)
                p_p0_d1l1_d2l1             += del_d1l1_d2l1 * inv_two_zeta * (F_mU[0] - F_mU[1]);
                scalar_type p_p1_d1l1_d2l1  = PmB[d2_l1] * (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]); // p_s0 (d1_l1)
                p_p1_d1l1_d2l1             -= PmC[d2_l1] * (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]); // p_s1 (d1_l1)
                p_p1_d1l1_d2l1             += del_d1l1_d2l1 * inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p2_d1l1_d2l1  = PmB[d2_l1] * (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]); // p_s0 (d1_l1)
                p_p2_d1l1_d2l1             -= PmC[d2_l1] * (PmA[d1_l1] * F_mU[3] - PmC[d1_l1] * F_mU[4]); // p_s1 (d1_l1)
                p_p2_d1l1_d2l1             += del_d1l1_d2l1 * inv_two_zeta * (F_mU[2] - F_mU[3]);
                scalar_type p_p0_d1l2_d2l1  = PmB[d2_l1] * (PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]); // p_s0 (d1_l2)
                p_p0_d1l2_d2l1             -= PmC[d2_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s1 (d1_l2)
                p_p0_d1l2_d2l1             += del_d1l2_d2l1 * inv_two_zeta * (F_mU[0] - F_mU[1]);
                scalar_type p_p1_d1l2_d2l1  = PmB[d2_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s0 (d1_l2)
                p_p1_d1l2_d2l1             -= PmC[d2_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s1 (d1_l2)
                p_p1_d1l2_d2l1             += del_d1l2_d2l1 * inv_two_zeta * (F_mU[1] - F_mU[2]);
                scalar_type p_p2_d1l2_d2l1  = PmB[d2_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s0 (d1_l2)
                p_p2_d1l2_d2l1             -= PmC[d2_l1] * (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]); // p_s1 (d1_l2)
                p_p2_d1l2_d2l1             += del_d1l2_d2l1 * inv_two_zeta * (F_mU[2] - F_mU[3]);

                scalar_type d1_p0_d2l1       = del_d1l1_d2l1 * ((PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]) - (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p0_d2l1                  += del_d1l2_d2l1 * ((PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]) - (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p0_d2l1                  *= inv_two_zeta;
                d1_p0_d2l1                  += PmB[d2_l1] * d1_s0 - PmC[d2_l1] * d1_s1;
                scalar_type d1_p1_d2l1       = del_d1l1_d2l1 * ((PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]) - (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p1_d2l1                  += del_d1l2_d2l1 * ((PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]) - (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p1_d2l1                  *= inv_two_zeta;
                d1_p1_d2l1                  += PmB[d2_l1] * d1_s1 - PmC[d2_l1] * d1_s2;
                scalar_type d1_p2_d2l1       = del_d1l1_d2l1 * ((PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]) - (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                d1_p2_d2l1                  += del_d1l2_d2l1 * ((PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]) - (PmA[d1_l1] * F_mU[3] - PmC[d1_l1] * F_mU[4]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                d1_p2_d2l1                  *= inv_two_zeta;
                d1_p2_d2l1                  += PmB[d2_l1] * d1_s2 - PmC[d2_l1] * d1_s3;
                for (int d2_l2 = 0; d2_l2 <= d2_l1; d2_l2++)
                {
                  bool del_d2 = d2_l1 == d2_l2, del_d1l1_d2l2 = d1_l1 == d2_l2, del_d1l2_d2l2 = d1_l2 == d2_l2;
                  scalar_type pre_term;
                  {
                    bool skip = same_func && (d2_l1 > d1_l1 || (d2_l1 == d1_l1 && d2_l2 > d1_l2));
                    pre_term  = !skip * clatom_charge_sh[j] * dens[dens_ind];
                    //if (valid_thread) printf("%.4e\n",dens[dens_ind]);
                    pre_term *= (!del_d1&&!del_d2)*1.0f + ((del_d1&&!del_d2)||(del_d2&&!del_d1))*G2G::gpu_normalization_factor + (del_d1&&del_d2)*G2G::gpu_normalization_factor*G2G::gpu_normalization_factor;
                    dens_ind += !skip;
                  }

                  scalar_type d1_p0_d2l2      = del_d1l1_d2l2 * ((PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]) - (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p0_d2l2                 += del_d1l2_d2l2 * ((PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]) - (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p0_d2l2                 *= inv_two_zeta;
                  d1_p0_d2l2                 += PmB[d2_l2] * d1_s0 - PmC[d2_l2] * d1_s1;
                  scalar_type d1_p1_d2l2      = del_d1l1_d2l2 * ((PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]) - (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  d1_p1_d2l2                 += del_d1l2_d2l2 * ((PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]) - (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  d1_p1_d2l2                 *= inv_two_zeta;
                  d1_p1_d2l2                 += PmB[d2_l2] * d1_s1 - PmC[d2_l2] * d1_s2;
                  scalar_type p_d2_0_d1l1     = del_d1l1_d2l2 * ((PmB[d2_l1] * F_mU[0] - PmC[d2_l1] * F_mU[1]) - (PmB[d2_l1] * F_mU[1] - PmC[d2_l1] * F_mU[2]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l1                += del_d2 * ((PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]) - (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l1                *= inv_two_zeta;
                  p_d2_0_d1l1                += PmB[d2_l2] * p_p0_d1l1_d2l1 - PmC[d2_l2] * p_p1_d1l1_d2l1;
                  scalar_type p_d2_1_d1l1     = del_d1l1_d2l2 * ((PmB[d2_l1] * F_mU[1] - PmC[d2_l1] * F_mU[2]) - (PmB[d2_l1] * F_mU[2] - PmC[d2_l1] * F_mU[3]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l1                += del_d2 * ((PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]) - (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l1                *= inv_two_zeta;
                  p_d2_1_d1l1                += PmB[d2_l2] * p_p1_d1l1_d2l1 - PmC[d2_l2] * p_p2_d1l1_d2l1;
                  scalar_type p_d2_0_d1l2     = del_d1l2_d2l2 * ((PmB[d2_l1] * F_mU[0] - PmC[d2_l1] * F_mU[1]) - (PmB[d2_l1] * F_mU[1] - PmC[d2_l1] * F_mU[2]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_0_d1l2                += del_d2 * ((PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]) - (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_0_d1l2                *= inv_two_zeta;
                  p_d2_0_d1l2                += PmB[d2_l2] * p_p0_d1l2_d2l1 - PmC[d2_l2] * p_p1_d1l2_d2l1;
                  scalar_type p_d2_1_d1l2     = del_d1l2_d2l2 * ((PmB[d2_l1] * F_mU[1] - PmC[d2_l1] * F_mU[2]) - (PmB[d2_l1] * F_mU[2] - PmC[d2_l1] * F_mU[3]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
                  p_d2_1_d1l2                += del_d2 * ((PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]) - (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
                  p_d2_1_d1l2                *= inv_two_zeta;
                  p_d2_1_d1l2                += PmB[d2_l2] * p_p1_d1l2_d2l1 - PmC[d2_l2] * p_p2_d1l2_d2l1;

                  scalar_type d_d0            = del_d1l1_d2l2 * (p_p0_d1l2_d2l1 - p_p1_d1l2_d2l1);
                  d_d0                       += del_d1l2_d2l2 * (p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1);
                  d_d0                       += del_d2 * (d1_s0 - d1_s1);
                  d_d0                       *= inv_two_zeta;
                  d_d0                       += PmB[d2_l2] * d1_p0_d2l1 - PmC[d2_l2] * d1_p1_d2l1;
                  scalar_type d_d1            = del_d1l1_d2l2 * (p_p1_d1l2_d2l1 - p_p2_d1l2_d2l1);
                  d_d1                       += del_d1l2_d2l2 * (p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1);
                  d_d1                       += del_d2 * (d1_s1 - d1_s2);
                  d_d1                       *= inv_two_zeta;
                  d_d1                       += PmB[d2_l2] * d1_p1_d2l1 - PmC[d2_l2] * d1_p2_d2l1;
                  for (int grad_l = 0; grad_l < 3; grad_l++)
                  {
                    bool del_d1l1g = d1_l1 == grad_l, del_d1l2g = d1_l2 == grad_l, del_d2l1g = d2_l1 == grad_l, del_d2l2g = d2_l2 == grad_l;
                    C_force_term  = del_d1l1g * p_d2_1_d1l2;
                    C_force_term += del_d1l2g * p_d2_1_d1l1;
                    C_force_term += del_d2l1g * d1_p1_d2l2;
                    C_force_term += del_d2l2g * d1_p1_d2l1;
                    C_force_term  = PmC[grad_l] * d_d1 + inv_two_zeta * C_force_term;
  
                    AB_common     = del_d1l1g * p_d2_0_d1l2;
                    AB_common    += del_d1l2g * p_d2_0_d1l1;
                    AB_common    += del_d2l1g * d1_p0_d2l2;
                    AB_common    += del_d2l2g * d1_p0_d2l1;

                    A_force_term  = inv_two_zeta * AB_common - C_force_term;
                    B_force_term  = PmB[grad_l] * d_d0 + A_force_term;
                    A_force_term += PmA[grad_l] * d_d0;
                    A_force_term *= 2.0f * ai;
                    A_force_term -= del_d1l1g * p_d2_0_d1l2;
                    A_force_term -= del_d1l2g * p_d2_0_d1l1;
                    B_force_term *= 2.0f * aj;
                    B_force_term -= del_d2l1g * d1_p0_d2l2;
                    B_force_term -= del_d2l2g * d1_p0_d2l1;

                    A_force[grad_l]     += pre_term * A_force_term;
                    B_force[grad_l]     += pre_term * B_force_term;
                    C_force[grad_l][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
                  }
                }
              }
            }
          }
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-D)-------------------------------------------
