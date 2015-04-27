//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (PS-D)-------------------------------------------
              {
                scalar_type F_mT[5];
                {
                  scalar_type PmQ[3];
                  PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
                  PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
                  PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
                  scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
                  lio_gamma<scalar_type,4>(F_mT,T);
                }
                {
                  scalar_type A_force_term, B_force_term, C_force_term;
                  #pragma unroll 3
                  for (uint p_l = 0; p_l < 3; p_l++) {

                    scalar_type ps_s0 = PmA[p_l] * F_mT[0] + WmP[p_l] * F_mT[1];
                    scalar_type ps_s1 = PmA[p_l] * F_mT[1] + WmP[p_l] * F_mT[2];
                    scalar_type ps_s2 = PmA[p_l] * F_mT[2] + WmP[p_l] * F_mT[3];
                    scalar_type ps_s3 = PmA[p_l] * F_mT[3] + WmP[p_l] * F_mT[4];

                    uint fit_dens_ind = 0;
                    //#pragma unroll 3
                    for (uint d_l1 = 0; d_l1 < 3; d_l1++) {

                      scalar_type ss_p1_1 = WmQ[d_l1] * F_mT[2];
                      scalar_type ss_p1_2 = WmQ[d_l1] * F_mT[3];

                      scalar_type ps_p1_0 = WmQ[d_l1] * ps_s1;
                      scalar_type ps_p1_1 = WmQ[d_l1] * ps_s2;
                      scalar_type ps_p1_2 = WmQ[d_l1] * ps_s3;
                      {
                        bool del_p1p3 = p_l == d_l1;
                        ps_p1_1 += del_p1p3 * inv_two_zeta_eta * F_mT[2];
                        ps_p1_2 += del_p1p3 * inv_two_zeta_eta * F_mT[3];
                      }

                      for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {

                        scalar_type ps_p2_0 = WmQ[d_l2] * ps_s1;
                        scalar_type ps_p2_1 = WmQ[d_l2] * ps_s2;

                        scalar_type ss_d0   = WmQ[d_l2] * ss_p1_1;
                        scalar_type ss_d1   = WmQ[d_l2] * ss_p1_2;

                        scalar_type ps_d0   = WmQ[d_l2] * ps_p1_1;
                        scalar_type ps_d1   = WmQ[d_l2] * ps_p1_2;

                        {
                          bool del_13 = p_l == d_l2;
                          ps_d0   += del_13 * inv_two_zeta_eta * ss_p1_1;
                          ps_d1   += del_13 * inv_two_zeta_eta * ss_p1_2;
                          ps_p2_0 += del_13 * inv_two_zeta_eta * F_mT[1];
                          ps_p2_1 += del_13 * inv_two_zeta_eta * F_mT[2];
                        }

                        scalar_type norm;
                        {
                          bool del_23 = d_l1 == d_l2;
                          ss_d0 += del_23 * inv_two_ak * (F_mT[0] - rho_ak * F_mT[1]);
                          ss_d1 += del_23 * inv_two_ak * (F_mT[1] - rho_ak * F_mT[2]);
                          ps_d0 += del_23 * inv_two_ak * (ps_s0 - rho_ak * ps_s1);
                          ps_d1 += del_23 * inv_two_ak * (ps_s1 - rho_ak * ps_s2);

                          norm = del_23 * G2G::gpu_normalization_factor + !del_23 * 1.0f;
                        }

                        scalar_type pre_term = norm * prefactor_dens * dens[p_l] * fit_dens_sh[j+fit_dens_ind];
                        fit_dens_ind++;

                        #pragma unroll 3
                        for (uint grad_l = 0; grad_l < 3; grad_l++) {
                          bool del_1g = p_l == grad_l, del_2g = d_l1 == grad_l, del_3g = d_l2 == grad_l;

                          A_force_term  = WmP[grad_l] * ps_d1;
                          A_force_term += del_1g * inv_two_zeta * (ss_d0 - rho_zeta * ss_d1);
                          A_force_term += del_2g * inv_two_zeta_eta * ps_p2_1;
                          A_force_term += del_3g * inv_two_zeta_eta * ps_p1_1;
                          B_force_term  = PmB[grad_l] * ps_d0 + A_force_term;
                          A_force_term  = PmA[grad_l] * ps_d0 + A_force_term;

                          C_force_term  = WmQ[grad_l] * ps_d1;
                          C_force_term += del_1g * inv_two_zeta_eta * ss_d1;
                          C_force_term += del_2g * inv_two_ak * (ps_p2_0 - rho_ak * ps_p2_1);
                          C_force_term += del_3g * inv_two_ak * (ps_p1_0 - rho_ak * ps_p1_1);

                          A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_1g * ss_d0);
                          B_force[grad_l]      += pre_term * 2.0f * aj * B_force_term;
                          C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_2g * ps_p2_0 - del_3g * ps_p1_0);
                        }
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo;
                C_force[1][tid] *= valid_thread * prefactor_mo;
                C_force[2][tid] *= valid_thread * prefactor_mo;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (PS-D)----------------------------------------------
