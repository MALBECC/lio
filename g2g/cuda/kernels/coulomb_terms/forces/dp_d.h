//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DP-D)-------------------------------------------
              {
                scalar_type F_mT[7];
                {
                  scalar_type PmQ[3];
                  PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
                  PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
                  PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
                  scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
                  lio_gamma<scalar_type,6>(F_mT,T);
                }
                {
                  scalar_type A_force_term, B_force_term, C_force_term;
                  uint mo_dens_ind = 0;
                  //#pragma unroll 3
                  for (uint d1_l1 = 0; d1_l1 < 3; d1_l1++) {

                    scalar_type p1s_s0 = PmA[d1_l1] * F_mT[0] + WmP[d1_l1] * F_mT[1];
                    scalar_type p1s_s1 = PmA[d1_l1] * F_mT[1] + WmP[d1_l1] * F_mT[2];
                    scalar_type p1s_s2 = PmA[d1_l1] * F_mT[2] + WmP[d1_l1] * F_mT[3];
                    scalar_type p1s_s3 = PmA[d1_l1] * F_mT[3] + WmP[d1_l1] * F_mT[4];
                    scalar_type p1s_s4 = PmA[d1_l1] * F_mT[4] + WmP[d1_l1] * F_mT[5];
                    scalar_type p1s_s5 = PmA[d1_l1] * F_mT[5] + WmP[d1_l1] * F_mT[6];

                    //#pragma unroll 3
                    for (uint d1_l2 = 0; d1_l2 <= d1_l1; d1_l2++) {

                      scalar_type p2s_s0 = PmA[d1_l2] * F_mT[0] + WmP[d1_l2] * F_mT[1];
                      scalar_type p2s_s1 = PmA[d1_l2] * F_mT[1] + WmP[d1_l2] * F_mT[2];
                      scalar_type p2s_s2 = PmA[d1_l2] * F_mT[2] + WmP[d1_l2] * F_mT[3];
                      scalar_type p2s_s3 = PmA[d1_l2] * F_mT[3] + WmP[d1_l2] * F_mT[4];
                      scalar_type p2s_s4 = PmA[d1_l2] * F_mT[4] + WmP[d1_l2] * F_mT[5];

                      scalar_type ds_s0 = PmA[d1_l2] * p1s_s0 + WmP[d1_l2] * p1s_s1;
                      scalar_type ds_s1 = PmA[d1_l2] * p1s_s1 + WmP[d1_l2] * p1s_s2;
                      scalar_type ds_s2 = PmA[d1_l2] * p1s_s2 + WmP[d1_l2] * p1s_s3;
                      scalar_type ds_s3 = PmA[d1_l2] * p1s_s3 + WmP[d1_l2] * p1s_s4;
                      scalar_type ds_s4 = PmA[d1_l2] * p1s_s4 + WmP[d1_l2] * p1s_s5;

                      scalar_type norm;
                      {
                        bool del_12 = d1_l1 == d1_l2;
                        ds_s0 += del_12 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        ds_s1 += del_12 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        ds_s2 += del_12 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                        ds_s3 += del_12 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                        ds_s4 += del_12 * inv_two_zeta * (F_mT[4] - rho_zeta * F_mT[5]);
                        norm = del_12 * gpu_normalization_factor + !del_12 * 1.0f;
                      }

                      #pragma unroll 3
                      for (uint p_l = 0; p_l < 3; p_l++) {

                        scalar_type sp_s2  = PmB[p_l] * F_mT[2] + WmP[p_l] * F_mT[3];
                        scalar_type sp_s3  = PmB[p_l] * F_mT[3] + WmP[p_l] * F_mT[4];

                        scalar_type p1p_s0 = PmB[p_l] * p1s_s0 + WmP[p_l] * p1s_s1;
                        scalar_type p1p_s1 = PmB[p_l] * p1s_s1 + WmP[p_l] * p1s_s2;
                        scalar_type p1p_s2 = PmB[p_l] * p1s_s2 + WmP[p_l] * p1s_s3;
                        scalar_type p1p_s3 = PmB[p_l] * p1s_s3 + WmP[p_l] * p1s_s4;

                        scalar_type p2p_s0 = PmB[p_l] * p2s_s0 + WmP[p_l] * p2s_s1;
                        scalar_type p2p_s1 = PmB[p_l] * p2s_s1 + WmP[p_l] * p2s_s2;
                        scalar_type p2p_s2 = PmB[p_l] * p2s_s2 + WmP[p_l] * p2s_s3;
                        scalar_type p2p_s3 = PmB[p_l] * p2s_s3 + WmP[p_l] * p2s_s4;

                        scalar_type dp_s0  = PmB[p_l] * ds_s0 + WmP[p_l] * ds_s1;
                        scalar_type dp_s1  = PmB[p_l] * ds_s1 + WmP[p_l] * ds_s2;
                        scalar_type dp_s2  = PmB[p_l] * ds_s2 + WmP[p_l] * ds_s3;
                        scalar_type dp_s3  = PmB[p_l] * ds_s3 + WmP[p_l] * ds_s4;

                        {
                          bool del_13 = d1_l1 == p_l;
                          p1p_s0 += del_13 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p1p_s1 += del_13 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p1p_s2 += del_13 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          p1p_s3 += del_13 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                          dp_s0  += del_13 * inv_two_zeta * (p2s_s0 - rho_zeta * p2s_s1);
                          dp_s1  += del_13 * inv_two_zeta * (p2s_s1 - rho_zeta * p2s_s2);
                          dp_s2  += del_13 * inv_two_zeta * (p2s_s2 - rho_zeta * p2s_s3);
                          dp_s3  += del_13 * inv_two_zeta * (p2s_s3 - rho_zeta * p2s_s4);
                        }
                        {
                          bool del_23 = d1_l2 == p_l;
                          p2p_s0 += del_23 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p2p_s1 += del_23 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p2p_s2 += del_23 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          p2p_s3 += del_23 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                          dp_s0  += del_23 * inv_two_zeta * (p1s_s0 - rho_zeta * p1s_s1);
                          dp_s1  += del_23 * inv_two_zeta * (p1s_s1 - rho_zeta * p1s_s2);
                          dp_s2  += del_23 * inv_two_zeta * (p1s_s2 - rho_zeta * p1s_s3);
                          dp_s3  += del_23 * inv_two_zeta * (p1s_s3 - rho_zeta * p1s_s4);
                        }

                        scalar_type mo_pre_term = norm * prefactor_dens * dens[mo_dens_ind];
                        mo_dens_ind++;

                        uint fit_dens_ind = 0;
                        //#pragma unroll 3
                        for (uint d3_l1 = 0; d3_l1 < 3; d3_l1++) {

                          scalar_type p1s_p1_1 = WmQ[d3_l1] * p1s_s2;
                          scalar_type p1s_p1_2 = WmQ[d3_l1] * p1s_s3;

                          scalar_type p2s_p1_1 = WmQ[d3_l1] * p2s_s2;
                          scalar_type p2s_p1_2 = WmQ[d3_l1] * p2s_s3;

                          scalar_type sp_p1_1  = WmQ[d3_l1] * sp_s2;
                          scalar_type sp_p1_2  = WmQ[d3_l1] * sp_s3;

                          scalar_type ds_p1_1  = WmQ[d3_l1] * ds_s2;
                          scalar_type ds_p1_2  = WmQ[d3_l1] * ds_s3;

                          scalar_type p1p_p1_1 = WmQ[d3_l1] * p1p_s2;
                          scalar_type p1p_p1_2 = WmQ[d3_l1] * p1p_s3;

                          scalar_type p2p_p1_1 = WmQ[d3_l1] * p2p_s2;
                          scalar_type p2p_p1_2 = WmQ[d3_l1] * p2p_s3;

                          scalar_type dp_p1_0  = WmQ[d3_l1] * dp_s1;
                          scalar_type dp_p1_1  = WmQ[d3_l1] * dp_s2;
                          scalar_type dp_p1_2  = WmQ[d3_l1] * dp_s3;
                          {
                            bool del_14 = d1_l1 == d3_l1;
                            p1s_p1_1 += del_14 * inv_two_zeta_eta * F_mT[2];
                            p1s_p1_2 += del_14 * inv_two_zeta_eta * F_mT[3];
                            ds_p1_1  += del_14 * inv_two_zeta_eta * p2s_s2;
                            ds_p1_2  += del_14 * inv_two_zeta_eta * p2s_s3;
                            p1p_p1_1 += del_14 * inv_two_zeta_eta * sp_s2;
                            p1p_p1_2 += del_14 * inv_two_zeta_eta * sp_s3;
                            dp_p1_0  += del_14 * inv_two_zeta_eta * p2p_s1;
                            dp_p1_1  += del_14 * inv_two_zeta_eta * p2p_s2;
                            dp_p1_2  += del_14 * inv_two_zeta_eta * p2p_s3;
                          }
                          {
                            bool del_24 = d1_l2 == d3_l1;
                            p2s_p1_1 += del_24 * inv_two_zeta_eta * F_mT[2];
                            p2s_p1_2 += del_24 * inv_two_zeta_eta * F_mT[3];
                            ds_p1_1  += del_24 * inv_two_zeta_eta * p1s_s2;
                            ds_p1_2  += del_24 * inv_two_zeta_eta * p1s_s3;
                            p2p_p1_1 += del_24 * inv_two_zeta_eta * sp_s2;
                            p2p_p1_2 += del_24 * inv_two_zeta_eta * sp_s3;
                            dp_p1_0  += del_24 * inv_two_zeta_eta * p1p_s1;
                            dp_p1_1  += del_24 * inv_two_zeta_eta * p1p_s2;
                            dp_p1_2  += del_24 * inv_two_zeta_eta * p1p_s3;
                          }
                          {
                            bool del_34 = p_l == d3_l1;
                            sp_p1_1  += del_34 * inv_two_zeta_eta * F_mT[2];
                            sp_p1_2  += del_34 * inv_two_zeta_eta * F_mT[3];
                            p1p_p1_1 += del_34 * inv_two_zeta_eta * p1s_s2;
                            p1p_p1_2 += del_34 * inv_two_zeta_eta * p1s_s3;
                            p2p_p1_1 += del_34 * inv_two_zeta_eta * p2s_s2;
                            p2p_p1_2 += del_34 * inv_two_zeta_eta * p2s_s3;
                            dp_p1_0  += del_34 * inv_two_zeta_eta * ds_s1;
                            dp_p1_1  += del_34 * inv_two_zeta_eta * ds_s2;
                            dp_p1_2  += del_34 * inv_two_zeta_eta * ds_s3;
                          }

                          for (uint d3_l2 = 0; d3_l2 <= d3_l1; d3_l2++) {

                            scalar_type ds_d0   = WmQ[d3_l2] * ds_p1_1;
                            scalar_type ds_d1   = WmQ[d3_l2] * ds_p1_2;

                            scalar_type dp_p2_0 = WmQ[d3_l2] * dp_s1;
                            scalar_type dp_p2_1 = WmQ[d3_l2] * dp_s2;

                            scalar_type p2p_d0  = WmQ[d3_l2] * p2p_p1_1;
                            scalar_type p2p_d1  = WmQ[d3_l2] * p2p_p1_2;

                            scalar_type p1p_d0  = WmQ[d3_l2] * p1p_p1_1;
                            scalar_type p1p_d1  = WmQ[d3_l2] * p1p_p1_2;

                            scalar_type dp_d0   = WmQ[d3_l2] * dp_p1_1;
                            scalar_type dp_d1   = WmQ[d3_l2] * dp_p1_2;

                            {
                              bool del_15 = d1_l1 == d3_l2;
                              ds_d0   += del_15 * inv_two_zeta_eta * p2s_p1_1;
                              ds_d1   += del_15 * inv_two_zeta_eta * p2s_p1_2;
                              dp_p2_0 += del_15 * inv_two_zeta_eta * p2p_s1;
                              dp_p2_1 += del_15 * inv_two_zeta_eta * p2p_s2;
                              p1p_d0  += del_15 * inv_two_zeta_eta * sp_p1_1;
                              p1p_d1  += del_15 * inv_two_zeta_eta * sp_p1_2;
                              dp_d0   += del_15 * inv_two_zeta_eta * p2p_p1_1;
                              dp_d1   += del_15 * inv_two_zeta_eta * p2p_p1_2;
                            }
                            {
                              bool del_25 = d1_l2 == d3_l2;
                              ds_d0   += del_25 * inv_two_zeta_eta * p1s_p1_1;
                              ds_d1   += del_25 * inv_two_zeta_eta * p1s_p1_2;
                              dp_p2_0 += del_25 * inv_two_zeta_eta * p1p_s1;
                              dp_p2_1 += del_25 * inv_two_zeta_eta * p1p_s2;
                              p2p_d0  += del_25 * inv_two_zeta_eta * sp_p1_1;
                              p2p_d1  += del_25 * inv_two_zeta_eta * sp_p1_2;
                              dp_d0   += del_25 * inv_two_zeta_eta * p1p_p1_1;
                              dp_d1   += del_25 * inv_two_zeta_eta * p1p_p1_2;
                            }
                            {
                              bool del_35 = p_l == d3_l2;
                              dp_p2_0 += del_35 * inv_two_zeta_eta * ds_s1;
                              dp_p2_1 += del_35 * inv_two_zeta_eta * ds_s2;
                              p2p_d0  += del_35 * inv_two_zeta_eta * p2s_p1_1;
                              p2p_d1  += del_35 * inv_two_zeta_eta * p2s_p1_2;
                              p1p_d0  += del_35 * inv_two_zeta_eta * p1s_p1_1;
                              p1p_d1  += del_35 * inv_two_zeta_eta * p1s_p1_2;
                              dp_d0   += del_35 * inv_two_zeta_eta * ds_p1_1;
                              dp_d1   += del_35 * inv_two_zeta_eta * ds_p1_2;
                            }
                            scalar_type pre_term;
                            {
                              bool del_45 = d3_l1 == d3_l2;
                              ds_d0  += del_45 * inv_two_ak * (ds_s0 - rho_ak * ds_s1);
                              ds_d1  += del_45 * inv_two_ak * (ds_s1 - rho_ak * ds_s2);
                              p2p_d0 += del_45 * inv_two_ak * (p2p_s0 - rho_ak * p2p_s1);
                              p2p_d1 += del_45 * inv_two_ak * (p2p_s1 - rho_ak * p2p_s2);
                              p1p_d0 += del_45 * inv_two_ak * (p1p_s0 - rho_ak * p1p_s1);
                              p1p_d1 += del_45 * inv_two_ak * (p1p_s1 - rho_ak * p1p_s2);
                              dp_d0  += del_45 * inv_two_ak * (dp_s0 - rho_ak * dp_s1);
                              dp_d1  += del_45 * inv_two_ak * (dp_s1 - rho_ak * dp_s2);

                              pre_term = del_45 * gpu_normalization_factor + !del_45 * 1.0f;
                            }
                            pre_term *= mo_pre_term * fit_dens_sh[j+fit_dens_ind];
                            fit_dens_ind++;

                            #pragma unroll 3
                            for (uint grad_l = 0; grad_l < 3; grad_l++) {
                              bool del_1g = d1_l1 == grad_l, del_2g = d1_l2 == grad_l, del_3g = p_l == grad_l, del_4g = d3_l1 == grad_l, del_5g = d3_l2 == grad_l;

                              A_force_term  = WmP[grad_l] * dp_d1;
                              A_force_term += del_1g * inv_two_zeta * (p2p_d0 - rho_zeta * p2p_d1);
                              A_force_term += del_2g * inv_two_zeta * (p1p_d0 - rho_zeta * p1p_d1);
                              A_force_term += del_3g * inv_two_zeta * (ds_d0 - rho_zeta * ds_d1);
                              A_force_term += del_4g * inv_two_zeta_eta * dp_p2_1;
                              A_force_term += del_5g * inv_two_zeta_eta * dp_p1_1;
                              B_force_term  = PmB[grad_l] * dp_d0 + A_force_term;
                              A_force_term  = PmA[grad_l] * dp_d0 + A_force_term;

                              C_force_term  = WmQ[grad_l] * dp_d1;
                              C_force_term += del_1g * inv_two_zeta_eta * p2p_d1;
                              C_force_term += del_2g * inv_two_zeta_eta * p1p_d1;
                              C_force_term += del_3g * inv_two_zeta_eta * ds_d1;
                              C_force_term += del_4g * inv_two_ak * (dp_p2_0 - rho_ak * dp_p2_1);
                              C_force_term += del_5g * inv_two_ak * (dp_p1_0 - rho_ak * dp_p1_1);
  
                              A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_1g * p2p_d0 - del_2g * p1p_d0);
                              B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - del_3g * ds_d0);
                              C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_4g * dp_p2_0 - del_5g * dp_p1_0);
                            }
                          }
                        }
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo;
                C_force[1][tid] *= valid_thread * prefactor_mo;
                C_force[2][tid] *= valid_thread * prefactor_mo;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (DP-D)----------------------------------------------
