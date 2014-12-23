//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DS-D)-------------------------------------------
              {
                scalar_type F_mT[6];
                {
                  scalar_type PmQ[3];
                  PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
                  PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
                  PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
                  scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
                  lio_gamma<scalar_type,5>(F_mT,T);
                }
                {
                  scalar_type A_force_term, B_force_term, C_force_term;
                  uint mo_dens_ind = 0;
                  #pragma unroll 3
                  for (uint d1_l1 = 0; d1_l1 < 3; d1_l1++) {

                    scalar_type p1s_s0 = PmA[d1_l1] * F_mT[0] + WmP[d1_l1] * F_mT[1];
                    scalar_type p1s_s1 = PmA[d1_l1] * F_mT[1] + WmP[d1_l1] * F_mT[2];
                    scalar_type p1s_s2 = PmA[d1_l1] * F_mT[2] + WmP[d1_l1] * F_mT[3];
                    scalar_type p1s_s3 = PmA[d1_l1] * F_mT[3] + WmP[d1_l1] * F_mT[4];
                    scalar_type p1s_s4 = PmA[d1_l1] * F_mT[4] + WmP[d1_l1] * F_mT[5];

                    for (uint d1_l2 = 0; d1_l2 <= d1_l1; d1_l2++) {

                      scalar_type p2s_s0 = PmA[d1_l2] * F_mT[0] + WmP[d1_l2] * F_mT[1];
                      scalar_type p2s_s1 = PmA[d1_l2] * F_mT[1] + WmP[d1_l2] * F_mT[2];
                      scalar_type p2s_s2 = PmA[d1_l2] * F_mT[2] + WmP[d1_l2] * F_mT[3];
                      scalar_type p2s_s3 = PmA[d1_l2] * F_mT[3] + WmP[d1_l2] * F_mT[4];

                      scalar_type ds_s0  = PmA[d1_l2] * p1s_s0 + WmP[d1_l2] * p1s_s1;
                      scalar_type ds_s1  = PmA[d1_l2] * p1s_s1 + WmP[d1_l2] * p1s_s2;
                      scalar_type ds_s2  = PmA[d1_l2] * p1s_s2 + WmP[d1_l2] * p1s_s3;
                      scalar_type ds_s3  = PmA[d1_l2] * p1s_s3 + WmP[d1_l2] * p1s_s4;

                      scalar_type mo_pre_term;
                      {
                        bool del_12 = d1_l1 == d1_l2;
                        ds_s0 += del_12 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        ds_s1 += del_12 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        ds_s2 += del_12 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                        ds_s3 += del_12 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);

                        mo_pre_term = del_12 * gpu_normalization_factor + !del_12 * 1.0f;
                      }
                      mo_pre_term *= prefactor_dens * dens[mo_dens_ind];
                      mo_dens_ind++;

                      uint fit_dens_ind = 0;
                      #pragma unroll 3
                      for (uint d3_l1 = 0; d3_l1 < 3; d3_l1++) {

                        scalar_type ss_p1_1  = WmQ[d3_l1] * F_mT[2];
                        scalar_type ss_p1_2  = WmQ[d3_l1] * F_mT[3];

                        scalar_type p1s_p1_1 = WmQ[d3_l1] * p1s_s2;
                        scalar_type p1s_p1_2 = WmQ[d3_l1] * p1s_s3;

                        scalar_type p2s_p1_1 = WmQ[d3_l1] * p2s_s2;
                        scalar_type p2s_p1_2 = WmQ[d3_l1] * p2s_s3;

                        scalar_type ds_p1_0  = WmQ[d3_l1] * ds_s1;
                        scalar_type ds_p1_1  = WmQ[d3_l1] * ds_s2;
                        scalar_type ds_p1_2  = WmQ[d3_l1] * ds_s3;
                        {
                          bool del_13 = d1_l1 == d3_l1;
                          p1s_p1_1 += del_13 * inv_two_zeta_eta * F_mT[2];
                          p1s_p1_2 += del_13 * inv_two_zeta_eta * F_mT[3];
                          ds_p1_0  += del_13 * inv_two_zeta_eta * p2s_s1;
                          ds_p1_1  += del_13 * inv_two_zeta_eta * p2s_s2;
                          ds_p1_2  += del_13 * inv_two_zeta_eta * p2s_s3;
                        }
                        {
                          bool del_23 = d1_l2 == d3_l1;
                          p2s_p1_1 += del_23 * inv_two_zeta_eta * F_mT[2];
                          p2s_p1_2 += del_23 * inv_two_zeta_eta * F_mT[3];
                          ds_p1_0 += del_23 * inv_two_zeta_eta * p1s_s1;
                          ds_p1_1 += del_23 * inv_two_zeta_eta * p1s_s2;
                          ds_p1_2 += del_23 * inv_two_zeta_eta * p1s_s3;
                        }

                        for (uint d3_l2 = 0; d3_l2 <= d3_l1; d3_l2++) {

                          scalar_type ds_p2_0 = WmQ[d3_l2] * ds_s1;
                          scalar_type ds_p2_1 = WmQ[d3_l2] * ds_s2;

                          scalar_type p2s_d0  = WmQ[d3_l2] * p2s_p1_1;
                          scalar_type p2s_d1  = WmQ[d3_l2] * p2s_p1_2;

                          scalar_type p1s_d0  = WmQ[d3_l2] * p1s_p1_1;
                          scalar_type p1s_d1  = WmQ[d3_l2] * p1s_p1_2;

                          scalar_type ds_d0   = WmQ[d3_l2] * ds_p1_1;
                          scalar_type ds_d1   = WmQ[d3_l2] * ds_p1_2;

                          {
                            bool del_14 = d1_l1 == d3_l2;
                            ds_p2_0 += del_14 * inv_two_zeta_eta * p2s_s1;
                            ds_p2_1 += del_14 * inv_two_zeta_eta * p2s_s2;
                            p1s_d0  += del_14 * inv_two_zeta_eta * ss_p1_1;
                            p1s_d1  += del_14 * inv_two_zeta_eta * ss_p1_2;
                            ds_d0   += del_14 * inv_two_zeta_eta * p2s_p1_1;
                            ds_d1   += del_14 * inv_two_zeta_eta * p2s_p1_2;
                          }
                          {
                            bool del_24 = d1_l2 == d3_l2;
                            ds_p2_0 += del_24 * inv_two_zeta_eta * p1s_s1;
                            ds_p2_1 += del_24 * inv_two_zeta_eta * p1s_s2;
                            p2s_d0  += del_24 * inv_two_zeta_eta * ss_p1_1;
                            p2s_d1  += del_24 * inv_two_zeta_eta * ss_p1_2;
                            ds_d0   += del_24 * inv_two_zeta_eta * p1s_p1_1;
                            ds_d1   += del_24 * inv_two_zeta_eta * p1s_p1_2;
                          }

                          scalar_type pre_term;
                          {
                            bool del_34 = d3_l1 == d3_l2;
                            p2s_d0  += del_34 * inv_two_ak * (p2s_s0 - rho_ak * p2s_s1);
                            p2s_d1  += del_34 * inv_two_ak * (p2s_s1 - rho_ak * p2s_s2);
                            p1s_d0  += del_34 * inv_two_ak * (p1s_s0 - rho_ak * p1s_s1);
                            p1s_d1  += del_34 * inv_two_ak * (p1s_s1 - rho_ak * p1s_s2);
                            ds_d0   += del_34 * inv_two_ak * (ds_s0 - rho_ak * ds_s1);
                            ds_d1   += del_34 * inv_two_ak * (ds_s1 - rho_ak * ds_s2);

                            pre_term = del_34 * gpu_normalization_factor + !del_34 * 1.0f;
                          }

                          pre_term *= mo_pre_term * fit_dens_sh[j+fit_dens_ind];
                          fit_dens_ind++;

                          #pragma unroll 3
                          for (uint grad_l = 0; grad_l < 3; grad_l++) {
                            bool del_1g = d1_l1 == grad_l, del_2g = d1_l2 == grad_l, del_3g = d3_l1 == grad_l, del_4g = d3_l2 == grad_l;

                            A_force_term  = WmP[grad_l] * ds_d1;
                            A_force_term += del_1g * inv_two_zeta * (p2s_d0 - rho_zeta * p2s_d1);
                            A_force_term += del_2g * inv_two_zeta * (p1s_d0 - rho_zeta * p1s_d1);
                            A_force_term += del_3g * inv_two_zeta_eta * ds_p2_1;
                            A_force_term += del_4g * inv_two_zeta_eta * ds_p1_1;
                            B_force_term  = PmB[grad_l] * ds_d0 + A_force_term;
                            A_force_term  = PmA[grad_l] * ds_d0 + A_force_term;

                            C_force_term  = WmQ[grad_l] * ds_d1;
                            C_force_term += del_1g * inv_two_zeta_eta * p2s_d1;
                            C_force_term += del_2g * inv_two_zeta_eta * p1s_d1;
                            C_force_term += del_3g * inv_two_ak * (ds_p2_0 - rho_ak * ds_p2_1);
                            C_force_term += del_4g * inv_two_ak * (ds_p1_0 - rho_ak * ds_p1_1);

                            A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_1g * p2s_d0 - del_2g * p1s_d0);
                            B_force[grad_l]      += pre_term * 2.0f * aj * B_force_term;
                            C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_3g * ds_p2_0 - del_4g * ds_p1_0);
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
//------------------------------------------END TERM-TYPE DEPENDENT PART (DS-D)----------------------------------------------
