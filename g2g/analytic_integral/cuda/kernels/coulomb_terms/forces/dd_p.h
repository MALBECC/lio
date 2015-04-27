//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DD-P)-------------------------------------------
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

                      scalar_type norm1;
                      {
                        bool del_12 = d1_l1 == d1_l2;
                        ds_s0 += del_12 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        ds_s1 += del_12 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        ds_s2 += del_12 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                        ds_s3 += del_12 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                        ds_s4 += del_12 * inv_two_zeta * (F_mT[4] - rho_zeta * F_mT[5]);
                        norm1 = del_12 * G2G::gpu_normalization_factor + !del_12 * 1.0f;
                      }

                      //#pragma unroll 3
                      for (uint d2_l1 = 0; d2_l1 < 3; d2_l1++) {

                        scalar_type sp1_s1  = PmB[d2_l1] * F_mT[1] + WmP[d2_l1] * F_mT[2];
                        scalar_type sp1_s2  = PmB[d2_l1] * F_mT[2] + WmP[d2_l1] * F_mT[3];
                        scalar_type sp1_s3  = PmB[d2_l1] * F_mT[3] + WmP[d2_l1] * F_mT[4];

                        scalar_type p1p1_s0 = PmB[d2_l1] * p1s_s0 + WmP[d2_l1] * p1s_s1;
                        scalar_type p1p1_s1 = PmB[d2_l1] * p1s_s1 + WmP[d2_l1] * p1s_s2;
                        scalar_type p1p1_s2 = PmB[d2_l1] * p1s_s2 + WmP[d2_l1] * p1s_s3;
                        scalar_type p1p1_s3 = PmB[d2_l1] * p1s_s3 + WmP[d2_l1] * p1s_s4;

                        scalar_type p2p1_s0 = PmB[d2_l1] * p2s_s0 + WmP[d2_l1] * p2s_s1;
                        scalar_type p2p1_s1 = PmB[d2_l1] * p2s_s1 + WmP[d2_l1] * p2s_s2;
                        scalar_type p2p1_s2 = PmB[d2_l1] * p2s_s2 + WmP[d2_l1] * p2s_s3;
                        scalar_type p2p1_s3 = PmB[d2_l1] * p2s_s3 + WmP[d2_l1] * p2s_s4;

                        scalar_type dp1_s0  = PmB[d2_l1] * ds_s0 + WmP[d2_l1] * ds_s1;
                        scalar_type dp1_s1  = PmB[d2_l1] * ds_s1 + WmP[d2_l1] * ds_s2;
                        scalar_type dp1_s2  = PmB[d2_l1] * ds_s2 + WmP[d2_l1] * ds_s3;
                        scalar_type dp1_s3  = PmB[d2_l1] * ds_s3 + WmP[d2_l1] * ds_s4;

                        {
                          bool del_13 = d1_l1 == d2_l1;
                          p1p1_s0 += del_13 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p1p1_s1 += del_13 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p1p1_s2 += del_13 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          p1p1_s3 += del_13 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                          dp1_s0  += del_13 * inv_two_zeta * (p2s_s0 - rho_zeta * p2s_s1);
                          dp1_s1  += del_13 * inv_two_zeta * (p2s_s1 - rho_zeta * p2s_s2);
                          dp1_s2  += del_13 * inv_two_zeta * (p2s_s2 - rho_zeta * p2s_s3);
                          dp1_s3  += del_13 * inv_two_zeta * (p2s_s3 - rho_zeta * p2s_s4);
                        }
                        {
                          bool del_23 = d1_l2 == d2_l1;
                          p2p1_s0 += del_23 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p2p1_s1 += del_23 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p2p1_s2 += del_23 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          p2p1_s3 += del_23 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                          dp1_s0  += del_23 * inv_two_zeta * (p1s_s0 - rho_zeta * p1s_s1);
                          dp1_s1  += del_23 * inv_two_zeta * (p1s_s1 - rho_zeta * p1s_s2);
                          dp1_s2  += del_23 * inv_two_zeta * (p1s_s2 - rho_zeta * p1s_s3);
                          dp1_s3  += del_23 * inv_two_zeta * (p1s_s3 - rho_zeta * p1s_s4);
                        }

                        for (uint d2_l2 = 0; d2_l2 <= d2_l1; d2_l2++) {

                          scalar_type p1p2_s1 = PmB[d2_l2] * p1s_s1 + WmP[d2_l2] * p1s_s2;
                          scalar_type p1p2_s2 = PmB[d2_l2] * p1s_s2 + WmP[d2_l2] * p1s_s3;

                          scalar_type p2p2_s1 = PmB[d2_l2] * p2s_s1 + WmP[d2_l2] * p2s_s2;
                          scalar_type p2p2_s2 = PmB[d2_l2] * p2s_s2 + WmP[d2_l2] * p2s_s3;

                          scalar_type sd_s1   = PmB[d2_l2] * sp1_s1 + WmP[d2_l2] * sp1_s2;
                          scalar_type sd_s2   = PmB[d2_l2] * sp1_s2 + WmP[d2_l2] * sp1_s3;

                          scalar_type dp2_s1  = PmB[d2_l2] * ds_s1 + WmP[d2_l2] * ds_s2;
                          scalar_type dp2_s2  = PmB[d2_l2] * ds_s2 + WmP[d2_l2] * ds_s3;

                          scalar_type p1d_s1  = PmB[d2_l2] * p1p1_s1 + WmP[d2_l2] * p1p1_s2;
                          scalar_type p1d_s2  = PmB[d2_l2] * p1p1_s2 + WmP[d2_l2] * p1p1_s3;

                          scalar_type p2d_s1  = PmB[d2_l2] * p2p1_s1 + WmP[d2_l2] * p2p1_s2;
                          scalar_type p2d_s2  = PmB[d2_l2] * p2p1_s2 + WmP[d2_l2] * p2p1_s3;

                          scalar_type dd_s0   = PmB[d2_l2] * dp1_s0 + WmP[d2_l2] * dp1_s1;
                          scalar_type dd_s1   = PmB[d2_l2] * dp1_s1 + WmP[d2_l2] * dp1_s2;
                          scalar_type dd_s2   = PmB[d2_l2] * dp1_s2 + WmP[d2_l2] * dp1_s3;

                          {
                            bool del_14 = d1_l1 == d2_l2;
                            p1p2_s1 += del_14 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                            p1p2_s2 += del_14 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                            dp2_s1  += del_14 * inv_two_zeta * (p2s_s1 - rho_zeta * p2s_s2);
                            dp2_s2  += del_14 * inv_two_zeta * (p2s_s2 - rho_zeta * p2s_s3);
                            p1d_s1  += del_14 * inv_two_zeta * (sp1_s1 - rho_zeta * sp1_s2);
                            p1d_s2  += del_14 * inv_two_zeta * (sp1_s2 - rho_zeta * sp1_s3);
                            dd_s0   += del_14 * inv_two_zeta * (p2p1_s0 - rho_zeta * p2p1_s1);
                            dd_s1   += del_14 * inv_two_zeta * (p2p1_s1 - rho_zeta * p2p1_s2);
                            dd_s2   += del_14 * inv_two_zeta * (p2p1_s2 - rho_zeta * p2p1_s3);
                          }
                          {
                            bool del_24 = d1_l2 == d2_l2;
                            p2p2_s1 += del_24 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                            p2p2_s2 += del_24 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                            dp2_s1  += del_24 * inv_two_zeta * (p1s_s1 - rho_zeta * p1s_s2);
                            dp2_s2  += del_24 * inv_two_zeta * (p1s_s2 - rho_zeta * p1s_s3);
                            p2d_s1  += del_24 * inv_two_zeta * (sp1_s1 - rho_zeta * sp1_s2);
                            p2d_s2  += del_24 * inv_two_zeta * (sp1_s2 - rho_zeta * sp1_s3);
                            dd_s0   += del_24 * inv_two_zeta * (p1p1_s0 - rho_zeta * p1p1_s1);
                            dd_s1   += del_24 * inv_two_zeta * (p1p1_s1 - rho_zeta * p1p1_s2);
                            dd_s2   += del_24 * inv_two_zeta * (p1p1_s2 - rho_zeta * p1p1_s3);
                          }
                          scalar_type norm2;
                          {
                            bool del_34 = d2_l1 == d2_l2;
                            sd_s1  += del_34 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                            sd_s2  += del_34 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                            p1d_s1 += del_34 * inv_two_zeta * (p1s_s1 - rho_zeta * p1s_s2);
                            p1d_s2 += del_34 * inv_two_zeta * (p1s_s2 - rho_zeta * p1s_s3);
                            p2d_s1 += del_34 * inv_two_zeta * (p2s_s1 - rho_zeta * p2s_s2);
                            p2d_s2 += del_34 * inv_two_zeta * (p2s_s2 - rho_zeta * p2s_s3);
                            dd_s0  += del_34 * inv_two_zeta * (ds_s0 - rho_zeta * ds_s1);
                            dd_s1  += del_34 * inv_two_zeta * (ds_s1 - rho_zeta * ds_s2);
                            dd_s2  += del_34 * inv_two_zeta * (ds_s2 - rho_zeta * ds_s3);

                            norm2 = del_34 * G2G::gpu_normalization_factor + !del_34 * 1.0f;
                          }

                          scalar_type mo_pre_term = norm1 * norm2 * prefactor_dens;
                          {
                            bool skip = same_func && (d2_l1 > d1_l1 || (d2_l1 == d1_l1 && d2_l2 > d1_l2));
                            mo_pre_term *= !skip * dens[mo_dens_ind];
                            mo_dens_ind += !skip;
                          }

                          #pragma unroll 3
                          for (uint p_l = 0; p_l < 3; p_l++) {

                            scalar_type dp1_p0 = WmQ[p_l] * dp1_s1;
                            scalar_type dp1_p1 = WmQ[p_l] * dp1_s2;

                            scalar_type dp2_p0 = WmQ[p_l] * dp2_s1;
                            scalar_type dp2_p1 = WmQ[p_l] * dp2_s2;

                            scalar_type p1d_p0 = WmQ[p_l] * p1d_s1;
                            scalar_type p1d_p1 = WmQ[p_l] * p1d_s2;

                            scalar_type p2d_p0 = WmQ[p_l] * p2d_s1;
                            scalar_type p2d_p1 = WmQ[p_l] * p2d_s2;

                            scalar_type dd_p0  = WmQ[p_l] * dd_s1;
                            scalar_type dd_p1  = WmQ[p_l] * dd_s2;

                            {
                              bool del_15 = d1_l1 == p_l;
                              dp1_p0 += del_15 * inv_two_zeta_eta * p2p1_s1;
                              dp1_p1 += del_15 * inv_two_zeta_eta * p2p1_s2;
                              dp2_p0 += del_15 * inv_two_zeta_eta * p2p2_s1;
                              dp2_p1 += del_15 * inv_two_zeta_eta * p2p2_s2;
                              p1d_p0 += del_15 * inv_two_zeta_eta * sd_s1;
                              p1d_p1 += del_15 * inv_two_zeta_eta * sd_s2;
                              dd_p0  += del_15 * inv_two_zeta_eta * p2d_s1;
                              dd_p1  += del_15 * inv_two_zeta_eta * p2d_s2;
                            }
                            {
                              bool del_25 = d1_l2 == p_l;
                              dp1_p0 += del_25 * inv_two_zeta_eta * p1p1_s1;
                              dp1_p1 += del_25 * inv_two_zeta_eta * p1p1_s2;
                              dp2_p0 += del_25 * inv_two_zeta_eta * p1p2_s1;
                              dp2_p1 += del_25 * inv_two_zeta_eta * p1p2_s2;
                              p2d_p0 += del_25 * inv_two_zeta_eta * sd_s1;
                              p2d_p1 += del_25 * inv_two_zeta_eta * sd_s2;
                              dd_p0  += del_25 * inv_two_zeta_eta * p1d_s1;
                              dd_p1  += del_25 * inv_two_zeta_eta * p1d_s2;
                            }
                            {
                              bool del_35 = d2_l1 == p_l;
                              dp1_p0 += del_35 * inv_two_zeta_eta * ds_s1;
                              dp1_p1 += del_35 * inv_two_zeta_eta * ds_s2;
                              p1d_p0 += del_35 * inv_two_zeta_eta * p1p2_s1;
                              p1d_p1 += del_35 * inv_two_zeta_eta * p1p2_s2;
                              p2d_p0 += del_35 * inv_two_zeta_eta * p2p2_s1;
                              p2d_p1 += del_35 * inv_two_zeta_eta * p2p2_s2;
                              dd_p0  += del_35 * inv_two_zeta_eta * dp2_s1;
                              dd_p1  += del_35 * inv_two_zeta_eta * dp2_s2;
                            }
                            {
                              bool del_45 = d2_l2 == p_l;
                              dp2_p0 += del_45 * inv_two_zeta_eta * ds_s1;
                              dp2_p1 += del_45 * inv_two_zeta_eta * ds_s2;
                              p1d_p0 += del_45 * inv_two_zeta_eta * p1p1_s1;
                              p1d_p1 += del_45 * inv_two_zeta_eta * p1p1_s2;
                              p2d_p0 += del_45 * inv_two_zeta_eta * p2p1_s1;
                              p2d_p1 += del_45 * inv_two_zeta_eta * p2p1_s2;
                              dd_p0  += del_45 * inv_two_zeta_eta * dp1_s1;
                              dd_p1  += del_45 * inv_two_zeta_eta * dp1_s2;
                            }

                            scalar_type pre_term = mo_pre_term * fit_dens_sh[j+p_l];
 
                            #pragma unroll 3
                            for (uint grad_l = 0; grad_l < 3; grad_l++) {
                              bool del_1g = d1_l1 == grad_l, del_2g = d1_l2 == grad_l, del_3g = d2_l1 == grad_l, del_4g = d2_l2 == grad_l, del_5g = p_l == grad_l;

                              A_force_term  = WmP[grad_l] * dd_p1;
                              A_force_term += del_1g * inv_two_zeta * (p2d_p0 - rho_zeta * p2d_p1);
                              A_force_term += del_2g * inv_two_zeta * (p1d_p0 - rho_zeta * p1d_p1);
                              A_force_term += del_3g * inv_two_zeta * (dp2_p0 - rho_zeta * dp2_p1);
                              A_force_term += del_4g * inv_two_zeta * (dp1_p0 - rho_zeta * dp1_p1);
                              A_force_term += del_5g * inv_two_zeta_eta * dd_s1;
                              B_force_term  = PmB[grad_l] * dd_p0 + A_force_term;
                              A_force_term  = PmA[grad_l] * dd_p0 + A_force_term;

                              C_force_term  = WmQ[grad_l] * dd_p1;
                              C_force_term += del_1g * inv_two_zeta_eta * p2d_p1;
                              C_force_term += del_2g * inv_two_zeta_eta * p1d_p1;
                              C_force_term += del_3g * inv_two_zeta_eta * dp2_p1;
                              C_force_term += del_4g * inv_two_zeta_eta * dp1_p1;
                              C_force_term += del_5g * inv_two_ak * (dd_s0 - rho_ak * dd_s1);
  
                              A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_1g * p2d_p0 - del_2g * p1d_p0);
                              B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - del_3g * dp2_p0 - del_4g * dp1_p0);
                              C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_5g * dd_s0);
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
//------------------------------------------END TERM-TYPE DEPENDENT PART (DD-P)----------------------------------------------
