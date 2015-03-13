//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DP-P)-------------------------------------------
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
                  //#pragma unroll 3
                  for (uint d_l1 = 0; d_l1 < 3; d_l1++) {

                    scalar_type p1s_s0 = PmA[d_l1] * F_mT[0] + WmP[d_l1] * F_mT[1];
                    scalar_type p1s_s1 = PmA[d_l1] * F_mT[1] + WmP[d_l1] * F_mT[2];
                    scalar_type p1s_s2 = PmA[d_l1] * F_mT[2] + WmP[d_l1] * F_mT[3];
                    scalar_type p1s_s3 = PmA[d_l1] * F_mT[3] + WmP[d_l1] * F_mT[4];
                    scalar_type p1s_s4 = PmA[d_l1] * F_mT[4] + WmP[d_l1] * F_mT[5];

                    //#pragma unroll 3
                    for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {

                      scalar_type p2s_s0 = PmA[d_l2] * F_mT[0] + WmP[d_l2] * F_mT[1];
                      scalar_type p2s_s1 = PmA[d_l2] * F_mT[1] + WmP[d_l2] * F_mT[2];
                      scalar_type p2s_s2 = PmA[d_l2] * F_mT[2] + WmP[d_l2] * F_mT[3];
                      scalar_type p2s_s3 = PmA[d_l2] * F_mT[3] + WmP[d_l2] * F_mT[4];

                      scalar_type ds_s0 = PmA[d_l2] * p1s_s0 + WmP[d_l2] * p1s_s1;
                      scalar_type ds_s1 = PmA[d_l2] * p1s_s1 + WmP[d_l2] * p1s_s2;
                      scalar_type ds_s2 = PmA[d_l2] * p1s_s2 + WmP[d_l2] * p1s_s3;
                      scalar_type ds_s3 = PmA[d_l2] * p1s_s3 + WmP[d_l2] * p1s_s4;

                      scalar_type norm;
                      {
                        bool del_12 = d_l1 == d_l2;
                        ds_s0 += del_12 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        ds_s1 += del_12 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        ds_s2 += del_12 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                        ds_s3 += del_12 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                        norm = del_12 * gpu_normalization_factor + !del_12 * 1.0f;
                      }

                      #pragma unroll 3
                      for (uint p2_l = 0; p2_l < 3; p2_l++) {

                        scalar_type sp_s1  = PmB[p2_l] * F_mT[1] + WmP[p2_l] * F_mT[2];
                        scalar_type sp_s2  = PmB[p2_l] * F_mT[2] + WmP[p2_l] * F_mT[3];

                        scalar_type p1p_s1 = PmB[p2_l] * p1s_s1 + WmP[p2_l] * p1s_s2;
                        scalar_type p1p_s2 = PmB[p2_l] * p1s_s2 + WmP[p2_l] * p1s_s3;

                        scalar_type p2p_s1 = PmB[p2_l] * p2s_s1 + WmP[p2_l] * p2s_s2;
                        scalar_type p2p_s2 = PmB[p2_l] * p2s_s2 + WmP[p2_l] * p2s_s3;

                        scalar_type dp_s0  = PmB[p2_l] * ds_s0 + WmP[p2_l] * ds_s1;
                        scalar_type dp_s1  = PmB[p2_l] * ds_s1 + WmP[p2_l] * ds_s2;
                        scalar_type dp_s2  = PmB[p2_l] * ds_s2 + WmP[p2_l] * ds_s3;

                        {
                          bool del_13 = d_l1 == p2_l;
                          p1p_s1 += del_13 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p1p_s2 += del_13 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          dp_s0  += del_13 * inv_two_zeta * (p2s_s0 - rho_zeta * p2s_s1);
                          dp_s1  += del_13 * inv_two_zeta * (p2s_s1 - rho_zeta * p2s_s2);
                          dp_s2  += del_13 * inv_two_zeta * (p2s_s2 - rho_zeta * p2s_s3);
                        }
                        {
                          bool del_23 = d_l2 == p2_l;
                          p2p_s1 += del_23 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p2p_s2 += del_23 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          dp_s0  += del_23 * inv_two_zeta * (p1s_s0 - rho_zeta * p1s_s1);
                          dp_s1  += del_23 * inv_two_zeta * (p1s_s1 - rho_zeta * p1s_s2);
                          dp_s2  += del_23 * inv_two_zeta * (p1s_s2 - rho_zeta * p1s_s3);
                        }

                        scalar_type mo_pre_term = norm * prefactor_dens * dens[mo_dens_ind];
                        mo_dens_ind++;

                        #pragma unroll 3
                        for (uint p3_l = 0; p3_l < 3; p3_l++) {

                          scalar_type ds_p0  = WmQ[p3_l] * ds_s1;
                          scalar_type ds_p1  = WmQ[p3_l] * ds_s2;

                          scalar_type p1p_p0 = WmQ[p3_l] * p1p_s1;
                          scalar_type p1p_p1 = WmQ[p3_l] * p1p_s2;

                          scalar_type p2p_p0 = WmQ[p3_l] * p2p_s1;
                          scalar_type p2p_p1 = WmQ[p3_l] * p2p_s2;

                          scalar_type dp_p0  = WmQ[p3_l] * dp_s1;
                          scalar_type dp_p1  = WmQ[p3_l] * dp_s2;
                          {
                            bool del_14 = d_l1 == p3_l;
                            ds_p0  += del_14 * inv_two_zeta_eta * p2s_s1;
                            ds_p1  += del_14 * inv_two_zeta_eta * p2s_s2;
                            p1p_p0 += del_14 * inv_two_zeta_eta * sp_s1;
                            p1p_p1 += del_14 * inv_two_zeta_eta * sp_s2;
                            dp_p0  += del_14 * inv_two_zeta_eta * p2p_s1;
                            dp_p1  += del_14 * inv_two_zeta_eta * p2p_s2;
                          }
                          {
                            bool del_24 = d_l2 == p3_l;
                            ds_p0  += del_24 * inv_two_zeta_eta * p1s_s1;
                            ds_p1  += del_24 * inv_two_zeta_eta * p1s_s2;
                            p2p_p0 += del_24 * inv_two_zeta_eta * sp_s1;
                            p2p_p1 += del_24 * inv_two_zeta_eta * sp_s2;
                            dp_p0  += del_24 * inv_two_zeta_eta * p1p_s1;
                            dp_p1  += del_24 * inv_two_zeta_eta * p1p_s2;
                          }
                          {
                            bool del_34 = p2_l == p3_l;
                            p1p_p0 += del_34 * inv_two_zeta_eta * p1s_s1;
                            p1p_p1 += del_34 * inv_two_zeta_eta * p1s_s2;
                            p2p_p0 += del_34 * inv_two_zeta_eta * p2s_s1;
                            p2p_p1 += del_34 * inv_two_zeta_eta * p2s_s2;
                            dp_p0  += del_34 * inv_two_zeta_eta * ds_s1;
                            dp_p1  += del_34 * inv_two_zeta_eta * ds_s2;
                          }

                          scalar_type pre_term = mo_pre_term * fit_dens_sh[j+p3_l];

                          #pragma unroll 3
                          for (uint grad_l = 0; grad_l < 3; grad_l++) {
                            bool del_d1g = d_l1 == grad_l, del_d2g = d_l2 == grad_l, del_p2g = p2_l == grad_l, del_p3g = p3_l == grad_l;

                            A_force_term  = WmP[grad_l] * dp_p1;
                            A_force_term += del_d1g * inv_two_zeta * (p2p_p0 - rho_zeta * p2p_p1);
                            A_force_term += del_d2g * inv_two_zeta * (p1p_p0 - rho_zeta * p1p_p1);
                            A_force_term += del_p2g * inv_two_zeta * (ds_p0 - rho_zeta * ds_p1);
                            A_force_term += del_p3g * inv_two_zeta_eta * dp_s1;
                            B_force_term  = PmB[grad_l] * dp_p0 + A_force_term;
                            A_force_term  = PmA[grad_l] * dp_p0 + A_force_term;

                            C_force_term  = WmQ[grad_l] * dp_p1;
                            C_force_term += del_d1g * inv_two_zeta_eta * p2p_p1;
                            C_force_term += del_d2g * inv_two_zeta_eta * p1p_p1;
                            C_force_term += del_p2g * inv_two_zeta_eta * ds_p1;
                            C_force_term += del_p3g * inv_two_ak * (dp_s0 - rho_ak * dp_s1);
  
                            A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_d1g * p2p_p0 - del_d2g * p1p_p0);
                            B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - del_p2g * ds_p0);
                            C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_p3g * dp_s0);
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
//------------------------------------------END TERM-TYPE DEPENDENT PART (DP-P)----------------------------------------------
