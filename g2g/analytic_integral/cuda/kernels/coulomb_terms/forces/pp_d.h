//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (PP-D)-------------------------------------------
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
                  for (uint p1_l = 0; p1_l < 3; p1_l++) {

                    scalar_type ps_s0 = PmA[p1_l] * F_mT[0] + WmP[p1_l] * F_mT[1];
                    scalar_type ps_s1 = PmA[p1_l] * F_mT[1] + WmP[p1_l] * F_mT[2];
                    scalar_type ps_s2 = PmA[p1_l] * F_mT[2] + WmP[p1_l] * F_mT[3];
                    scalar_type ps_s3 = PmA[p1_l] * F_mT[3] + WmP[p1_l] * F_mT[4];
                    scalar_type ps_s4 = PmA[p1_l] * F_mT[4] + WmP[p1_l] * F_mT[5];

                    #pragma unroll 3
                    for (uint p2_l = 0; p2_l < 3; p2_l++) {

                      scalar_type sp_s0 = PmB[p2_l] * F_mT[0] + WmP[p2_l] * F_mT[1];
                      scalar_type sp_s1 = PmB[p2_l] * F_mT[1] + WmP[p2_l] * F_mT[2];
                      scalar_type sp_s2 = PmB[p2_l] * F_mT[2] + WmP[p2_l] * F_mT[3];
                      scalar_type sp_s3 = PmB[p2_l] * F_mT[3] + WmP[p2_l] * F_mT[4];

                      scalar_type pp_s0 = PmB[p2_l] * ps_s0 + WmP[p2_l] * ps_s1;
                      scalar_type pp_s1 = PmB[p2_l] * ps_s1 + WmP[p2_l] * ps_s2;
                      scalar_type pp_s2 = PmB[p2_l] * ps_s2 + WmP[p2_l] * ps_s3;
                      scalar_type pp_s3 = PmB[p2_l] * ps_s3 + WmP[p2_l] * ps_s4;
                      {
                        bool del_12 = p1_l == p2_l;
                        pp_s0 += del_12 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        pp_s1 += del_12 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        pp_s2 += del_12 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                        pp_s3 += del_12 * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);
                      }
                      scalar_type mo_pre_term;
                      {
                        bool skip = same_func && (p2_l > p1_l);
                        mo_pre_term = !skip * prefactor_dens * dens[mo_dens_ind];
                        mo_dens_ind += !skip;
                      }

                      uint fit_dens_ind = 0;
                      //#pragma unroll 3
                      for (uint d_l1 = 0; d_l1 < 3; d_l1++) {

                        scalar_type ss_p1_1 = WmQ[d_l1] * F_mT[2];
                        scalar_type ss_p1_2 = WmQ[d_l1] * F_mT[3];

                        //scalar_type ps_p1_0 = WmQ[d_l1] * ps_s1;
                        scalar_type ps_p1_1 = WmQ[d_l1] * ps_s2;
                        scalar_type ps_p1_2 = WmQ[d_l1] * ps_s3;

                        //scalar_type sp_p1_0 = WmQ[d_l1] * sp_s1;
                        scalar_type sp_p1_1 = WmQ[d_l1] * sp_s2;
                        scalar_type sp_p1_2 = WmQ[d_l1] * sp_s3;

                        scalar_type pp_p1_0 = WmQ[d_l1] * pp_s1;
                        scalar_type pp_p1_1 = WmQ[d_l1] * pp_s2;
                        scalar_type pp_p1_2 = WmQ[d_l1] * pp_s3;
                        {
                          bool del_13 = p1_l == d_l1;
                          //ps_p1_0 += del_13 * inv_two_zeta_eta * F_mT[1];
                          ps_p1_1 += del_13 * inv_two_zeta_eta * F_mT[2];
                          ps_p1_2 += del_13 * inv_two_zeta_eta * F_mT[3];
                          pp_p1_0 += del_13 * inv_two_zeta_eta * sp_s1;
                          pp_p1_1 += del_13 * inv_two_zeta_eta * sp_s2;
                          pp_p1_2 += del_13 * inv_two_zeta_eta * sp_s3;
                        }
                        {
                          bool del_23 = p2_l == d_l1;
                          //sp_p1_0 += del_23 * inv_two_zeta_eta * F_mT[1];
                          sp_p1_1 += del_23 * inv_two_zeta_eta * F_mT[2];
                          sp_p1_2 += del_23 * inv_two_zeta_eta * F_mT[3];
                          pp_p1_0 += del_23 * inv_two_zeta_eta * ps_s1;
                          pp_p1_1 += del_23 * inv_two_zeta_eta * ps_s2;
                          pp_p1_2 += del_23 * inv_two_zeta_eta * ps_s3;
                        }

                        for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {

                          scalar_type pp_p2_0 = WmQ[d_l2] * pp_s1;
                          scalar_type pp_p2_1 = WmQ[d_l2] * pp_s2;

                          scalar_type sp_d0   = WmQ[d_l2] * sp_p1_1;
                          scalar_type sp_d1   = WmQ[d_l2] * sp_p1_2;

                          scalar_type ps_d0   = WmQ[d_l2] * ps_p1_1;
                          scalar_type ps_d1   = WmQ[d_l2] * ps_p1_2;

                          scalar_type pp_d0   = WmQ[d_l2] * pp_p1_1;
                          scalar_type pp_d1   = WmQ[d_l2] * pp_p1_2;

                          {
                            bool del_14 = p1_l == d_l2;
                            pp_p2_0 += del_14 * inv_two_zeta_eta * sp_s1;
                            pp_p2_1 += del_14 * inv_two_zeta_eta * sp_s2;
                            ps_d0   += del_14 * inv_two_zeta_eta * ss_p1_1;
                            ps_d1   += del_14 * inv_two_zeta_eta * ss_p1_2;
                            pp_d0   += del_14 * inv_two_zeta_eta * sp_p1_1;
                            pp_d1   += del_14 * inv_two_zeta_eta * sp_p1_2;
                          }
                          {
                            bool del_24 = p2_l == d_l2;
                            pp_p2_0 += del_24 * inv_two_zeta_eta * ps_s1;
                            pp_p2_1 += del_24 * inv_two_zeta_eta * ps_s2;
                            sp_d0   += del_24 * inv_two_zeta_eta * ss_p1_1;
                            sp_d1   += del_24 * inv_two_zeta_eta * ss_p1_2;
                            pp_d0   += del_24 * inv_two_zeta_eta * ps_p1_1;
                            pp_d1   += del_24 * inv_two_zeta_eta * ps_p1_2;
                          }
                          scalar_type norm;
                          {
                            bool del_34 = d_l1 == d_l2;
                            sp_d0   += del_34 * inv_two_ak * (sp_s0 - rho_ak * sp_s1);
                            sp_d1   += del_34 * inv_two_ak * (sp_s1 - rho_ak * sp_s2);
                            ps_d0   += del_34 * inv_two_ak * (ps_s0 - rho_ak * ps_s1);
                            ps_d1   += del_34 * inv_two_ak * (ps_s1 - rho_ak * ps_s2);
                            pp_d0   += del_34 * inv_two_ak * (pp_s0 - rho_ak * pp_s1);
                            pp_d1   += del_34 * inv_two_ak * (pp_s1 - rho_ak * pp_s2);

                            norm = del_34 * G2G::gpu_normalization_factor + !del_34 * 1.0f;
                          }

                          scalar_type pre_term = norm * mo_pre_term * fit_dens_sh[j+fit_dens_ind];
                          fit_dens_ind++;

                          #pragma unroll 3
                          for (uint grad_l = 0; grad_l < 3; grad_l++) {
                            bool del_1g = p1_l == grad_l, del_2g = p2_l == grad_l, del_3g = d_l1 == grad_l, del_4g = d_l2 == grad_l;

                            A_force_term  = WmP[grad_l] * pp_d1;
                            A_force_term += del_1g * inv_two_zeta * (sp_d0 - rho_zeta * sp_d1);
                            A_force_term += del_2g * inv_two_zeta * (ps_d0 - rho_zeta * ps_d1);
                            A_force_term += del_3g * inv_two_zeta_eta * pp_p2_1;
                            A_force_term += del_4g * inv_two_zeta_eta * pp_p1_1;
                            B_force_term  = PmB[grad_l] * pp_d0 + A_force_term;
                            A_force_term  = PmA[grad_l] * pp_d0 + A_force_term;

                            C_force_term  = WmQ[grad_l] * pp_d1;
                            C_force_term += del_1g * inv_two_zeta_eta * sp_d1;
                            C_force_term += del_2g * inv_two_zeta_eta * ps_d1;
                            C_force_term += del_3g * inv_two_ak * (pp_p2_0 - rho_ak * pp_p2_1);
                            C_force_term += del_4g * inv_two_ak * (pp_p1_0 - rho_ak * pp_p1_1);

                            A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_1g * sp_d0);
                            B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - del_2g * ps_d0);
                            C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_3g * pp_p2_0 - del_4g * pp_p1_0);
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
//------------------------------------------END TERM-TYPE DEPENDENT PART (PP-D)----------------------------------------------
