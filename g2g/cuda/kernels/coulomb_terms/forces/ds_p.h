//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DS-P)-------------------------------------------
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
                  uint mo_dens_ind = 0;
                  //#pragma unroll 3
                  for (uint d_l1 = 0; d_l1 < 3; d_l1++) {

                    scalar_type p1s_s0 = PmA[d_l1] * F_mT[0] + WmP[d_l1] * F_mT[1];
                    scalar_type p1s_s1 = PmA[d_l1] * F_mT[1] + WmP[d_l1] * F_mT[2];
                    scalar_type p1s_s2 = PmA[d_l1] * F_mT[2] + WmP[d_l1] * F_mT[3];
                    scalar_type p1s_s3 = PmA[d_l1] * F_mT[3] + WmP[d_l1] * F_mT[4];

                    //#pragma unroll 3
                    for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {

                      scalar_type p2s_s1 = PmA[d_l2] * F_mT[1] + WmP[d_l2] * F_mT[2];
                      scalar_type p2s_s2 = PmA[d_l2] * F_mT[2] + WmP[d_l2] * F_mT[3];

                      scalar_type ds_s0 = PmA[d_l2] * p1s_s0 + WmP[d_l2] * p1s_s1;
                      scalar_type ds_s1 = PmA[d_l2] * p1s_s1 + WmP[d_l2] * p1s_s2;
                      scalar_type ds_s2 = PmA[d_l2] * p1s_s2 + WmP[d_l2] * p1s_s3;

                      scalar_type mo_pre_term;
                      {
                        bool del_12 = d_l1 == d_l2;
                        ds_s0 += del_12 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        ds_s1 += del_12 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        ds_s2 += del_12 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                        mo_pre_term = del_12 * gpu_normalization_factor + !del_12 * 1.0f;
                      }

                      mo_pre_term *= prefactor_dens * dens[mo_dens_ind];
                      mo_dens_ind++;

                      #pragma unroll 3
                      for (uint p_l = 0; p_l < 3; p_l++) {

                        scalar_type p1s_p0 = WmQ[p_l] * p1s_s1;
                        scalar_type p1s_p1 = WmQ[p_l] * p1s_s2;

                        scalar_type p2s_p0 = WmQ[p_l] * p2s_s1;
                        scalar_type p2s_p1 = WmQ[p_l] * p2s_s2;

                        scalar_type ds_p0 = WmQ[p_l] * ds_s1;
                        scalar_type ds_p1 = WmQ[p_l] * ds_s2;
                        {
                          bool del_13 = d_l1 == p_l;
                          p1s_p0 += del_13 * inv_two_zeta_eta * F_mT[1];
                          p1s_p1 += del_13 * inv_two_zeta_eta * F_mT[2];
                          ds_p0  += del_13 * inv_two_zeta_eta * p2s_s1;
                          ds_p1  += del_13 * inv_two_zeta_eta * p2s_s2;
                        }
                        {
                          bool del_23 = d_l2 == p_l;
                          p2s_p0 += del_23 * inv_two_zeta_eta * F_mT[1];
                          p2s_p1 += del_23 * inv_two_zeta_eta * F_mT[2];
                          ds_p0  += del_23 * inv_two_zeta_eta * p1s_s1;
                          ds_p1  += del_23 * inv_two_zeta_eta * p1s_s2;
                        }

                        scalar_type pre_term = mo_pre_term * fit_dens_sh[j+p_l];

                        #pragma unroll 3
                        for (uint grad_l = 0; grad_l < 3; grad_l++) {
                          bool del_d1g = d_l1 == grad_l, del_d2g = d_l2 == grad_l, del_pg = p_l == grad_l;

                          A_force_term  = WmP[grad_l] * ds_p1;
                          A_force_term += del_d1g * inv_two_zeta * (p2s_p0 - rho_zeta * p2s_p1);
                          A_force_term += del_d2g * inv_two_zeta * (p1s_p0 - rho_zeta * p1s_p1);
                          A_force_term += del_pg * inv_two_zeta_eta * ds_s1;
                          B_force_term  = PmB[grad_l] * ds_p0 + A_force_term;
                          A_force_term  = PmA[grad_l] * ds_p0 + A_force_term;

                          C_force_term  = WmQ[grad_l] * ds_p1;
                          C_force_term += del_d1g * inv_two_zeta_eta * p2s_p1;
                          C_force_term += del_d2g * inv_two_zeta_eta * p1s_p1;
                          C_force_term += del_pg * inv_two_ak * (ds_s0 - rho_ak * ds_s1);

                          A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_d1g * p2s_p0 - del_d2g * p1s_p0);
                          B_force[grad_l]      += pre_term * 2.0f * aj * B_force_term;
                          C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_pg * ds_s0);
                        }
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo;
                C_force[1][tid] *= valid_thread * prefactor_mo;
                C_force[2][tid] *= valid_thread * prefactor_mo;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (DS-P)----------------------------------------------
