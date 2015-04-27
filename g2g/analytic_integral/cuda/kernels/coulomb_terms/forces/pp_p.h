//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (PP-P)-------------------------------------------
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
                  #pragma unroll 3
                  for (uint p1_l = 0; p1_l < 3; p1_l++) {

                    scalar_type ps_s0 = PmA[p1_l] * F_mT[0] + WmP[p1_l] * F_mT[1];
                    scalar_type ps_s1 = PmA[p1_l] * F_mT[1] + WmP[p1_l] * F_mT[2];
                    scalar_type ps_s2 = PmA[p1_l] * F_mT[2] + WmP[p1_l] * F_mT[3];
                    scalar_type ps_s3 = PmA[p1_l] * F_mT[3] + WmP[p1_l] * F_mT[4];

                    #pragma unroll 3
                    for (uint p2_l = 0; p2_l < 3; p2_l++) {

                      scalar_type sp_s1 = PmB[p2_l] * F_mT[1] + WmP[p2_l] * F_mT[2];
                      scalar_type sp_s2 = PmB[p2_l] * F_mT[2] + WmP[p2_l] * F_mT[3];

                      scalar_type pp_s0 = PmB[p2_l] * ps_s0 + WmP[p2_l] * ps_s1;
                      scalar_type pp_s1 = PmB[p2_l] * ps_s1 + WmP[p2_l] * ps_s2;
                      scalar_type pp_s2 = PmB[p2_l] * ps_s2 + WmP[p2_l] * ps_s3;
                      {
                        bool del_12 = p1_l == p2_l;
                        pp_s0 += del_12 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        pp_s1 += del_12 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        pp_s2 += del_12 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                      }
                      scalar_type mo_pre_term;
                      {
                        bool skip = same_func && (p2_l > p1_l);
                        mo_pre_term = !skip * prefactor_dens * dens[mo_dens_ind];
                        mo_dens_ind += !skip;
                      }

                      #pragma unroll 3
                      for (uint p3_l = 0; p3_l < 3; p3_l++) {

                        scalar_type ps_p0 = WmQ[p3_l] * ps_s1;
                        scalar_type ps_p1 = WmQ[p3_l] * ps_s2;

                        scalar_type sp_p0 = WmQ[p3_l] * sp_s1;
                        scalar_type sp_p1 = WmQ[p3_l] * sp_s2;

                        scalar_type pp_p0 = WmQ[p3_l] * pp_s1;
                        scalar_type pp_p1 = WmQ[p3_l] * pp_s2;
                        {
                          bool del_13 = p1_l == p3_l;
                          ps_p0 += del_13 * inv_two_zeta_eta * F_mT[1];
                          ps_p1 += del_13 * inv_two_zeta_eta * F_mT[2];
                          pp_p0 += del_13 * inv_two_zeta_eta * sp_s1;
                          pp_p1 += del_13 * inv_two_zeta_eta * sp_s2;
                        }
                        {
                          bool del_23 = p2_l == p3_l;
                          sp_p0 += del_23 * inv_two_zeta_eta * F_mT[1];
                          sp_p1 += del_23 * inv_two_zeta_eta * F_mT[2];
                          pp_p0 += del_23 * inv_two_zeta_eta * ps_s1;
                          pp_p1 += del_23 * inv_two_zeta_eta * ps_s2;
                        }

                        scalar_type pre_term = mo_pre_term * fit_dens_sh[j+p3_l];

                        #pragma unroll 3
                        for (uint grad_l = 0; grad_l < 3; grad_l++) {
                          bool del_p1g = p1_l == grad_l, del_p2g = p2_l == grad_l, del_p3g = p3_l == grad_l;

                          A_force_term  = WmP[grad_l] * pp_p1;
                          A_force_term += del_p1g * inv_two_zeta * (sp_p0 - rho_zeta * sp_p1);
                          A_force_term += del_p2g * inv_two_zeta * (ps_p0 - rho_zeta * ps_p1);
                          A_force_term += del_p3g * inv_two_zeta_eta * pp_s1;
                          B_force_term  = PmB[grad_l] * pp_p0 + A_force_term;
                          A_force_term  = PmA[grad_l] * pp_p0 + A_force_term;

                          C_force_term  = WmQ[grad_l] * pp_p1;
                          C_force_term += del_p1g * inv_two_zeta_eta * sp_p1;
                          C_force_term += del_p2g * inv_two_zeta_eta * ps_p1;
                          C_force_term += del_p3g * inv_two_ak * (pp_s0 - rho_ak * pp_s1);

                          A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_p1g * sp_p0);
                          B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - del_p2g * ps_p0);
                          C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_p3g * pp_s0);
                        }
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo;
                C_force[1][tid] *= valid_thread * prefactor_mo;
                C_force[2][tid] *= valid_thread * prefactor_mo;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (PP-P)----------------------------------------------
