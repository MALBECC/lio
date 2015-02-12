//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DP-S)-------------------------------------------
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
                  uint dens_ind = 0;
                  //#pragma unroll 3
                  for (uint d_l1 = 0; d_l1 < 3; d_l1++) {

                    scalar_type p1_s0 = PmA[d_l1] * F_mT[0] + WmP[d_l1] * F_mT[1];
		    scalar_type p1_s1 = PmA[d_l1] * F_mT[1] + WmP[d_l1] * F_mT[2];
                    scalar_type p1_s2 = PmA[d_l1] * F_mT[2] + WmP[d_l1] * F_mT[3];
                    scalar_type p1_s3 = PmA[d_l1] * F_mT[3] + WmP[d_l1] * F_mT[4];

                    for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {

                      scalar_type p2_s0 = PmA[d_l2] * F_mT[0] + WmP[d_l2] * F_mT[1];
		      scalar_type p2_s1 = PmA[d_l2] * F_mT[1] + WmP[d_l2] * F_mT[2];
		      scalar_type p2_s2 = PmA[d_l2] * F_mT[2] + WmP[d_l2] * F_mT[3];

                      scalar_type d_s0 = PmA[d_l2] * p1_s0 + WmP[d_l2] * p1_s1;
                      scalar_type d_s1 = PmA[d_l2] * p1_s1 + WmP[d_l2] * p1_s2;
                      scalar_type d_s2 = PmA[d_l2] * p1_s2 + WmP[d_l2] * p1_s3;

                      scalar_type norm;
                      {
                        bool del_d = d_l1 == d_l2;
                        d_s0            += del_d * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        d_s1            += del_d * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        d_s2            += del_d * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);

                        norm = del_d * gpu_normalization_factor + !del_d * 1.0f;
                      }

                      #pragma unroll 3
                      for (uint p_l = 0; p_l < 3; p_l++) {

                        scalar_type d_p0 = 0.0f, d_p1 = 0.0f;

                        scalar_type p1_p0 = PmB[p_l] * p1_s0 + WmP[p_l] * p1_s1;
                        scalar_type p1_p1 = PmB[p_l] * p1_s1 + WmP[p_l] * p1_s2;
                        scalar_type p2_p0 = PmB[p_l] * p2_s0 + WmP[p_l] * p2_s1;
                        scalar_type p2_p1 = PmB[p_l] * p2_s1 + WmP[p_l] * p2_s2;

                        {
                          bool del_dp1 = d_l1 == p_l;
                          p1_p0            += del_dp1 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p1_p1            += del_dp1 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          d_p0             += del_dp1 * (p2_s0 - rho_zeta * p2_s1);
                          d_p1             += del_dp1 * (p2_s1 - rho_zeta * p2_s2);
                        }
                        {
                          bool del_dp2 = d_l2 == p_l;
                          p2_p0            += del_dp2 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p2_p1            += del_dp2 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          d_p0             += del_dp2 * (p1_s0 - rho_zeta * p1_s1);
                          d_p1             += del_dp2 * (p1_s1 - rho_zeta * p1_s2);
                        }

                        d_p0             *= inv_two_zeta;
                        d_p0             += PmB[p_l] * d_s0 + WmP[p_l] * d_s1;
                        d_p1             *= inv_two_zeta;
                        d_p1             += PmB[p_l] * d_s1 + WmP[p_l] * d_s2;

                        scalar_type pre_term = norm * prefactor_dens * dens[dens_ind];
                        dens_ind++;

                        #pragma unroll 3
                        for (uint grad_l = 0; grad_l < 3; grad_l++) {
                          bool del_d1g = d_l1 == grad_l, del_d2g = d_l2 == grad_l, del_pg = p_l == grad_l;

                          C_force_term  = del_d1g * p2_p1;
                          C_force_term += del_d2g * p1_p1;
                          C_force_term += del_pg * d_s1;
                          C_force_term *= inv_two_zeta_eta;
                          C_force_term += WmQ[grad_l] * d_p1;

                          A_force_term  = del_d1g * (p2_p0 - rho_zeta * p2_p1);
                          A_force_term += del_d2g * (p1_p0 - rho_zeta * p1_p1);
                          A_force_term += del_pg * (d_s0 - rho_zeta * d_s1);
                          A_force_term *= inv_two_zeta;
                          A_force_term += WmP[grad_l] * d_p1;

                          B_force_term  = PmB[grad_l] * d_p0 + A_force_term;
                          A_force_term  = PmA[grad_l] * d_p0 + A_force_term;

                          A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - (del_d1g * p2_p0 + del_d2g * p1_p0));
                          B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - del_pg * d_s0);
                          C_force[grad_l][tid] += pre_term * C_force_term;
                        }
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[1][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[2][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (DP-S)----------------------------------------------
