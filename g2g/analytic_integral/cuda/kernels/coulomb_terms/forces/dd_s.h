//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DD-S)-------------------------------------------
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
                  uint dens_ind = 0;
                  //#pragma unroll 3
                  for (uint d1_l1 = 0; d1_l1 < 3; d1_l1++) {

                    scalar_type p1_s0 = PmA[d1_l1] * F_mT[0] + WmP[d1_l1] * F_mT[1];
		    scalar_type p1_s1 = PmA[d1_l1] * F_mT[1] + WmP[d1_l1] * F_mT[2];
                    scalar_type p1_s2 = PmA[d1_l1] * F_mT[2] + WmP[d1_l1] * F_mT[3];
                    scalar_type p1_s3 = PmA[d1_l1] * F_mT[3] + WmP[d1_l1] * F_mT[4];
                    scalar_type p1_s4 = PmA[d1_l1] * F_mT[4] + WmP[d1_l1] * F_mT[5];

                    for (uint d1_l2 = 0; d1_l2 <= d1_l1; d1_l2++) {

                      scalar_type p2_s0 = PmA[d1_l2] * F_mT[0] + WmP[d1_l2] * F_mT[1];
		      scalar_type p2_s1 = PmA[d1_l2] * F_mT[1] + WmP[d1_l2] * F_mT[2];
		      scalar_type p2_s2 = PmA[d1_l2] * F_mT[2] + WmP[d1_l2] * F_mT[3];
		      scalar_type p2_s3 = PmA[d1_l2] * F_mT[3] + WmP[d1_l2] * F_mT[4];

                      scalar_type d_s0 = PmA[d1_l2] * p1_s0 + WmP[d1_l2] * p1_s1;
                      scalar_type d_s1 = PmA[d1_l2] * p1_s1 + WmP[d1_l2] * p1_s2;
                      scalar_type d_s2 = PmA[d1_l2] * p1_s2 + WmP[d1_l2] * p1_s3;
                      scalar_type d_s3 = PmA[d1_l2] * p1_s3 + WmP[d1_l2] * p1_s4;

                      scalar_type norm1;
                      {
                        bool del_d = d1_l1 == d1_l2;
                        d_s0 += del_d * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                        d_s1 += del_d * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                        d_s2 += del_d * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                        d_s3 += del_d * inv_two_zeta * (F_mT[3] - rho_zeta * F_mT[4]);

                        norm1 = del_d * G2G::gpu_normalization_factor + !del_d * 1.0f;
                      }

                      //#pragma unroll 3
                      for (uint d2_l1 = 0; d2_l1 < 3; d2_l1++) {

                        scalar_type s_p1_0 = PmB[d2_l1] * F_mT[0] + WmP[d2_l1] * F_mT[1];
                        scalar_type s_p1_1 = PmB[d2_l1] * F_mT[1] + WmP[d2_l1] * F_mT[2];
                        scalar_type s_p1_2 = PmB[d2_l1] * F_mT[2] + WmP[d2_l1] * F_mT[3];

                        scalar_type d_p1_0 = 0.0f, d_p1_1 = 0.0f, d_p1_2 = 0.0f;

                        scalar_type p1_p1_0 = PmB[d2_l1] * p1_s0 + WmP[d2_l1] * p1_s1;
                        scalar_type p1_p1_1 = PmB[d2_l1] * p1_s1 + WmP[d2_l1] * p1_s2;
                        scalar_type p1_p1_2 = PmB[d2_l1] * p1_s2 + WmP[d2_l1] * p1_s3;
                        scalar_type p2_p1_0 = PmB[d2_l1] * p2_s0 + WmP[d2_l1] * p2_s1;
                        scalar_type p2_p1_1 = PmB[d2_l1] * p2_s1 + WmP[d2_l1] * p2_s2;
                        scalar_type p2_p1_2 = PmB[d2_l1] * p2_s2 + WmP[d2_l1] * p2_s3;

                        {
                          bool del_dp1 = d1_l1 == d2_l1;
                          p1_p1_0 += del_dp1 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p1_p1_1 += del_dp1 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p1_p1_2 += del_dp1 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          d_p1_0  += del_dp1 * (p2_s0 - rho_zeta * p2_s1);
                          d_p1_1  += del_dp1 * (p2_s1 - rho_zeta * p2_s2);
                          d_p1_2  += del_dp1 * (p2_s2 - rho_zeta * p2_s3);
                        }
                        {
                          bool del_dp2 = d1_l2 == d2_l1;
                          p2_p1_0 += del_dp2 * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                          p2_p1_1 += del_dp2 * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
                          p2_p1_2 += del_dp2 * inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
                          d_p1_0  += del_dp2 * (p1_s0 - rho_zeta * p1_s1);
                          d_p1_1  += del_dp2 * (p1_s1 - rho_zeta * p1_s2);
                          d_p1_2  += del_dp2 * (p1_s2 - rho_zeta * p1_s3);
                        }

                        d_p1_0 *= inv_two_zeta;
                        d_p1_1 *= inv_two_zeta;
                        d_p1_2 *= inv_two_zeta;
                        d_p1_0 += PmB[d2_l1] * d_s0 + WmP[d2_l1] * d_s1;
                        d_p1_1 += PmB[d2_l1] * d_s1 + WmP[d2_l1] * d_s2;
                        d_p1_2 += PmB[d2_l1] * d_s2 + WmP[d2_l1] * d_s3;

                        for (uint d2_l2 = 0; d2_l2 <= d2_l1; d2_l2++) {

                          scalar_type d_p2_0 = 0.0f, d_p2_1 = 0.0f;
                          scalar_type p1_d_0 = 0.0f, p1_d_1 = 0.0f;
                          scalar_type p2_d_0 = 0.0f, p2_d_1 = 0.0f;
                          scalar_type d_d0 = 0.0f, d_d1 = 0.0f;
                          scalar_type norm2;
                          {
                            bool del_dp1 = d1_l1 == d2_l2;
                            d_p2_0  += del_dp1 * (p2_s0 - rho_zeta * p2_s1);
                            d_p2_1  += del_dp1 * (p2_s1 - rho_zeta * p2_s2);
                            p1_d_0  += del_dp1 * (s_p1_0 - rho_zeta * s_p1_1);
                            p1_d_1  += del_dp1 * (s_p1_1 - rho_zeta * s_p1_2);
                            d_d0    += del_dp1 * (p2_p1_0 - rho_zeta * p2_p1_1);
                            d_d1    += del_dp1 * (p2_p1_1 - rho_zeta * p2_p1_2);
                          }
                          {
                            bool del_dp2 = d1_l2 == d2_l2;
                            d_p2_0  += del_dp2 * (p1_s0 - rho_zeta * p1_s1);
                            d_p2_1  += del_dp2 * (p1_s1 - rho_zeta * p1_s2);
                            p2_d_0  += del_dp2 * (s_p1_0 - rho_zeta * s_p1_1);
                            p2_d_1  += del_dp2 * (s_p1_1 - rho_zeta * s_p1_2);
                            d_d0    += del_dp2 * (p1_p1_0 - rho_zeta * p1_p1_1);
                            d_d1    += del_dp2 * (p1_p1_1 - rho_zeta * p1_p1_2);
                          }
                          {
                            bool del_d = d2_l1 == d2_l2;
                            p1_d_0  += del_d * (p1_s0 - rho_zeta * p1_s1);
                            p1_d_1  += del_d * (p1_s1 - rho_zeta * p1_s2);
                            p2_d_0  += del_d * (p2_s0 - rho_zeta * p2_s1);
                            p2_d_1  += del_d * (p2_s1 - rho_zeta * p2_s2);
                            d_d0    += del_d * (d_s0 - rho_zeta * d_s1);
                            d_d1    += del_d * (d_s1 - rho_zeta * d_s2);
                            norm2    = del_d * G2G::gpu_normalization_factor + !del_d * 1.0f;
                          }
                          d_p2_0 *= inv_two_zeta;
                          d_p2_1 *= inv_two_zeta;
                          p1_d_0 *= inv_two_zeta;
                          p1_d_1 *= inv_two_zeta;
                          p2_d_0 *= inv_two_zeta;
                          p2_d_1 *= inv_two_zeta;
                          d_d0   *= inv_two_zeta;
                          d_d1   *= inv_two_zeta;
                          d_p2_0 += PmB[d2_l2] * d_s0 + WmP[d2_l2] * d_s1;
                          d_p2_1 += PmB[d2_l2] * d_s1 + WmP[d2_l2] * d_s2;
                          p1_d_0 += PmB[d2_l2] * p1_p1_0 + WmP[d2_l2] * p1_p1_1;
                          p1_d_1 += PmB[d2_l2] * p1_p1_1 + WmP[d2_l2] * p1_p1_2;
                          p2_d_0 += PmB[d2_l2] * p2_p1_0 + WmP[d2_l2] * p2_p1_1;
                          p2_d_1 += PmB[d2_l2] * p2_p1_1 + WmP[d2_l2] * p2_p1_2;
                          d_d0   += PmB[d2_l2] * d_p1_0 + WmP[d2_l2] * d_p1_1;
                          d_d1   += PmB[d2_l2] * d_p1_1 + WmP[d2_l2] * d_p1_2;

                          bool skip = same_func && (d2_l1 > d1_l1 || ((d2_l1 == d1_l1) && d2_l2 > d1_l2));
                          scalar_type pre_term = !skip * norm1 * norm2 * prefactor_dens * dens[dens_ind];
                          dens_ind += !skip;

                          #pragma unroll 3
                          for (uint grad_l = 0; grad_l < 3; grad_l++) {
                            bool del_d11g = d1_l1 == grad_l, del_d12g = d1_l2 == grad_l, del_d21g = d2_l1 == grad_l, del_d22g = d2_l2 == grad_l;

                            C_force_term  = del_d11g * p2_d_1;
                            C_force_term += del_d12g * p1_d_1;
                            C_force_term += del_d21g * d_p2_1;
                            C_force_term += del_d22g * d_p1_1;
                            C_force_term *= inv_two_zeta_eta;
                            C_force_term += WmQ[grad_l] * d_d1;

                            A_force_term  = del_d11g * (p2_d_0 - rho_zeta * p2_d_1);
                            A_force_term += del_d12g * (p1_d_0 - rho_zeta * p1_d_1);
                            A_force_term += del_d21g * (d_p2_0 - rho_zeta * d_p2_1);
                            A_force_term += del_d22g * (d_p1_0 - rho_zeta * d_p1_1);
                            A_force_term *= inv_two_zeta;
                            A_force_term += WmP[grad_l] * d_d1;

                            B_force_term  = PmB[grad_l] * d_d0 + A_force_term;
                            A_force_term  = PmA[grad_l] * d_d0 + A_force_term;

                            A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - (del_d11g * p2_d_0 + del_d12g * p1_d_0));
                            B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - (del_d21g * d_p2_0 + del_d22g * d_p1_0));
                            C_force[grad_l][tid] += pre_term * C_force_term;
                          }
                        }
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[1][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[2][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (DD-S)----------------------------------------------
