//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (PS-P)-------------------------------------------
              {
                scalar_type F_mT[4];
                {
                  scalar_type PmQ[3];
                  PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
                  PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
                  PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
                  scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
                  lio_gamma<scalar_type,3>(F_mT,T);
                }
                {
                  scalar_type A_force_term, B_force_term, C_force_term;
                  #pragma unroll 3
                  for (uint p1_l = 0; p1_l < 3; p1_l++) {

                    scalar_type ps_s0 = PmA[p1_l] * F_mT[0] + WmP[p1_l] * F_mT[1];
                    scalar_type ps_s1 = PmA[p1_l] * F_mT[1] + WmP[p1_l] * F_mT[2];
                    scalar_type ps_s2 = PmA[p1_l] * F_mT[2] + WmP[p1_l] * F_mT[3];

                    #pragma unroll 3
                    for (uint p3_l = 0; p3_l < 3; p3_l++) {

                      scalar_type ss_p0 = WmQ[p3_l] * F_mT[1];
                      scalar_type ss_p1 = WmQ[p3_l] * F_mT[2];

                      scalar_type ps_p0 = WmQ[p3_l] * ps_s1;
                      scalar_type ps_p1 = WmQ[p3_l] * ps_s2;
                      {
                        bool del_p1p3 = p1_l == p3_l;
                        ps_p0 += del_p1p3 * inv_two_zeta_eta * F_mT[1];
                        ps_p1 += del_p1p3 * inv_two_zeta_eta * F_mT[2];
                      }

                      scalar_type pre_term = prefactor_dens * dens[p1_l] * fit_dens_sh[j+p3_l];

                      #pragma unroll 3
                      for (uint grad_l = 0; grad_l < 3; grad_l++) {
                        bool del_p1g = p1_l == grad_l, del_p3g = p3_l == grad_l;

                        A_force_term  = WmP[grad_l] * ps_p1;
                        A_force_term += del_p1g * inv_two_zeta * (ss_p0 - rho_zeta * ss_p1);
                        A_force_term += del_p3g * inv_two_zeta_eta * ps_s1;
                        B_force_term  = PmB[grad_l] * ps_p0 + A_force_term;
                        A_force_term  = PmA[grad_l] * ps_p0 + A_force_term;

                        C_force_term  = WmQ[grad_l] * ps_p1;
                        C_force_term += del_p1g * inv_two_zeta_eta * ss_p1;
                        C_force_term += del_p3g * inv_two_ak * (ps_s0 - rho_ak * ps_s1);

                        A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_p1g * ss_p0);
                        B_force[grad_l]      += pre_term * 2.0f * aj * B_force_term;
                        C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_p3g * ps_s0);
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo;
                C_force[1][tid] *= valid_thread * prefactor_mo;
                C_force[2][tid] *= valid_thread * prefactor_mo;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (PS-P)----------------------------------------------
