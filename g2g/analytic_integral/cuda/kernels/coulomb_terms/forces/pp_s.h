//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (PP-S)-------------------------------------------
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
                  uint dens_ind = 0;
                  #pragma unroll 3
                  for (uint p1_l = 0; p1_l < 3; p1_l++) {

                    scalar_type p_s0 = PmA[p1_l] * F_mT[0] + WmP[p1_l] * F_mT[1];
		    scalar_type p_s1 = PmA[p1_l] * F_mT[1] + WmP[p1_l] * F_mT[2];
                    scalar_type p_s2 = PmA[p1_l] * F_mT[2] + WmP[p1_l] * F_mT[3];

                    #pragma unroll 3
                    for (uint p2_l = 0; p2_l < 3; p2_l++) {
                      bool del_p = p1_l == p2_l;

                      scalar_type s_p0 = PmB[p2_l] * F_mT[0] + WmP[p2_l] * F_mT[1];
                      scalar_type s_p1 = PmB[p2_l] * F_mT[1] + WmP[p2_l] * F_mT[2];

                      scalar_type p_p0 = PmB[p2_l] * p_s0 + WmP[p2_l] * p_s1;
                      p_p0            += del_p * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                      scalar_type p_p1 = PmB[p2_l] * p_s1 + WmP[p2_l] * p_s2;
                      p_p1            += del_p * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);

                      bool skip = same_func && (p2_l > p1_l);
                      scalar_type pre_term = !skip * prefactor_dens * dens[dens_ind];
                      dens_ind += !skip;

                      #pragma unroll 3
                      for (uint grad_l = 0; grad_l < 3; grad_l++) {
                        bool del_p1g = p1_l == grad_l, del_p2g = p2_l == grad_l;

                        C_force_term  = del_p1g * s_p1;
                        C_force_term += del_p2g * p_s1;
                        C_force_term *= inv_two_zeta_eta;
                        C_force_term += WmQ[grad_l] * p_p1;

                        A_force_term  = del_p1g * (s_p0 - rho_zeta * s_p1);
                        A_force_term += del_p2g * (p_s0 - rho_zeta * p_s1);
                        A_force_term *= inv_two_zeta;
                        A_force_term += WmP[grad_l] * p_p1;

                        B_force_term  = PmB[grad_l] * p_p0 + A_force_term;
                        A_force_term  = PmA[grad_l] * p_p0 + A_force_term;

                        A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - del_p1g * s_p0);
                        B_force[grad_l]      += pre_term * (2.0f * aj * B_force_term - del_p2g * p_s0);
                        C_force[grad_l][tid] += pre_term * C_force_term;
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[1][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[2][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (PP-S)----------------------------------------------
