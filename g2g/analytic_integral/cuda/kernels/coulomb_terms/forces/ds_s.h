//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (DS-S)-------------------------------------------
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
                  //#pragma unroll 3
                  for (uint d_l1 = 0; d_l1 < 3; d_l1++) {

                    scalar_type p1_s0 = PmA[d_l1] * F_mT[0] + WmP[d_l1] * F_mT[1];
		    scalar_type p1_s1 = PmA[d_l1] * F_mT[1] + WmP[d_l1] * F_mT[2];
                    scalar_type p1_s2 = PmA[d_l1] * F_mT[2] + WmP[d_l1] * F_mT[3];

                    for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {
                      bool del_d = d_l1 == d_l2;

                      scalar_type p2_s0 = PmA[d_l2] * F_mT[0] + WmP[d_l2] * F_mT[1];
		      scalar_type p2_s1 = PmA[d_l2] * F_mT[1] + WmP[d_l2] * F_mT[2];

                      scalar_type d_s0 = PmA[d_l2] * p1_s0 + WmP[d_l2] * p1_s1;
                      d_s0            += del_d * inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
                      scalar_type d_s1 = PmA[d_l2] * p1_s1 + WmP[d_l2] * p1_s2;
                      d_s1            += del_d * inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);

                      scalar_type pre_term = prefactor_dens * dens[dens_ind];
                      pre_term *= del_d * G2G::gpu_normalization_factor + !del_d * 1.0f;
                      dens_ind++;

                      #pragma unroll 3
                      for (uint grad_l = 0; grad_l < 3; grad_l++) {
                        bool del_d1g = d_l1 == grad_l, del_d2g = d_l2 == grad_l;

                        C_force_term  = del_d1g * p2_s1;
                        C_force_term += del_d2g * p1_s1;
                        C_force_term *= inv_two_zeta_eta;
                        C_force_term += WmQ[grad_l] * d_s1;

                        A_force_term  = del_d1g * (p2_s0 - rho_zeta * p2_s1);
                        A_force_term += del_d2g * (p1_s0 - rho_zeta * p1_s1);
                        A_force_term *= inv_two_zeta;
                        A_force_term += WmP[grad_l] * d_s1;

                        B_force_term  = PmB[grad_l] * d_s0 + A_force_term;
                        A_force_term  = PmA[grad_l] * d_s0 + A_force_term;

                        A_force[grad_l]      += pre_term * (2.0f * ai * A_force_term - (del_d1g * p2_s0 + del_d2g * p1_s0));
                        B_force[grad_l]      += pre_term * 2.0f * aj * B_force_term;
                        C_force[grad_l][tid] += pre_term * C_force_term;
                      }
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[1][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
                C_force[2][tid] *= valid_thread * prefactor_mo * 2.0f * ac_val_dens_sh[j].x;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (DS-S)----------------------------------------------
