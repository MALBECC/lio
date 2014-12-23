//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (SS-P)-------------------------------------------
              {
                scalar_type F_mT[3];
                {
                  scalar_type PmQ[3];
                  PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
                  PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
                  PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
                  scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
                  lio_gamma<scalar_type,2>(F_mT,T);
                }
                {
                  scalar_type A_force_term, B_force_term, C_force_term;
                  #pragma unroll 3
                  for (uint p_l = 0; p_l < 3; p_l++) {
                    scalar_type ss_p0 = WmQ[p_l] * F_mT[1];
                    scalar_type ss_p1 = WmQ[p_l] * F_mT[2];
                    scalar_type pre_term = prefactor_dens * fit_dens_sh[j+p_l];
                    //if (valid_thread) printf("%.4e\n",fit_dens_sh[j+p_l]);
                    #pragma unroll 3
                    for (uint grad_l = 0; grad_l < 3; grad_l++) {
                      bool del_pg = p_l == grad_l;

                      A_force_term  = WmP[grad_l] * ss_p1;
                      A_force_term += del_pg * inv_two_zeta_eta * F_mT[1];
                      B_force_term  = PmB[grad_l] * ss_p0 + A_force_term;
                      A_force_term  = PmA[grad_l] * ss_p0 + A_force_term;

                      C_force_term  = WmQ[grad_l] * ss_p1;
                      C_force_term += del_pg * inv_two_ak * (F_mT[0] - rho_ak * F_mT[1]);

                      A_force[grad_l]      += pre_term * 2.0f * ai * A_force_term;
                      B_force[grad_l]      += pre_term * 2.0f * aj * B_force_term;
                      C_force[grad_l][tid] += pre_term * (2.0f * ac_val_dens_sh[j].x * C_force_term - del_pg * F_mT[0]);
                    }
                  }
                }
                C_force[0][tid] *= valid_thread * prefactor_mo;
                C_force[1][tid] *= valid_thread * prefactor_mo;
                C_force[2][tid] *= valid_thread * prefactor_mo;
              }
//------------------------------------------END TERM-TYPE DEPENDENT PART (SS-P)----------------------------------------------
