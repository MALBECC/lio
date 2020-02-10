{
  scalar_type F_mT[5];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type, 4>(F_mT, T);
  }
  {
    // START INDEX i1=0, CENTER 1
    {
      scalar_type p1ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
      scalar_type p1ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
      scalar_type p1ss_3 = PmA[0] * F_mT[3] + WmP[0] * F_mT[4];
      scalar_type p1ss_0 = PmA[0] * F_mT[0] + WmP[0] * F_mT[1];
      // START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type sp2s_1 = PmB[0] * F_mT[1] + WmP[0] * F_mT[2];
        scalar_type p1p2s_2 = PmB[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type sp2s_2 = PmB[0] * F_mT[2] + WmP[0] * F_mT[3];
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        p1p2s_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[0];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[0];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[0];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
      // START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type sp2s_1 = PmB[1] * F_mT[1] + WmP[1] * F_mT[2];
        scalar_type p1p2s_2 = PmB[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type sp2s_2 = PmB[1] * F_mT[2] + WmP[1] * F_mT[3];
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[1];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[1];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[1];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
      // START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type sp2s_1 = PmB[2] * F_mT[1] + WmP[2] * F_mT[2];
        scalar_type p1p2s_2 = PmB[2] * p1ss_2 + WmP[2] * p1ss_3;
        scalar_type sp2s_2 = PmB[2] * F_mT[2] + WmP[2] * F_mT[3];
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[2];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[2];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[2];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
    }
    // START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      scalar_type p1ss_2 = PmA[1] * F_mT[2] + WmP[1] * F_mT[3];
      scalar_type p1ss_3 = PmA[1] * F_mT[3] + WmP[1] * F_mT[4];
      scalar_type p1ss_0 = PmA[1] * F_mT[0] + WmP[1] * F_mT[1];
      // START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type p1p2s_2 = PmB[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type sp2s_1 = PmB[0] * F_mT[1] + WmP[0] * F_mT[2];
        scalar_type sp2s_2 = PmB[0] * F_mT[2] + WmP[0] * F_mT[3];
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[3];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[3];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[3];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
      // START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type p1p2s_2 = PmB[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type sp2s_1 = PmB[1] * F_mT[1] + WmP[1] * F_mT[2];
        scalar_type sp2s_2 = PmB[1] * F_mT[2] + WmP[1] * F_mT[3];
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        p1p2s_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[4];
          // START INDEX igrad=0
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[4];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[4];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
      // START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type p1p2s_2 = PmB[2] * p1ss_2 + WmP[2] * p1ss_3;
        scalar_type sp2s_1 = PmB[2] * F_mT[1] + WmP[2] * F_mT[2];
        scalar_type sp2s_2 = PmB[2] * F_mT[2] + WmP[2] * F_mT[3];
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[5];
          // START INDEX igrad=0
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[5];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[5];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
    }
    // START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      scalar_type p1ss_2 = PmA[2] * F_mT[2] + WmP[2] * F_mT[3];
      scalar_type p1ss_3 = PmA[2] * F_mT[3] + WmP[2] * F_mT[4];
      scalar_type p1ss_0 = PmA[2] * F_mT[0] + WmP[2] * F_mT[1];
      // START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type p1p2s_2 = PmB[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type sp2s_1 = PmB[0] * F_mT[1] + WmP[0] * F_mT[2];
        scalar_type sp2s_2 = PmB[0] * F_mT[2] + WmP[0] * F_mT[3];
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[6];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[6];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[6];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
      // START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type p1p2s_2 = PmB[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type sp2s_1 = PmB[1] * F_mT[1] + WmP[1] * F_mT[2];
        scalar_type sp2s_2 = PmB[1] * F_mT[2] + WmP[1] * F_mT[3];
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[7];
          // START INDEX igrad=0
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[7];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[7];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            B_force_term -= p1sp3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
      // START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type p1p2s_2 = PmB[2] * p1ss_2 + WmP[2] * p1ss_3;
        scalar_type sp2s_1 = PmB[2] * F_mT[1] + WmP[2] * F_mT[2];
        scalar_type sp2s_2 = PmB[2] * F_mT[2] + WmP[2] * F_mT[3];
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        p1p2s_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[0] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[0] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[8];
          // START INDEX igrad=0
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[1] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[1] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          scalar_type preterm = fit_dens_sh[j + 1] * prefactor_dens * dens[8];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term =
                inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            C_force_term -= p1p2s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_0 = WmQ[2] * sp2s_1;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_0 = WmQ[2] * p1ss_1;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_0 += inv_two_zeta_eta * F_mT[1];
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_0 += inv_two_zeta_eta * F_mT[1];
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          scalar_type preterm = fit_dens_sh[j + 2] * prefactor_dens * dens[8];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = WmQ[0] * p1p2p3_1;
            scalar_type A_force_term = WmP[0] * p1p2p3_1;
            scalar_type B_force_term = PmB[0] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[0] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = WmQ[1] * p1p2p3_1;
            scalar_type A_force_term = WmP[1] * p1p2p3_1;
            scalar_type B_force_term = PmB[1] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[1] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = inv_two_zeta_eta * sp2p3_1;
            C_force_term += inv_two_zeta_eta * p1sp3_1;
            C_force_term += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            C_force_term += WmQ[2] * p1p2p3_1;
            scalar_type A_force_term =
                inv_two_zeta * (sp2p3_0 - rho_zeta * sp2p3_1);
            A_force_term += inv_two_zeta * (p1sp3_0 - rho_zeta * p1sp3_1);
            A_force_term += inv_two_zeta_eta * p1p2s_1;
            A_force_term += WmP[2] * p1p2p3_1;
            scalar_type B_force_term = PmB[2] * p1p2p3_0 + A_force_term;
            A_force_term += PmA[2] * p1p2p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            C_force_term *= 2.0 * ac_val_dens_sh[j].x;
            A_force_term -= sp2p3_0;
            B_force_term -= p1sp3_0;
            C_force_term -= p1p2s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
    }
  }
}
