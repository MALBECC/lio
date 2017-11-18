{
  scalar_type F_mT[3];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type, 2>(F_mT, T);
  }
  {
    // START INDEX i1=0, CENTER 1
    {
      scalar_type p1ss_0 = PmA[0] * F_mT[0] + WmP[0] * F_mT[1];
      scalar_type p1ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
      scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[0];
      // START INDEX igrad=0
      {
        scalar_type C_force_term = inv_two_zeta_eta * F_mT[1];
        C_force_term += WmQ[0] * p1ss_1;
        scalar_type A_force_term =
            inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        A_force_term += WmP[0] * p1ss_1;
        scalar_type B_force_term = PmB[0] * p1ss_0 + A_force_term;
        A_force_term += PmA[0] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force_term -= F_mT[0];
        A_force[0] += preterm * A_force_term;
        B_force[0] += preterm * B_force_term;
        C_force[0][tid] += preterm * C_force_term;
      }
      // START INDEX igrad=1
      {
        scalar_type C_force_term = WmQ[1] * p1ss_1;
        scalar_type A_force_term = WmP[1] * p1ss_1;
        scalar_type B_force_term = PmB[1] * p1ss_0 + A_force_term;
        A_force_term += PmA[1] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force[1] += preterm * A_force_term;
        B_force[1] += preterm * B_force_term;
        C_force[1][tid] += preterm * C_force_term;
      }
      // START INDEX igrad=2
      {
        scalar_type C_force_term = WmQ[2] * p1ss_1;
        scalar_type A_force_term = WmP[2] * p1ss_1;
        scalar_type B_force_term = PmB[2] * p1ss_0 + A_force_term;
        A_force_term += PmA[2] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force[2] += preterm * A_force_term;
        B_force[2] += preterm * B_force_term;
        C_force[2][tid] += preterm * C_force_term;
      }
    }
    // START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_0 = PmA[1] * F_mT[0] + WmP[1] * F_mT[1];
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[1];
      // START INDEX igrad=0
      {
        scalar_type C_force_term = WmQ[0] * p1ss_1;
        scalar_type A_force_term = WmP[0] * p1ss_1;
        scalar_type B_force_term = PmB[0] * p1ss_0 + A_force_term;
        A_force_term += PmA[0] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force[0] += preterm * A_force_term;
        B_force[0] += preterm * B_force_term;
        C_force[0][tid] += preterm * C_force_term;
      }
      // START INDEX igrad=1
      {
        scalar_type C_force_term = inv_two_zeta_eta * F_mT[1];
        C_force_term += WmQ[1] * p1ss_1;
        scalar_type A_force_term =
            inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        A_force_term += WmP[1] * p1ss_1;
        scalar_type B_force_term = PmB[1] * p1ss_0 + A_force_term;
        A_force_term += PmA[1] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force_term -= F_mT[0];
        A_force[1] += preterm * A_force_term;
        B_force[1] += preterm * B_force_term;
        C_force[1][tid] += preterm * C_force_term;
      }
      // START INDEX igrad=2
      {
        scalar_type C_force_term = WmQ[2] * p1ss_1;
        scalar_type A_force_term = WmP[2] * p1ss_1;
        scalar_type B_force_term = PmB[2] * p1ss_0 + A_force_term;
        A_force_term += PmA[2] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force[2] += preterm * A_force_term;
        B_force[2] += preterm * B_force_term;
        C_force[2][tid] += preterm * C_force_term;
      }
    }
    // START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_0 = PmA[2] * F_mT[0] + WmP[2] * F_mT[1];
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      scalar_type preterm = fit_dens_sh[j + 0] * prefactor_dens * dens[2];
      // START INDEX igrad=0
      {
        scalar_type C_force_term = WmQ[0] * p1ss_1;
        scalar_type A_force_term = WmP[0] * p1ss_1;
        scalar_type B_force_term = PmB[0] * p1ss_0 + A_force_term;
        A_force_term += PmA[0] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force[0] += preterm * A_force_term;
        B_force[0] += preterm * B_force_term;
        C_force[0][tid] += preterm * C_force_term;
      }
      // START INDEX igrad=1
      {
        scalar_type C_force_term = WmQ[1] * p1ss_1;
        scalar_type A_force_term = WmP[1] * p1ss_1;
        scalar_type B_force_term = PmB[1] * p1ss_0 + A_force_term;
        A_force_term += PmA[1] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force[1] += preterm * A_force_term;
        B_force[1] += preterm * B_force_term;
        C_force[1][tid] += preterm * C_force_term;
      }
      // START INDEX igrad=2
      {
        scalar_type C_force_term = inv_two_zeta_eta * F_mT[1];
        C_force_term += WmQ[2] * p1ss_1;
        scalar_type A_force_term =
            inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        A_force_term += WmP[2] * p1ss_1;
        scalar_type B_force_term = PmB[2] * p1ss_0 + A_force_term;
        A_force_term += PmA[2] * p1ss_0;
        A_force_term *= 2.0 * ai;
        B_force_term *= 2.0 * aj;
        C_force_term *= 2.0 * ac_val_dens_sh[j].x;
        A_force_term -= F_mT[0];
        A_force[2] += preterm * A_force_term;
        B_force[2] += preterm * B_force_term;
        C_force[2][tid] += preterm * C_force_term;
      }
    }
  }
}
