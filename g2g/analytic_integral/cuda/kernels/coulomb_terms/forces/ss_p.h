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
    //START INDEX i1=0, CENTER 3
    {
      scalar_type ssp1_0 = WmQ[0] * F_mT[1];
      scalar_type ssp1_1 = WmQ[0] * F_mT[2];
      scalar_type preterm = fit_dens_sh[j+0] * prefactor_dens * dens[0];
      //START INDEX igrad=0
      {
        scalar_type C_force_term = inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        C_force_term += WmQ[0] * ssp1_1;
        scalar_type A_force_term = inv_two_zeta_eta * F_mT[1];
        A_force_term += WmP[0] * ssp1_1;
        scalar_type B_force_term = PmB[0] * ssp1_0 + A_force_term;
        A_force_term += PmA[0] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        C_force_term -= F_mT[0];
        A_force[0]      += preterm * A_force_term;
        B_force[0]      += preterm * B_force_term;
        C_force[0][tid] += preterm * C_force_term;
      }
      //START INDEX igrad=1
      {
        scalar_type C_force_term = WmQ[1] * ssp1_1;
        scalar_type A_force_term = WmP[1] * ssp1_1;
        scalar_type B_force_term = PmB[1] * ssp1_0 + A_force_term;
        A_force_term += PmA[1] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        A_force[1]      += preterm * A_force_term;
        B_force[1]      += preterm * B_force_term;
        C_force[1][tid] += preterm * C_force_term;
      }
      //START INDEX igrad=2
      {
        scalar_type C_force_term = WmQ[2] * ssp1_1;
        scalar_type A_force_term = WmP[2] * ssp1_1;
        scalar_type B_force_term = PmB[2] * ssp1_0 + A_force_term;
        A_force_term += PmA[2] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        A_force[2]      += preterm * A_force_term;
        B_force[2]      += preterm * B_force_term;
        C_force[2][tid] += preterm * C_force_term;
      }
    }
    //START INDEX i1=1, CENTER 3
    {
      scalar_type ssp1_0 = WmQ[1] * F_mT[1];
      scalar_type ssp1_1 = WmQ[1] * F_mT[2];
      scalar_type preterm = fit_dens_sh[j+1] * prefactor_dens * dens[0];
      //START INDEX igrad=0
      {
        scalar_type C_force_term = WmQ[0] * ssp1_1;
        scalar_type A_force_term = WmP[0] * ssp1_1;
        scalar_type B_force_term = PmB[0] * ssp1_0 + A_force_term;
        A_force_term += PmA[0] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        A_force[0]      += preterm * A_force_term;
        B_force[0]      += preterm * B_force_term;
        C_force[0][tid] += preterm * C_force_term;
      }
      //START INDEX igrad=1
      {
        scalar_type C_force_term = inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        C_force_term += WmQ[1] * ssp1_1;
        scalar_type A_force_term = inv_two_zeta_eta * F_mT[1];
        A_force_term += WmP[1] * ssp1_1;
        scalar_type B_force_term = PmB[1] * ssp1_0 + A_force_term;
        A_force_term += PmA[1] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        C_force_term -= F_mT[0];
        A_force[1]      += preterm * A_force_term;
        B_force[1]      += preterm * B_force_term;
        C_force[1][tid] += preterm * C_force_term;
      }
      //START INDEX igrad=2
      {
        scalar_type C_force_term = WmQ[2] * ssp1_1;
        scalar_type A_force_term = WmP[2] * ssp1_1;
        scalar_type B_force_term = PmB[2] * ssp1_0 + A_force_term;
        A_force_term += PmA[2] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        A_force[2]      += preterm * A_force_term;
        B_force[2]      += preterm * B_force_term;
        C_force[2][tid] += preterm * C_force_term;
      }
    }
    //START INDEX i1=2, CENTER 3
    {
      scalar_type ssp1_0 = WmQ[2] * F_mT[1];
      scalar_type ssp1_1 = WmQ[2] * F_mT[2];
      scalar_type preterm = fit_dens_sh[j+2] * prefactor_dens * dens[0];
      //START INDEX igrad=0
      {
        scalar_type C_force_term = WmQ[0] * ssp1_1;
        scalar_type A_force_term = WmP[0] * ssp1_1;
        scalar_type B_force_term = PmB[0] * ssp1_0 + A_force_term;
        A_force_term += PmA[0] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        A_force[0]      += preterm * A_force_term;
        B_force[0]      += preterm * B_force_term;
        C_force[0][tid] += preterm * C_force_term;
      }
      //START INDEX igrad=1
      {
        scalar_type C_force_term = WmQ[1] * ssp1_1;
        scalar_type A_force_term = WmP[1] * ssp1_1;
        scalar_type B_force_term = PmB[1] * ssp1_0 + A_force_term;
        A_force_term += PmA[1] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        A_force[1]      += preterm * A_force_term;
        B_force[1]      += preterm * B_force_term;
        C_force[1][tid] += preterm * C_force_term;
      }
      //START INDEX igrad=2
      {
        scalar_type C_force_term = inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        C_force_term += WmQ[2] * ssp1_1;
        scalar_type A_force_term = inv_two_zeta_eta * F_mT[1];
        A_force_term += WmP[2] * ssp1_1;
        scalar_type B_force_term = PmB[2] * ssp1_0 + A_force_term;
        A_force_term += PmA[2] * ssp1_0;
        A_force_term *= 2.0f * ai;
        B_force_term *= 2.0f * aj;
        C_force_term *= 2.0f * ac_val_dens_sh[j].x;
        C_force_term -= F_mT[0];
        A_force[2]      += preterm * A_force_term;
        B_force[2]      += preterm * B_force_term;
        C_force[2][tid] += preterm * C_force_term;
      }
    }
  }
}
