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
    //START INDEX i1=0, CENTER 3
    {
      scalar_type ssp1_1 = WmQ[0] * F_mT[2];
      scalar_type ssp1_2 = WmQ[0] * F_mT[3];
      scalar_type ssp1_0 = WmQ[0] * F_mT[1];
      //START INDEX i2=0, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[0] * ssp1_1;
        scalar_type ssd12_1 = WmQ[0] * ssp1_2;
        scalar_type ssp2_0 = WmQ[0] * F_mT[1];
        scalar_type ssp2_1 = WmQ[0] * F_mT[2];
        scalar_type norm2 = 1.0;
        ssd12_0 += inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        ssd12_1 += inv_two_eta * (F_mT[1] - rho_eta * F_mT[2]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
        preterm *= fit_dens_sh[j+0] * prefactor_dens * dens[0];
        //START INDEX igrad=0
        {
          scalar_type C_force_term = inv_two_eta * (ssp2_0 - rho_eta * ssp2_1);
          C_force_term += inv_two_eta * (ssp1_0 - rho_eta * ssp1_1);
          C_force_term += WmQ[0] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp2_1;
          A_force_term += inv_two_zeta_eta * ssp1_1;
          A_force_term += WmP[0] * ssd12_1;
          scalar_type B_force_term = PmB[0] * ssd12_0 + A_force_term;
          A_force_term += PmA[0] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp2_0;
          C_force_term -= ssp1_0;
          A_force[0]      += preterm * A_force_term;
          B_force[0]      += preterm * B_force_term;
          C_force[0][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=1
        {
          scalar_type C_force_term = WmQ[1] * ssd12_1;
          scalar_type A_force_term = WmP[1] * ssd12_1;
          scalar_type B_force_term = PmB[1] * ssd12_0 + A_force_term;
          A_force_term += PmA[1] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[1]      += preterm * A_force_term;
          B_force[1]      += preterm * B_force_term;
          C_force[1][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=2
        {
          scalar_type C_force_term = WmQ[2] * ssd12_1;
          scalar_type A_force_term = WmP[2] * ssd12_1;
          scalar_type B_force_term = PmB[2] * ssd12_0 + A_force_term;
          A_force_term += PmA[2] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[2]      += preterm * A_force_term;
          B_force[2]      += preterm * B_force_term;
          C_force[2][tid] += preterm * C_force_term;
        }
      }
    }
    //START INDEX i1=1, CENTER 3
    {
      scalar_type ssp1_1 = WmQ[1] * F_mT[2];
      scalar_type ssp1_2 = WmQ[1] * F_mT[3];
      scalar_type ssp1_0 = WmQ[1] * F_mT[1];
      //START INDEX i2=0, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[0] * ssp1_1;
        scalar_type ssd12_1 = WmQ[0] * ssp1_2;
        scalar_type ssp2_0 = WmQ[0] * F_mT[1];
        scalar_type ssp2_1 = WmQ[0] * F_mT[2];
        scalar_type norm2 = 1.0;
        scalar_type preterm = norm2;
        preterm *= fit_dens_sh[j+1] * prefactor_dens * dens[0];
        //START INDEX igrad=0
        {
          scalar_type C_force_term = inv_two_eta * (ssp1_0 - rho_eta * ssp1_1);
          C_force_term += WmQ[0] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp1_1;
          A_force_term += WmP[0] * ssd12_1;
          scalar_type B_force_term = PmB[0] * ssd12_0 + A_force_term;
          A_force_term += PmA[0] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp1_0;
          A_force[0]      += preterm * A_force_term;
          B_force[0]      += preterm * B_force_term;
          C_force[0][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=1
        {
          scalar_type C_force_term = inv_two_eta * (ssp2_0 - rho_eta * ssp2_1);
          C_force_term += WmQ[1] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp2_1;
          A_force_term += WmP[1] * ssd12_1;
          scalar_type B_force_term = PmB[1] * ssd12_0 + A_force_term;
          A_force_term += PmA[1] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp2_0;
          A_force[1]      += preterm * A_force_term;
          B_force[1]      += preterm * B_force_term;
          C_force[1][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=2
        {
          scalar_type C_force_term = WmQ[2] * ssd12_1;
          scalar_type A_force_term = WmP[2] * ssd12_1;
          scalar_type B_force_term = PmB[2] * ssd12_0 + A_force_term;
          A_force_term += PmA[2] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[2]      += preterm * A_force_term;
          B_force[2]      += preterm * B_force_term;
          C_force[2][tid] += preterm * C_force_term;
        }
      }
      //START INDEX i2=1, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[1] * ssp1_1;
        scalar_type ssd12_1 = WmQ[1] * ssp1_2;
        scalar_type ssp2_0 = WmQ[1] * F_mT[1];
        scalar_type ssp2_1 = WmQ[1] * F_mT[2];
        scalar_type norm2 = 1.0;
        ssd12_0 += inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        ssd12_1 += inv_two_eta * (F_mT[1] - rho_eta * F_mT[2]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
        preterm *= fit_dens_sh[j+2] * prefactor_dens * dens[0];
        //START INDEX igrad=0
        {
          scalar_type C_force_term = WmQ[0] * ssd12_1;
          scalar_type A_force_term = WmP[0] * ssd12_1;
          scalar_type B_force_term = PmB[0] * ssd12_0 + A_force_term;
          A_force_term += PmA[0] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[0]      += preterm * A_force_term;
          B_force[0]      += preterm * B_force_term;
          C_force[0][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=1
        {
          scalar_type C_force_term = inv_two_eta * (ssp2_0 - rho_eta * ssp2_1);
          C_force_term += inv_two_eta * (ssp1_0 - rho_eta * ssp1_1);
          C_force_term += WmQ[1] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp2_1;
          A_force_term += inv_two_zeta_eta * ssp1_1;
          A_force_term += WmP[1] * ssd12_1;
          scalar_type B_force_term = PmB[1] * ssd12_0 + A_force_term;
          A_force_term += PmA[1] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp2_0;
          C_force_term -= ssp1_0;
          A_force[1]      += preterm * A_force_term;
          B_force[1]      += preterm * B_force_term;
          C_force[1][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=2
        {
          scalar_type C_force_term = WmQ[2] * ssd12_1;
          scalar_type A_force_term = WmP[2] * ssd12_1;
          scalar_type B_force_term = PmB[2] * ssd12_0 + A_force_term;
          A_force_term += PmA[2] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[2]      += preterm * A_force_term;
          B_force[2]      += preterm * B_force_term;
          C_force[2][tid] += preterm * C_force_term;
        }
      }
    }
    //START INDEX i1=2, CENTER 3
    {
      scalar_type ssp1_1 = WmQ[2] * F_mT[2];
      scalar_type ssp1_2 = WmQ[2] * F_mT[3];
      scalar_type ssp1_0 = WmQ[2] * F_mT[1];
      //START INDEX i2=0, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[0] * ssp1_1;
        scalar_type ssd12_1 = WmQ[0] * ssp1_2;
        scalar_type ssp2_0 = WmQ[0] * F_mT[1];
        scalar_type ssp2_1 = WmQ[0] * F_mT[2];
        scalar_type norm2 = 1.0;
        scalar_type preterm = norm2;
        preterm *= fit_dens_sh[j+3] * prefactor_dens * dens[0];
        //START INDEX igrad=0
        {
          scalar_type C_force_term = inv_two_eta * (ssp1_0 - rho_eta * ssp1_1);
          C_force_term += WmQ[0] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp1_1;
          A_force_term += WmP[0] * ssd12_1;
          scalar_type B_force_term = PmB[0] * ssd12_0 + A_force_term;
          A_force_term += PmA[0] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp1_0;
          A_force[0]      += preterm * A_force_term;
          B_force[0]      += preterm * B_force_term;
          C_force[0][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=1
        {
          scalar_type C_force_term = WmQ[1] * ssd12_1;
          scalar_type A_force_term = WmP[1] * ssd12_1;
          scalar_type B_force_term = PmB[1] * ssd12_0 + A_force_term;
          A_force_term += PmA[1] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[1]      += preterm * A_force_term;
          B_force[1]      += preterm * B_force_term;
          C_force[1][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=2
        {
          scalar_type C_force_term = inv_two_eta * (ssp2_0 - rho_eta * ssp2_1);
          C_force_term += WmQ[2] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp2_1;
          A_force_term += WmP[2] * ssd12_1;
          scalar_type B_force_term = PmB[2] * ssd12_0 + A_force_term;
          A_force_term += PmA[2] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp2_0;
          A_force[2]      += preterm * A_force_term;
          B_force[2]      += preterm * B_force_term;
          C_force[2][tid] += preterm * C_force_term;
        }
      }
      //START INDEX i2=1, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[1] * ssp1_1;
        scalar_type ssd12_1 = WmQ[1] * ssp1_2;
        scalar_type ssp2_0 = WmQ[1] * F_mT[1];
        scalar_type ssp2_1 = WmQ[1] * F_mT[2];
        scalar_type norm2 = 1.0;
        scalar_type preterm = norm2;
        preterm *= fit_dens_sh[j+4] * prefactor_dens * dens[0];
        //START INDEX igrad=0
        {
          scalar_type C_force_term = WmQ[0] * ssd12_1;
          scalar_type A_force_term = WmP[0] * ssd12_1;
          scalar_type B_force_term = PmB[0] * ssd12_0 + A_force_term;
          A_force_term += PmA[0] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[0]      += preterm * A_force_term;
          B_force[0]      += preterm * B_force_term;
          C_force[0][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=1
        {
          scalar_type C_force_term = inv_two_eta * (ssp1_0 - rho_eta * ssp1_1);
          C_force_term += WmQ[1] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp1_1;
          A_force_term += WmP[1] * ssd12_1;
          scalar_type B_force_term = PmB[1] * ssd12_0 + A_force_term;
          A_force_term += PmA[1] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp1_0;
          A_force[1]      += preterm * A_force_term;
          B_force[1]      += preterm * B_force_term;
          C_force[1][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=2
        {
          scalar_type C_force_term = inv_two_eta * (ssp2_0 - rho_eta * ssp2_1);
          C_force_term += WmQ[2] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp2_1;
          A_force_term += WmP[2] * ssd12_1;
          scalar_type B_force_term = PmB[2] * ssd12_0 + A_force_term;
          A_force_term += PmA[2] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp2_0;
          A_force[2]      += preterm * A_force_term;
          B_force[2]      += preterm * B_force_term;
          C_force[2][tid] += preterm * C_force_term;
        }
      }
      //START INDEX i2=2, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[2] * ssp1_1;
        scalar_type ssd12_1 = WmQ[2] * ssp1_2;
        scalar_type ssp2_0 = WmQ[2] * F_mT[1];
        scalar_type ssp2_1 = WmQ[2] * F_mT[2];
        scalar_type norm2 = 1.0;
        ssd12_0 += inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        ssd12_1 += inv_two_eta * (F_mT[1] - rho_eta * F_mT[2]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
        preterm *= fit_dens_sh[j+5] * prefactor_dens * dens[0];
        //START INDEX igrad=0
        {
          scalar_type C_force_term = WmQ[0] * ssd12_1;
          scalar_type A_force_term = WmP[0] * ssd12_1;
          scalar_type B_force_term = PmB[0] * ssd12_0 + A_force_term;
          A_force_term += PmA[0] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[0]      += preterm * A_force_term;
          B_force[0]      += preterm * B_force_term;
          C_force[0][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=1
        {
          scalar_type C_force_term = WmQ[1] * ssd12_1;
          scalar_type A_force_term = WmP[1] * ssd12_1;
          scalar_type B_force_term = PmB[1] * ssd12_0 + A_force_term;
          A_force_term += PmA[1] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          A_force[1]      += preterm * A_force_term;
          B_force[1]      += preterm * B_force_term;
          C_force[1][tid] += preterm * C_force_term;
        }
        //START INDEX igrad=2
        {
          scalar_type C_force_term = inv_two_eta * (ssp2_0 - rho_eta * ssp2_1);
          C_force_term += inv_two_eta * (ssp1_0 - rho_eta * ssp1_1);
          C_force_term += WmQ[2] * ssd12_1;
          scalar_type A_force_term = inv_two_zeta_eta * ssp2_1;
          A_force_term += inv_two_zeta_eta * ssp1_1;
          A_force_term += WmP[2] * ssd12_1;
          scalar_type B_force_term = PmB[2] * ssd12_0 + A_force_term;
          A_force_term += PmA[2] * ssd12_0;
          A_force_term *= 2.0 * ai;
          B_force_term *= 2.0 * aj;
          C_force_term *= 2.0 * ac_val_dens_sh[j].x;
          C_force_term -= ssp2_0;
          C_force_term -= ssp1_0;
          A_force[2]      += preterm * A_force_term;
          B_force[2]      += preterm * B_force_term;
          C_force[2][tid] += preterm * C_force_term;
        }
      }
    }
  }
}
