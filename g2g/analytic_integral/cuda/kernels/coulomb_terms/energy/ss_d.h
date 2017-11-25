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
    // START INDEX i1=0, CENTER 3
    {
      scalar_type ssp1_1 = WmQ[0] * F_mT[2];
      // START INDEX i2=0, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[0] * ssp1_1;
        scalar_type norm2 = 1.0;
        ssd12_0 += inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[0] +=
            (double)(preterm * fit_dens_sh[j + 0] * prefactor_dens * ssd12_0);
#else
        rc_sh[0][tid] += (double)(preterm * dens[0] * prefactor_dens * ssd12_0);
#endif
      }
    }
    // START INDEX i1=1, CENTER 3
    {
      scalar_type ssp1_1 = WmQ[1] * F_mT[2];
      // START INDEX i2=0, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[0] * ssp1_1;
        scalar_type norm2 = 1.0;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[0] +=
            (double)(preterm * fit_dens_sh[j + 1] * prefactor_dens * ssd12_0);
#else
        rc_sh[1][tid] += (double)(preterm * dens[0] * prefactor_dens * ssd12_0);
#endif
      }
      // START INDEX i2=1, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[1] * ssp1_1;
        scalar_type norm2 = 1.0;
        ssd12_0 += inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[0] +=
            (double)(preterm * fit_dens_sh[j + 2] * prefactor_dens * ssd12_0);
#else
        rc_sh[2][tid] += (double)(preterm * dens[0] * prefactor_dens * ssd12_0);
#endif
      }
    }
    // START INDEX i1=2, CENTER 3
    {
      scalar_type ssp1_1 = WmQ[2] * F_mT[2];
      // START INDEX i2=0, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[0] * ssp1_1;
        scalar_type norm2 = 1.0;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[0] +=
            (double)(preterm * fit_dens_sh[j + 3] * prefactor_dens * ssd12_0);
#else
        rc_sh[3][tid] += (double)(preterm * dens[0] * prefactor_dens * ssd12_0);
#endif
      }
      // START INDEX i2=1, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[1] * ssp1_1;
        scalar_type norm2 = 1.0;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[0] +=
            (double)(preterm * fit_dens_sh[j + 4] * prefactor_dens * ssd12_0);
#else
        rc_sh[4][tid] += (double)(preterm * dens[0] * prefactor_dens * ssd12_0);
#endif
      }
      // START INDEX i2=2, CENTER 3
      {
        scalar_type ssd12_0 = WmQ[2] * ssp1_1;
        scalar_type norm2 = 1.0;
        ssd12_0 += inv_two_eta * (F_mT[0] - rho_eta * F_mT[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[0] +=
            (double)(preterm * fit_dens_sh[j + 5] * prefactor_dens * ssd12_0);
#else
        rc_sh[5][tid] += (double)(preterm * dens[0] * prefactor_dens * ssd12_0);
#endif
      }
    }
  }
}
