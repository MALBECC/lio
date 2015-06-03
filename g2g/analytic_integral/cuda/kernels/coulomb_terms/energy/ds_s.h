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
    //START INDEX i1=0, CENTER 1
    {
      scalar_type p1ss_0 = PmA[0] * F_mT[0] + WmP[0] * F_mT[1];
      scalar_type p1ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
      //START INDEX i2=0, CENTER 1
      {
        scalar_type d12ss_0 = PmA[0] * p1ss_0 + WmP[0] * p1ss_1;
        scalar_type norm2 = 1.0f;
        d12ss_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[0] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  d12ss_0 );
#else
        rc_sh[0][tid] += (double)( preterm * dens[0] * prefactor_dens *  d12ss_0 );
#endif
      }
    }
    //START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_0 = PmA[1] * F_mT[0] + WmP[1] * F_mT[1];
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      //START INDEX i2=0, CENTER 1
      {
        scalar_type d12ss_0 = PmA[0] * p1ss_0 + WmP[0] * p1ss_1;
        scalar_type norm2 = 1.0f;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[1] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  d12ss_0 );
#else
        rc_sh[0][tid] += (double)( preterm * dens[1] * prefactor_dens *  d12ss_0 );
#endif
      }
      //START INDEX i2=1, CENTER 1
      {
        scalar_type d12ss_0 = PmA[1] * p1ss_0 + WmP[1] * p1ss_1;
        scalar_type norm2 = 1.0f;
        d12ss_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[2] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  d12ss_0 );
#else
        rc_sh[0][tid] += (double)( preterm * dens[2] * prefactor_dens *  d12ss_0 );
#endif
      }
    }
    //START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_0 = PmA[2] * F_mT[0] + WmP[2] * F_mT[1];
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      //START INDEX i2=0, CENTER 1
      {
        scalar_type d12ss_0 = PmA[0] * p1ss_0 + WmP[0] * p1ss_1;
        scalar_type norm2 = 1.0f;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[3] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  d12ss_0 );
#else
        rc_sh[0][tid] += (double)( preterm * dens[3] * prefactor_dens *  d12ss_0 );
#endif
      }
      //START INDEX i2=1, CENTER 1
      {
        scalar_type d12ss_0 = PmA[1] * p1ss_0 + WmP[1] * p1ss_1;
        scalar_type norm2 = 1.0f;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[4] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  d12ss_0 );
#else
        rc_sh[0][tid] += (double)( preterm * dens[4] * prefactor_dens *  d12ss_0 );
#endif
      }
      //START INDEX i2=2, CENTER 1
      {
        scalar_type d12ss_0 = PmA[2] * p1ss_0 + WmP[2] * p1ss_1;
        scalar_type norm2 = 1.0f;
        d12ss_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
#ifdef FOCK_CALC
        my_fock[5] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  d12ss_0 );
#else
        rc_sh[0][tid] += (double)( preterm * dens[5] * prefactor_dens *  d12ss_0 );
#endif
      }
    }
  }
}
