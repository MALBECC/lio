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
      //START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        my_fock[0] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
      //START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        my_fock[1] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
      //START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        my_fock[2] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
    }
    //START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_0 = PmA[1] * F_mT[0] + WmP[1] * F_mT[1];
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      //START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        my_fock[3] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
      //START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        my_fock[4] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
      //START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        my_fock[5] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
    }
    //START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_0 = PmA[2] * F_mT[0] + WmP[2] * F_mT[1];
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      //START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        my_fock[6] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
      //START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        my_fock[7] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
      //START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        my_fock[8] += (double)( fit_dens_sh[j+0] * prefactor_dens *  p1p2s_0 );
      }
    }
  }
}
