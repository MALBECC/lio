{
  scalar_type F_mT[2];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type,1>(F_mT,T);
  }
  {
    //START INDEX i1=0, CENTER 3
    {
      scalar_type ssp1_0 = WmQ[0] * F_mT[1];
      my_fock[0] += (double)( fit_dens_sh[j+0] * prefactor_dens *  ssp1_0 );
    }
    //START INDEX i1=1, CENTER 3
    {
      scalar_type ssp1_0 = WmQ[1] * F_mT[1];
      my_fock[0] += (double)( fit_dens_sh[j+1] * prefactor_dens *  ssp1_0 );
    }
    //START INDEX i1=2, CENTER 3
    {
      scalar_type ssp1_0 = WmQ[2] * F_mT[1];
      my_fock[0] += (double)( fit_dens_sh[j+2] * prefactor_dens *  ssp1_0 );
    }
  }
}
