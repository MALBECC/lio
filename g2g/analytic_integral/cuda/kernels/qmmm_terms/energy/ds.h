{
  scalar_type F_mU[3];
  {
    scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
    // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
    lio_gamma<scalar_type,2>(F_mU,U);
  }
  {
    //START INDEX i1=0, CENTER 1
    {
      scalar_type p1s_0 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
      scalar_type p1s_1 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
      //START INDEX i2=0, CENTER 1
      {
        scalar_type d12s_0 = PmA[0] * p1s_0 - PmC[0] * p1s_1;
        scalar_type norm2 = 1.0f;
        d12s_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
        my_fock[0] += (double)( preterm * clatom_charge_sh[j] * d12s_0 );
      }
    }
    //START INDEX i1=1, CENTER 1
    {
      scalar_type p1s_0 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
      scalar_type p1s_1 = PmA[1] * F_mU[1] - PmC[1] * F_mU[2];
      //START INDEX i2=0, CENTER 1
      {
        scalar_type d12s_0 = PmA[0] * p1s_0 - PmC[0] * p1s_1;
        scalar_type norm2 = 1.0f;
        scalar_type preterm = norm2;
        my_fock[1] += (double)( preterm * clatom_charge_sh[j] * d12s_0 );
      }
      //START INDEX i2=1, CENTER 1
      {
        scalar_type d12s_0 = PmA[1] * p1s_0 - PmC[1] * p1s_1;
        scalar_type norm2 = 1.0f;
        d12s_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
        my_fock[2] += (double)( preterm * clatom_charge_sh[j] * d12s_0 );
      }
    }
    //START INDEX i1=2, CENTER 1
    {
      scalar_type p1s_0 = PmA[2] * F_mU[0] - PmC[2] * F_mU[1];
      scalar_type p1s_1 = PmA[2] * F_mU[1] - PmC[2] * F_mU[2];
      //START INDEX i2=0, CENTER 1
      {
        scalar_type d12s_0 = PmA[0] * p1s_0 - PmC[0] * p1s_1;
        scalar_type norm2 = 1.0f;
        scalar_type preterm = norm2;
        my_fock[3] += (double)( preterm * clatom_charge_sh[j] * d12s_0 );
      }
      //START INDEX i2=1, CENTER 1
      {
        scalar_type d12s_0 = PmA[1] * p1s_0 - PmC[1] * p1s_1;
        scalar_type norm2 = 1.0f;
        scalar_type preterm = norm2;
        my_fock[4] += (double)( preterm * clatom_charge_sh[j] * d12s_0 );
      }
      //START INDEX i2=2, CENTER 1
      {
        scalar_type d12s_0 = PmA[2] * p1s_0 - PmC[2] * p1s_1;
        scalar_type norm2 = 1.0f;
        d12s_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
        norm2 = G2G::gpu_normalization_factor;
        scalar_type preterm = norm2;
        my_fock[5] += (double)( preterm * clatom_charge_sh[j] * d12s_0 );
      }
    }
  }
}
