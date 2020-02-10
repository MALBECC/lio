{
  scalar_type F_mU[5];
  {
    scalar_type U =
        (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
    // TODO (maybe): test out storing F(m,U) values in texture and doing a
    // texture fetch here rather than the function calculation
    lio_gamma<scalar_type, 4>(F_mU, U);
  }
  {
    // START INDEX i1=0, CENTER 1
    {
      scalar_type p1s_0 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
      scalar_type p1s_1 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
      scalar_type p1s_2 = PmA[0] * F_mU[2] - PmC[0] * F_mU[3];
      scalar_type p1s_3 = PmA[0] * F_mU[3] - PmC[0] * F_mU[4];
      // START INDEX i2=0, CENTER 1
      {
        scalar_type d12s_0 = PmA[0] * p1s_0 - PmC[0] * p1s_1;
        scalar_type d12s_1 = PmA[0] * p1s_1 - PmC[0] * p1s_2;
        scalar_type p2s_0 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
        scalar_type p2s_1 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
        scalar_type d12s_2 = PmA[0] * p1s_2 - PmC[0] * p1s_3;
        scalar_type p2s_2 = PmA[0] * F_mU[2] - PmC[0] * F_mU[3];
        scalar_type norm2 = 1.0;
        d12s_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
        d12s_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
        d12s_2 += inv_two_zeta * (F_mU[2] - F_mU[3]);
        norm2 = G2G::gpu_normalization_factor;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3_0 = PmB[0] * d12s_0 - PmC[0] * d12s_1;
          scalar_type d12p3_1 = PmB[0] * d12s_1 - PmC[0] * d12s_2;
          scalar_type p2p3_0 = PmB[0] * p2s_0 - PmC[0] * p2s_1;
          scalar_type p2p3_1 = PmB[0] * p2s_1 - PmC[0] * p2s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          d12p3_0 += inv_two_zeta * (p2s_0 - p2s_1);
          d12p3_1 += inv_two_zeta * (p2s_1 - p2s_2);
          p1p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p1p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          p2p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p2p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[0] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3_0 = PmB[1] * d12s_0 - PmC[1] * d12s_1;
          scalar_type d12p3_1 = PmB[1] * d12s_1 - PmC[1] * d12s_2;
          scalar_type p2p3_0 = PmB[1] * p2s_0 - PmC[1] * p2s_1;
          scalar_type p2p3_1 = PmB[1] * p2s_1 - PmC[1] * p2s_2;
          scalar_type p1p3_0 = PmB[1] * p1s_0 - PmC[1] * p1s_1;
          scalar_type p1p3_1 = PmB[1] * p1s_1 - PmC[1] * p1s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[1] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[2] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3_0 = PmB[2] * d12s_0 - PmC[2] * d12s_1;
          scalar_type d12p3_1 = PmB[2] * d12s_1 - PmC[2] * d12s_2;
          scalar_type p2p3_0 = PmB[2] * p2s_0 - PmC[2] * p2s_1;
          scalar_type p2p3_1 = PmB[2] * p2s_1 - PmC[2] * p2s_2;
          scalar_type p1p3_0 = PmB[2] * p1s_0 - PmC[2] * p1s_1;
          scalar_type p1p3_1 = PmB[2] * p1s_1 - PmC[2] * p1s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[3] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[4] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=2, CENTER 2
          {
            scalar_type d12d34_0 = PmB[2] * d12p3_0 - PmC[2] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[5] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
      }
    }
    // START INDEX i1=1, CENTER 1
    {
      scalar_type p1s_0 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
      scalar_type p1s_1 = PmA[1] * F_mU[1] - PmC[1] * F_mU[2];
      scalar_type p1s_2 = PmA[1] * F_mU[2] - PmC[1] * F_mU[3];
      scalar_type p1s_3 = PmA[1] * F_mU[3] - PmC[1] * F_mU[4];
      // START INDEX i2=0, CENTER 1
      {
        scalar_type d12s_0 = PmA[0] * p1s_0 - PmC[0] * p1s_1;
        scalar_type d12s_1 = PmA[0] * p1s_1 - PmC[0] * p1s_2;
        scalar_type d12s_2 = PmA[0] * p1s_2 - PmC[0] * p1s_3;
        scalar_type p2s_0 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
        scalar_type p2s_1 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
        scalar_type p2s_2 = PmA[0] * F_mU[2] - PmC[0] * F_mU[3];
        scalar_type norm2 = 1.0;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3_0 = PmB[0] * d12s_0 - PmC[0] * d12s_1;
          scalar_type d12p3_1 = PmB[0] * d12s_1 - PmC[0] * d12s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[6] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3_0 = PmB[1] * d12s_0 - PmC[1] * d12s_1;
          scalar_type d12p3_1 = PmB[1] * d12s_1 - PmC[1] * d12s_2;
          scalar_type p1p3_0 = PmB[1] * p1s_0 - PmC[1] * p1s_1;
          scalar_type p1p3_1 = PmB[1] * p1s_1 - PmC[1] * p1s_2;
          scalar_type p2p3_0 = PmB[1] * p2s_0 - PmC[1] * p2s_1;
          scalar_type p2p3_1 = PmB[1] * p2s_1 - PmC[1] * p2s_2;
          d12p3_0 += inv_two_zeta * (p2s_0 - p2s_1);
          d12p3_1 += inv_two_zeta * (p2s_1 - p2s_2);
          p1p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p1p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[7] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[8] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3_0 = PmB[2] * d12s_0 - PmC[2] * d12s_1;
          scalar_type d12p3_1 = PmB[2] * d12s_1 - PmC[2] * d12s_2;
          scalar_type p1p3_0 = PmB[2] * p1s_0 - PmC[2] * p1s_1;
          scalar_type p1p3_1 = PmB[2] * p1s_1 - PmC[2] * p1s_2;
          scalar_type p2p3_0 = PmB[2] * p2s_0 - PmC[2] * p2s_1;
          scalar_type p2p3_1 = PmB[2] * p2s_1 - PmC[2] * p2s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[9] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[10] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=2, CENTER 2
          {
            scalar_type d12d34_0 = PmB[2] * d12p3_0 - PmC[2] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[11] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
      }
      // START INDEX i2=1, CENTER 1
      {
        scalar_type d12s_0 = PmA[1] * p1s_0 - PmC[1] * p1s_1;
        scalar_type d12s_1 = PmA[1] * p1s_1 - PmC[1] * p1s_2;
        scalar_type d12s_2 = PmA[1] * p1s_2 - PmC[1] * p1s_3;
        scalar_type p2s_0 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
        scalar_type p2s_1 = PmA[1] * F_mU[1] - PmC[1] * F_mU[2];
        scalar_type p2s_2 = PmA[1] * F_mU[2] - PmC[1] * F_mU[3];
        scalar_type norm2 = 1.0;
        d12s_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
        d12s_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
        d12s_2 += inv_two_zeta * (F_mU[2] - F_mU[3]);
        norm2 = G2G::gpu_normalization_factor;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3_0 = PmB[0] * d12s_0 - PmC[0] * d12s_1;
          scalar_type d12p3_1 = PmB[0] * d12s_1 - PmC[0] * d12s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[12] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3_0 = PmB[1] * d12s_0 - PmC[1] * d12s_1;
          scalar_type d12p3_1 = PmB[1] * d12s_1 - PmC[1] * d12s_2;
          scalar_type p2p3_0 = PmB[1] * p2s_0 - PmC[1] * p2s_1;
          scalar_type p2p3_1 = PmB[1] * p2s_1 - PmC[1] * p2s_2;
          scalar_type p1p3_0 = PmB[1] * p1s_0 - PmC[1] * p1s_1;
          scalar_type p1p3_1 = PmB[1] * p1s_1 - PmC[1] * p1s_2;
          d12p3_0 += inv_two_zeta * (p2s_0 - p2s_1);
          d12p3_1 += inv_two_zeta * (p2s_1 - p2s_2);
          p1p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p1p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          p2p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p2p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[13] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[14] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3_0 = PmB[2] * d12s_0 - PmC[2] * d12s_1;
          scalar_type d12p3_1 = PmB[2] * d12s_1 - PmC[2] * d12s_2;
          scalar_type p2p3_0 = PmB[2] * p2s_0 - PmC[2] * p2s_1;
          scalar_type p2p3_1 = PmB[2] * p2s_1 - PmC[2] * p2s_2;
          scalar_type p1p3_0 = PmB[2] * p1s_0 - PmC[2] * p1s_1;
          scalar_type p1p3_1 = PmB[2] * p1s_1 - PmC[2] * p1s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[15] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[16] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=2, CENTER 2
          {
            scalar_type d12d34_0 = PmB[2] * d12p3_0 - PmC[2] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[17] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
      }
    }
    // START INDEX i1=2, CENTER 1
    {
      scalar_type p1s_0 = PmA[2] * F_mU[0] - PmC[2] * F_mU[1];
      scalar_type p1s_1 = PmA[2] * F_mU[1] - PmC[2] * F_mU[2];
      scalar_type p1s_2 = PmA[2] * F_mU[2] - PmC[2] * F_mU[3];
      scalar_type p1s_3 = PmA[2] * F_mU[3] - PmC[2] * F_mU[4];
      // START INDEX i2=0, CENTER 1
      {
        scalar_type d12s_0 = PmA[0] * p1s_0 - PmC[0] * p1s_1;
        scalar_type d12s_1 = PmA[0] * p1s_1 - PmC[0] * p1s_2;
        scalar_type d12s_2 = PmA[0] * p1s_2 - PmC[0] * p1s_3;
        scalar_type p2s_0 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
        scalar_type p2s_1 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
        scalar_type p2s_2 = PmA[0] * F_mU[2] - PmC[0] * F_mU[3];
        scalar_type norm2 = 1.0;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3_0 = PmB[0] * d12s_0 - PmC[0] * d12s_1;
          scalar_type d12p3_1 = PmB[0] * d12s_1 - PmC[0] * d12s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[18] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3_0 = PmB[1] * d12s_0 - PmC[1] * d12s_1;
          scalar_type d12p3_1 = PmB[1] * d12s_1 - PmC[1] * d12s_2;
          scalar_type p1p3_0 = PmB[1] * p1s_0 - PmC[1] * p1s_1;
          scalar_type p1p3_1 = PmB[1] * p1s_1 - PmC[1] * p1s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[19] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[20] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3_0 = PmB[2] * d12s_0 - PmC[2] * d12s_1;
          scalar_type d12p3_1 = PmB[2] * d12s_1 - PmC[2] * d12s_2;
          scalar_type p1p3_0 = PmB[2] * p1s_0 - PmC[2] * p1s_1;
          scalar_type p1p3_1 = PmB[2] * p1s_1 - PmC[2] * p1s_2;
          scalar_type p2p3_0 = PmB[2] * p2s_0 - PmC[2] * p2s_1;
          scalar_type p2p3_1 = PmB[2] * p2s_1 - PmC[2] * p2s_2;
          d12p3_0 += inv_two_zeta * (p2s_0 - p2s_1);
          d12p3_1 += inv_two_zeta * (p2s_1 - p2s_2);
          p1p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p1p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[21] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[22] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=2, CENTER 2
          {
            scalar_type d12d34_0 = PmB[2] * d12p3_0 - PmC[2] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[23] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
      }
      // START INDEX i2=1, CENTER 1
      {
        scalar_type d12s_0 = PmA[1] * p1s_0 - PmC[1] * p1s_1;
        scalar_type d12s_1 = PmA[1] * p1s_1 - PmC[1] * p1s_2;
        scalar_type d12s_2 = PmA[1] * p1s_2 - PmC[1] * p1s_3;
        scalar_type p2s_0 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
        scalar_type p2s_1 = PmA[1] * F_mU[1] - PmC[1] * F_mU[2];
        scalar_type p2s_2 = PmA[1] * F_mU[2] - PmC[1] * F_mU[3];
        scalar_type norm2 = 1.0;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3_0 = PmB[0] * d12s_0 - PmC[0] * d12s_1;
          scalar_type d12p3_1 = PmB[0] * d12s_1 - PmC[0] * d12s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[24] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3_0 = PmB[1] * d12s_0 - PmC[1] * d12s_1;
          scalar_type d12p3_1 = PmB[1] * d12s_1 - PmC[1] * d12s_2;
          scalar_type p1p3_0 = PmB[1] * p1s_0 - PmC[1] * p1s_1;
          scalar_type p1p3_1 = PmB[1] * p1s_1 - PmC[1] * p1s_2;
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[25] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[26] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3_0 = PmB[2] * d12s_0 - PmC[2] * d12s_1;
          scalar_type d12p3_1 = PmB[2] * d12s_1 - PmC[2] * d12s_2;
          scalar_type p1p3_0 = PmB[2] * p1s_0 - PmC[2] * p1s_1;
          scalar_type p1p3_1 = PmB[2] * p1s_1 - PmC[2] * p1s_2;
          scalar_type p2p3_0 = PmB[2] * p2s_0 - PmC[2] * p2s_1;
          scalar_type p2p3_1 = PmB[2] * p2s_1 - PmC[2] * p2s_2;
          d12p3_0 += inv_two_zeta * (p2s_0 - p2s_1);
          d12p3_1 += inv_two_zeta * (p2s_1 - p2s_2);
          p1p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p1p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[27] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[28] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=2, CENTER 2
          {
            scalar_type d12d34_0 = PmB[2] * d12p3_0 - PmC[2] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[29] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
      }
      // START INDEX i2=2, CENTER 1
      {
        scalar_type d12s_0 = PmA[2] * p1s_0 - PmC[2] * p1s_1;
        scalar_type d12s_1 = PmA[2] * p1s_1 - PmC[2] * p1s_2;
        scalar_type d12s_2 = PmA[2] * p1s_2 - PmC[2] * p1s_3;
        scalar_type p2s_0 = PmA[2] * F_mU[0] - PmC[2] * F_mU[1];
        scalar_type p2s_1 = PmA[2] * F_mU[1] - PmC[2] * F_mU[2];
        scalar_type p2s_2 = PmA[2] * F_mU[2] - PmC[2] * F_mU[3];
        scalar_type norm2 = 1.0;
        d12s_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
        d12s_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
        d12s_2 += inv_two_zeta * (F_mU[2] - F_mU[3]);
        norm2 = G2G::gpu_normalization_factor;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3_0 = PmB[0] * d12s_0 - PmC[0] * d12s_1;
          scalar_type d12p3_1 = PmB[0] * d12s_1 - PmC[0] * d12s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[30] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3_0 = PmB[1] * d12s_0 - PmC[1] * d12s_1;
          scalar_type d12p3_1 = PmB[1] * d12s_1 - PmC[1] * d12s_2;
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[31] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[32] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3_0 = PmB[2] * d12s_0 - PmC[2] * d12s_1;
          scalar_type d12p3_1 = PmB[2] * d12s_1 - PmC[2] * d12s_2;
          scalar_type p2p3_0 = PmB[2] * p2s_0 - PmC[2] * p2s_1;
          scalar_type p2p3_1 = PmB[2] * p2s_1 - PmC[2] * p2s_2;
          scalar_type p1p3_0 = PmB[2] * p1s_0 - PmC[2] * p1s_1;
          scalar_type p1p3_1 = PmB[2] * p1s_1 - PmC[2] * p1s_2;
          d12p3_0 += inv_two_zeta * (p2s_0 - p2s_1);
          d12p3_1 += inv_two_zeta * (p2s_1 - p2s_2);
          p1p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p1p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          p2p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p2p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          // START INDEX i4=0, CENTER 2
          {
            scalar_type d12d34_0 = PmB[0] * d12p3_0 - PmC[0] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[33] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=1, CENTER 2
          {
            scalar_type d12d34_0 = PmB[1] * d12p3_0 - PmC[1] * d12p3_1;
            scalar_type norm4 = 1.0;
            scalar_type preterm = norm2 * norm4;
            my_fock[34] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
          // START INDEX i4=2, CENTER 2
          {
            scalar_type d12d34_0 = PmB[2] * d12p3_0 - PmC[2] * d12p3_1;
            d12d34_0 += inv_two_zeta * (p2p3_0 - p2p3_1);
            d12d34_0 += inv_two_zeta * (p1p3_0 - p1p3_1);
            scalar_type norm4 = 1.0;
            d12d34_0 += inv_two_zeta * (d12s_0 - d12s_1);
            norm4 = G2G::gpu_normalization_factor;
            scalar_type preterm = norm2 * norm4;
            my_fock[35] += (double)(preterm * clatom_charge_sh[j] * d12d34_0);
          }
        }
      }
    }
  }
}
