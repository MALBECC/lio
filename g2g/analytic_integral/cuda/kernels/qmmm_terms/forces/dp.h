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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[0];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            B_force_term -= d12s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = PmC[1] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = PmC[2] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[1];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = PmC[2] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[2];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = PmC[1] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type p2p3_0 = PmB[0] * p2s_0 - PmC[0] * p2s_1;
          scalar_type p2p3_1 = PmB[0] * p2s_1 - PmC[0] * p2s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          p2p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p2p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[3];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            B_force_term -= d12s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = PmC[2] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[4];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            B_force_term -= d12s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = PmC[2] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[5];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type p2p3_0 = PmB[0] * p2s_0 - PmC[0] * p2s_1;
          scalar_type p2p3_1 = PmB[0] * p2s_1 - PmC[0] * p2s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[6];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = PmC[2] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[7];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = PmC[0] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            B_force_term -= d12s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = PmC[2] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[8];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = PmC[0] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type p2p3_0 = PmB[0] * p2s_0 - PmC[0] * p2s_1;
          scalar_type p2p3_1 = PmB[0] * p2s_1 - PmC[0] * p2s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          p2p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p2p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[9];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            B_force_term -= d12s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = PmC[1] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[10];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[11];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = PmC[1] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            B_force_term -= d12s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type p2p3_0 = PmB[0] * p2s_0 - PmC[0] * p2s_1;
          scalar_type p2p3_1 = PmB[0] * p2s_1 - PmC[0] * p2s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[12];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          d12p3_0 += inv_two_zeta * (p1s_0 - p1s_1);
          d12p3_1 += inv_two_zeta * (p1s_1 - p1s_2);
          p2p3_0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
          p2p3_1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[13];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = PmC[0] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            B_force_term -= d12s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[14];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = PmC[0] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p1p3_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            B_force_term -= d12s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type p2p3_0 = PmB[0] * p2s_0 - PmC[0] * p2s_1;
          scalar_type p2p3_1 = PmB[0] * p2s_1 - PmC[0] * p2s_2;
          scalar_type p1p3_0 = PmB[0] * p1s_0 - PmC[0] * p1s_1;
          scalar_type p1p3_1 = PmB[0] * p1s_1 - PmC[0] * p1s_2;
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[15];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[0] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = PmC[1] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[16];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = PmC[0] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[1] * d12p3_1;
            scalar_type A_force_term = d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            B_force_term -= d12s_0;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
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
          scalar_type preterm = norm2;
          preterm *= clatom_charge_sh[j] * dens[17];
          // START INDEX igrad=0
          {
            scalar_type C_force_term = PmC[0] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[0] * d12p3_0 + A_force_term;
            A_force_term += PmA[0] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[0] += preterm * A_force_term;
            B_force[0] += preterm * B_force_term;
            C_force[0][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=1
          {
            scalar_type C_force_term = PmC[1] * d12p3_1;
            scalar_type A_force_term = -C_force_term;
            scalar_type B_force_term = PmB[1] * d12p3_0 + A_force_term;
            A_force_term += PmA[1] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force[1] += preterm * A_force_term;
            B_force[1] += preterm * B_force_term;
            C_force[1][tid] += preterm * C_force_term;
          }
          // START INDEX igrad=2
          {
            scalar_type C_force_term = p2p3_1;
            C_force_term += p1p3_1;
            C_force_term += d12s_1;
            C_force_term *= inv_two_zeta;
            C_force_term += PmC[2] * d12p3_1;
            scalar_type A_force_term = p2p3_0;
            A_force_term += p1p3_0;
            A_force_term += d12s_0;
            A_force_term *= inv_two_zeta;
            A_force_term -= C_force_term;
            scalar_type B_force_term = PmB[2] * d12p3_0 + A_force_term;
            A_force_term += PmA[2] * d12p3_0;
            A_force_term *= 2.0 * ai;
            B_force_term *= 2.0 * aj;
            A_force_term -= p2p3_0;
            A_force_term -= p1p3_0;
            B_force_term -= d12s_0;
            A_force[2] += preterm * A_force_term;
            B_force[2] += preterm * B_force_term;
            C_force[2][tid] += preterm * C_force_term;
          }
        }
      }
    }
  }
}
