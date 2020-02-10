{
  scalar_type F_mU[2];
  {
    scalar_type U =
        (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
    // TODO (maybe): test out storing F(m,U) values in texture and doing a
    // texture fetch here rather than the function calculation
    lio_gamma<scalar_type, 1>(F_mU, U);
  }
  {
    scalar_type preterm = clatom_charge_sh[j] * dens[0];
    // START INDEX igrad=0
    {
      scalar_type C_force_term = PmC[0] * F_mU[1];
      scalar_type A_force_term = -C_force_term;
      scalar_type B_force_term = PmB[0] * F_mU[0] + A_force_term;
      A_force_term += PmA[0] * F_mU[0];
      A_force_term *= 2.0 * ai;
      B_force_term *= 2.0 * aj;
      A_force[0] += preterm * A_force_term;
      B_force[0] += preterm * B_force_term;
      C_force[0][tid] += preterm * C_force_term;
    }
    // START INDEX igrad=1
    {
      scalar_type C_force_term = PmC[1] * F_mU[1];
      scalar_type A_force_term = -C_force_term;
      scalar_type B_force_term = PmB[1] * F_mU[0] + A_force_term;
      A_force_term += PmA[1] * F_mU[0];
      A_force_term *= 2.0 * ai;
      B_force_term *= 2.0 * aj;
      A_force[1] += preterm * A_force_term;
      B_force[1] += preterm * B_force_term;
      C_force[1][tid] += preterm * C_force_term;
    }
    // START INDEX igrad=2
    {
      scalar_type C_force_term = PmC[2] * F_mU[1];
      scalar_type A_force_term = -C_force_term;
      scalar_type B_force_term = PmB[2] * F_mU[0] + A_force_term;
      A_force_term += PmA[2] * F_mU[0];
      A_force_term *= 2.0 * ai;
      B_force_term *= 2.0 * aj;
      A_force[2] += preterm * A_force_term;
      B_force[2] += preterm * B_force_term;
      C_force[2][tid] += preterm * C_force_term;
    }
  }
}
