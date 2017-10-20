{
  scalar_type F_mU[1];
  {
    scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
    // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
    lio_gamma<scalar_type,0>(F_mU,U);
  }
  {
    my_fock[0] += (double)( clatom_charge_sh[j] * F_mU[0] );
  }
}
