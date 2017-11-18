{
  scalar_type F_mT[1];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type, 0>(F_mT, T);
  }
  {
#ifdef FOCK_CALC
    my_fock[0] += (double)(fit_dens_sh[j + 0] * prefactor_dens * F_mT[0]);
#else
    rc_sh[0][tid] += (double)(dens[0] * prefactor_dens * F_mT[0]);
#endif
  }
}
