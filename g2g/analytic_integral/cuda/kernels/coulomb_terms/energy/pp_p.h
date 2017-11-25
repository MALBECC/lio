{
  scalar_type F_mT[4];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type, 3>(F_mT, T);
  }
  {
    // START INDEX i1=0, CENTER 1
    {
      scalar_type p1ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
      scalar_type p1ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
      // START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type sp2s_1 = PmB[0] * F_mT[1] + WmP[0] * F_mT[2];
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[0] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[0] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[0] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[0] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[0] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[0] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
      // START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type sp2s_1 = PmB[1] * F_mT[1] + WmP[1] * F_mT[2];
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
#ifdef FOCK_CALC
          my_fock[1] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[1] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[1] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[1] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[1] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[1] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
      // START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type sp2s_1 = PmB[2] * F_mT[1] + WmP[2] * F_mT[2];
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
#ifdef FOCK_CALC
          my_fock[2] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[2] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[2] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[2] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[2] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[2] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
    }
    // START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      scalar_type p1ss_2 = PmA[1] * F_mT[2] + WmP[1] * F_mT[3];
      // START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type sp2s_1 = PmB[0] * F_mT[1] + WmP[0] * F_mT[2];
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[3] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[3] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
#ifdef FOCK_CALC
          my_fock[3] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[3] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[3] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[3] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
      // START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type sp2s_1 = PmB[1] * F_mT[1] + WmP[1] * F_mT[2];
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[4] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[4] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[4] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[4] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[4] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[4] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
      // START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type sp2s_1 = PmB[2] * F_mT[1] + WmP[2] * F_mT[2];
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[5] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[5] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
#ifdef FOCK_CALC
          my_fock[5] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[5] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[5] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[5] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
    }
    // START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      scalar_type p1ss_2 = PmA[2] * F_mT[2] + WmP[2] * F_mT[3];
      // START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type sp2s_1 = PmB[0] * F_mT[1] + WmP[0] * F_mT[2];
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[6] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[6] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[6] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[6] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
#ifdef FOCK_CALC
          my_fock[6] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[6] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
      // START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type sp2s_1 = PmB[1] * F_mT[1] + WmP[1] * F_mT[2];
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[7] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[7] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[7] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[7] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
#ifdef FOCK_CALC
          my_fock[7] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[7] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
      // START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type sp2s_1 = PmB[2] * F_mT[1] + WmP[2] * F_mT[2];
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        // START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[0] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[8] +=
              (double)(fit_dens_sh[j + 0] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[0][tid] += (double)(dens[8] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[1] * p1p2s_1;
#ifdef FOCK_CALC
          my_fock[8] +=
              (double)(fit_dens_sh[j + 1] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[1][tid] += (double)(dens[8] * prefactor_dens * p1p2p3_0);
#endif
        }
        // START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_0 = WmQ[2] * p1p2s_1;
          p1p2p3_0 += inv_two_zeta_eta * sp2s_1;
          p1p2p3_0 += inv_two_zeta_eta * p1ss_1;
#ifdef FOCK_CALC
          my_fock[8] +=
              (double)(fit_dens_sh[j + 2] * prefactor_dens * p1p2p3_0);
#else
          rc_sh[2][tid] += (double)(dens[8] * prefactor_dens * p1p2p3_0);
#endif
        }
      }
    }
  }
}
