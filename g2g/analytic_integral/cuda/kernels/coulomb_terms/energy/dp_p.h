{
  scalar_type F_mT[5];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type, 4>(F_mT, T);
  }
  {
    // START INDEX i1=0, CENTER 1
    {
      scalar_type p1ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
      scalar_type p1ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
      scalar_type p1ss_3 = PmA[0] * F_mT[3] + WmP[0] * F_mT[4];
      // START INDEX i2=0, CENTER 1
      {
        scalar_type d12ss_1 = PmA[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type d12ss_2 = PmA[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type p2ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
        scalar_type p2ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
        scalar_type norm2 = 1.0;
        d12ss_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        d12ss_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        norm2 = G2G::gpu_normalization_factor;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[0] * d12ss_1 + WmP[0] * d12ss_2;
          scalar_type p2p3s_1 = PmB[0] * p2ss_1 + WmP[0] * p2ss_2;
          scalar_type p1p3s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
          d12p3s_1 += inv_two_zeta * (p2ss_1 - rho_zeta * p2ss_2);
          p1p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          d12p3s_1 += inv_two_zeta * (p1ss_1 - rho_zeta * p1ss_2);
          p2p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[0] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[0] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[0] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[0] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[0] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[0] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[1] * d12ss_1 + WmP[1] * d12ss_2;
          scalar_type p2p3s_1 = PmB[1] * p2ss_1 + WmP[1] * p2ss_2;
          scalar_type p1p3s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[1] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[1] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[1] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[1] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[1] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[1] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[2] * d12ss_1 + WmP[2] * d12ss_2;
          scalar_type p2p3s_1 = PmB[2] * p2ss_1 + WmP[2] * p2ss_2;
          scalar_type p1p3s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[2] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[2] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[2] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[2] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[2] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[2] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
      }
    }
    // START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      scalar_type p1ss_2 = PmA[1] * F_mT[2] + WmP[1] * F_mT[3];
      scalar_type p1ss_3 = PmA[1] * F_mT[3] + WmP[1] * F_mT[4];
      // START INDEX i2=0, CENTER 1
      {
        scalar_type d12ss_1 = PmA[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type d12ss_2 = PmA[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type p2ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
        scalar_type p2ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
        scalar_type norm2 = 1.0;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[0] * d12ss_1 + WmP[0] * d12ss_2;
          scalar_type p1p3s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
          scalar_type p2p3s_1 = PmB[0] * p2ss_1 + WmP[0] * p2ss_2;
          d12p3s_1 += inv_two_zeta * (p1ss_1 - rho_zeta * p1ss_2);
          p2p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[3] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[3] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[3] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[3] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[3] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[3] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[1] * d12ss_1 + WmP[1] * d12ss_2;
          scalar_type p1p3s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
          scalar_type p2p3s_1 = PmB[1] * p2ss_1 + WmP[1] * p2ss_2;
          d12p3s_1 += inv_two_zeta * (p2ss_1 - rho_zeta * p2ss_2);
          p1p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[4] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[4] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[4] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[4] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[4] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[4] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[2] * d12ss_1 + WmP[2] * d12ss_2;
          scalar_type p1p3s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
          scalar_type p2p3s_1 = PmB[2] * p2ss_1 + WmP[2] * p2ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[5] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[5] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[5] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[5] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[5] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[5] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
      }
      // START INDEX i2=1, CENTER 1
      {
        scalar_type d12ss_1 = PmA[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type d12ss_2 = PmA[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type p2ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
        scalar_type p2ss_2 = PmA[1] * F_mT[2] + WmP[1] * F_mT[3];
        scalar_type norm2 = 1.0;
        d12ss_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        d12ss_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        norm2 = G2G::gpu_normalization_factor;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[0] * d12ss_1 + WmP[0] * d12ss_2;
          scalar_type p2p3s_1 = PmB[0] * p2ss_1 + WmP[0] * p2ss_2;
          scalar_type p1p3s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[6] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[6] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[6] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[6] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[6] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[6] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[1] * d12ss_1 + WmP[1] * d12ss_2;
          scalar_type p2p3s_1 = PmB[1] * p2ss_1 + WmP[1] * p2ss_2;
          scalar_type p1p3s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
          d12p3s_1 += inv_two_zeta * (p2ss_1 - rho_zeta * p2ss_2);
          p1p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          d12p3s_1 += inv_two_zeta * (p1ss_1 - rho_zeta * p1ss_2);
          p2p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[7] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[7] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[7] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[7] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[7] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[7] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[2] * d12ss_1 + WmP[2] * d12ss_2;
          scalar_type p2p3s_1 = PmB[2] * p2ss_1 + WmP[2] * p2ss_2;
          scalar_type p1p3s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[8] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[8] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[8] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[8] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[8] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[8] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
      }
    }
    // START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      scalar_type p1ss_2 = PmA[2] * F_mT[2] + WmP[2] * F_mT[3];
      scalar_type p1ss_3 = PmA[2] * F_mT[3] + WmP[2] * F_mT[4];
      // START INDEX i2=0, CENTER 1
      {
        scalar_type d12ss_1 = PmA[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type d12ss_2 = PmA[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type p2ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
        scalar_type p2ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
        scalar_type norm2 = 1.0;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[0] * d12ss_1 + WmP[0] * d12ss_2;
          scalar_type p1p3s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
          scalar_type p2p3s_1 = PmB[0] * p2ss_1 + WmP[0] * p2ss_2;
          d12p3s_1 += inv_two_zeta * (p1ss_1 - rho_zeta * p1ss_2);
          p2p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[9] += (double)(preterm * fit_dens_sh[j + 0] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[9] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[9] += (double)(preterm * fit_dens_sh[j + 1] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[9] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[9] += (double)(preterm * fit_dens_sh[j + 2] *
                                   prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[9] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[1] * d12ss_1 + WmP[1] * d12ss_2;
          scalar_type p1p3s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
          scalar_type p2p3s_1 = PmB[1] * p2ss_1 + WmP[1] * p2ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[10] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[10] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[10] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[10] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[10] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[10] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[2] * d12ss_1 + WmP[2] * d12ss_2;
          scalar_type p1p3s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
          scalar_type p2p3s_1 = PmB[2] * p2ss_1 + WmP[2] * p2ss_2;
          d12p3s_1 += inv_two_zeta * (p2ss_1 - rho_zeta * p2ss_2);
          p1p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[11] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[11] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[11] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[11] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[11] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[11] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
      }
      // START INDEX i2=1, CENTER 1
      {
        scalar_type d12ss_1 = PmA[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type d12ss_2 = PmA[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type p2ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
        scalar_type p2ss_2 = PmA[1] * F_mT[2] + WmP[1] * F_mT[3];
        scalar_type norm2 = 1.0;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[0] * d12ss_1 + WmP[0] * d12ss_2;
          scalar_type p1p3s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
          scalar_type p2p3s_1 = PmB[0] * p2ss_1 + WmP[0] * p2ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[12] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[12] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[12] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[12] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[12] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[12] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[1] * d12ss_1 + WmP[1] * d12ss_2;
          scalar_type p1p3s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
          scalar_type p2p3s_1 = PmB[1] * p2ss_1 + WmP[1] * p2ss_2;
          d12p3s_1 += inv_two_zeta * (p1ss_1 - rho_zeta * p1ss_2);
          p2p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[13] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[13] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[13] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[13] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[13] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[13] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[2] * d12ss_1 + WmP[2] * d12ss_2;
          scalar_type p1p3s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
          scalar_type p2p3s_1 = PmB[2] * p2ss_1 + WmP[2] * p2ss_2;
          d12p3s_1 += inv_two_zeta * (p2ss_1 - rho_zeta * p2ss_2);
          p1p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[14] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[14] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[14] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[14] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[14] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[14] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
      }
      // START INDEX i2=2, CENTER 1
      {
        scalar_type d12ss_1 = PmA[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type d12ss_2 = PmA[2] * p1ss_2 + WmP[2] * p1ss_3;
        scalar_type p2ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
        scalar_type p2ss_2 = PmA[2] * F_mT[2] + WmP[2] * F_mT[3];
        scalar_type norm2 = 1.0;
        d12ss_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        d12ss_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        norm2 = G2G::gpu_normalization_factor;
        // START INDEX i3=0, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[0] * d12ss_1 + WmP[0] * d12ss_2;
          scalar_type p2p3s_1 = PmB[0] * p2ss_1 + WmP[0] * p2ss_2;
          scalar_type p1p3s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[15] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[15] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[15] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[15] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[15] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[15] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=1, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[1] * d12ss_1 + WmP[1] * d12ss_2;
          scalar_type p2p3s_1 = PmB[1] * p2ss_1 + WmP[1] * p2ss_2;
          scalar_type p1p3s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[16] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[16] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[16] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[16] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[16] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[16] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
        // START INDEX i3=2, CENTER 2
        {
          scalar_type d12p3s_1 = PmB[2] * d12ss_1 + WmP[2] * d12ss_2;
          scalar_type p2p3s_1 = PmB[2] * p2ss_1 + WmP[2] * p2ss_2;
          scalar_type p1p3s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
          d12p3s_1 += inv_two_zeta * (p2ss_1 - rho_zeta * p2ss_2);
          p1p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          d12p3s_1 += inv_two_zeta * (p1ss_1 - rho_zeta * p1ss_2);
          p2p3s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
          // START INDEX i4=0, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[0] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[17] += (double)(preterm * fit_dens_sh[j + 0] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[0][tid] +=
                (double)(preterm * dens[17] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=1, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[1] * d12p3s_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[17] += (double)(preterm * fit_dens_sh[j + 1] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[1][tid] +=
                (double)(preterm * dens[17] * prefactor_dens * d12p3p4_0);
#endif
          }
          // START INDEX i4=2, CENTER 3
          {
            scalar_type d12p3p4_0 = WmQ[2] * d12p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p2p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * p1p3s_1;
            d12p3p4_0 += inv_two_zeta_eta * d12ss_1;
            scalar_type preterm = norm2;
#ifdef FOCK_CALC
            my_fock[17] += (double)(preterm * fit_dens_sh[j + 2] *
                                    prefactor_dens * d12p3p4_0);
#else
            rc_sh[2][tid] +=
                (double)(preterm * dens[17] * prefactor_dens * d12p3p4_0);
#endif
          }
        }
      }
    }
  }
}
