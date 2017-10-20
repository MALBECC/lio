{
  scalar_type F_mT[4];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type,3>(F_mT,T);
  }
  {
    //START INDEX i1=0, CENTER 1
    {
      scalar_type p1ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
      scalar_type p1ss_0 = PmA[0] * F_mT[0] + WmP[0] * F_mT[1];
      scalar_type p1ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
      //START INDEX i2=0, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[0] * p1ss_2;
        scalar_type ssp2_1 = WmQ[0] * F_mT[2];
        p1sp2_1 += inv_two_zeta_eta * F_mT[2];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          p1sd23_0 += inv_two_zeta_eta * ssp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[0] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[0][tid] += (double)( preterm * dens[0] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
      //START INDEX i2=1, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[1] * p1ss_2;
        scalar_type ssp2_1 = WmQ[1] * F_mT[2];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          p1sd23_0 += inv_two_zeta_eta * ssp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[0] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[1][tid] += (double)( preterm * dens[0] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[1] * p1sp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[0] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[2][tid] += (double)( preterm * dens[0] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
      //START INDEX i2=2, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[2] * p1ss_2;
        scalar_type ssp2_1 = WmQ[2] * F_mT[2];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          p1sd23_0 += inv_two_zeta_eta * ssp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[0] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[3][tid] += (double)( preterm * dens[0] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[1] * p1sp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[0] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[4][tid] += (double)( preterm * dens[0] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[2] * p1sp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[0] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[5][tid] += (double)( preterm * dens[0] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
    }
    //START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_2 = PmA[1] * F_mT[2] + WmP[1] * F_mT[3];
      scalar_type p1ss_0 = PmA[1] * F_mT[0] + WmP[1] * F_mT[1];
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      //START INDEX i2=0, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[0] * p1ss_2;
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[1] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[0][tid] += (double)( preterm * dens[1] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
      //START INDEX i2=1, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[1] * p1ss_2;
        scalar_type ssp2_1 = WmQ[1] * F_mT[2];
        p1sp2_1 += inv_two_zeta_eta * F_mT[2];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[1] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[1][tid] += (double)( preterm * dens[1] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[1] * p1sp2_1;
          p1sd23_0 += inv_two_zeta_eta * ssp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[1] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[2][tid] += (double)( preterm * dens[1] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
      //START INDEX i2=2, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[2] * p1ss_2;
        scalar_type ssp2_1 = WmQ[2] * F_mT[2];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[1] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[3][tid] += (double)( preterm * dens[1] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[1] * p1sp2_1;
          p1sd23_0 += inv_two_zeta_eta * ssp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[1] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[4][tid] += (double)( preterm * dens[1] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[2] * p1sp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[1] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[5][tid] += (double)( preterm * dens[1] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
    }
    //START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_2 = PmA[2] * F_mT[2] + WmP[2] * F_mT[3];
      scalar_type p1ss_0 = PmA[2] * F_mT[0] + WmP[2] * F_mT[1];
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      //START INDEX i2=0, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[0] * p1ss_2;
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[2] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[0][tid] += (double)( preterm * dens[2] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
      //START INDEX i2=1, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[1] * p1ss_2;
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[2] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[1][tid] += (double)( preterm * dens[2] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[1] * p1sp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[2] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[2][tid] += (double)( preterm * dens[2] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
      //START INDEX i2=2, CENTER 3
      {
        scalar_type p1sp2_1 = WmQ[2] * p1ss_2;
        scalar_type ssp2_1 = WmQ[2] * F_mT[2];
        p1sp2_1 += inv_two_zeta_eta * F_mT[2];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[0] * p1sp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[2] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[3][tid] += (double)( preterm * dens[2] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[1] * p1sp2_1;
          scalar_type norm3 = 1.0;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[2] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[4][tid] += (double)( preterm * dens[2] * prefactor_dens *  p1sd23_0 );
#endif
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1sd23_0 = WmQ[2] * p1sp2_1;
          p1sd23_0 += inv_two_zeta_eta * ssp2_1;
          scalar_type norm3 = 1.0;
          p1sd23_0 += inv_two_eta * (p1ss_0 - rho_eta * p1ss_1);
          norm3 = G2G::gpu_normalization_factor;
          scalar_type preterm = norm3;
#ifdef FOCK_CALC
          my_fock[2] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1sd23_0 );
#else
          rc_sh[5][tid] += (double)( preterm * dens[2] * prefactor_dens *  p1sd23_0 );
#endif
        }
      }
    }
  }
}
