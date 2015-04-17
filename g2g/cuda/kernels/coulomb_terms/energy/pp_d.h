{
  scalar_type F_mT[5];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type,4>(F_mT,T);
  }
  {
    //START INDEX i1=0, CENTER 1
    {
      scalar_type p1ss_2 = PmA[0] * F_mT[2] + WmP[0] * F_mT[3];
      scalar_type p1ss_3 = PmA[0] * F_mT[3] + WmP[0] * F_mT[4];
      scalar_type p1ss_0 = PmA[0] * F_mT[0] + WmP[0] * F_mT[1];
      scalar_type p1ss_1 = PmA[0] * F_mT[1] + WmP[0] * F_mT[2];
      //START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type sp2s_2 = PmB[0] * F_mT[2] + WmP[0] * F_mT[3];
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        p1p2s_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[0] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[0] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[0] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[0] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[0] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[0] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
      //START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type sp2s_2 = PmB[1] * F_mT[2] + WmP[1] * F_mT[3];
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[1] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[1] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[1] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[1] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[1] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[1] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
      //START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[2] * p1ss_2 + WmP[2] * p1ss_3;
        scalar_type sp2s_2 = PmB[2] * F_mT[2] + WmP[2] * F_mT[3];
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[0] * sp2s_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[2] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[2] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[2] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[2] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[2] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[2] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
    }
    //START INDEX i1=1, CENTER 1
    {
      scalar_type p1ss_2 = PmA[1] * F_mT[2] + WmP[1] * F_mT[3];
      scalar_type p1ss_3 = PmA[1] * F_mT[3] + WmP[1] * F_mT[4];
      scalar_type p1ss_0 = PmA[1] * F_mT[0] + WmP[1] * F_mT[1];
      scalar_type p1ss_1 = PmA[1] * F_mT[1] + WmP[1] * F_mT[2];
      //START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type sp2s_2 = PmB[0] * F_mT[2] + WmP[0] * F_mT[3];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[3] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[3] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[3] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[3] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[3] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[3] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
      //START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type sp2s_2 = PmB[1] * F_mT[2] + WmP[1] * F_mT[3];
        p1p2s_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[4] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[4] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[4] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[4] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[4] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[4] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
      //START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[2] * p1ss_2 + WmP[2] * p1ss_3;
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type sp2s_2 = PmB[2] * F_mT[2] + WmP[2] * F_mT[3];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[5] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[1] * sp2s_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[5] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[5] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[5] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[5] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[5] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
    }
    //START INDEX i1=2, CENTER 1
    {
      scalar_type p1ss_2 = PmA[2] * F_mT[2] + WmP[2] * F_mT[3];
      scalar_type p1ss_3 = PmA[2] * F_mT[3] + WmP[2] * F_mT[4];
      scalar_type p1ss_0 = PmA[2] * F_mT[0] + WmP[2] * F_mT[1];
      scalar_type p1ss_1 = PmA[2] * F_mT[1] + WmP[2] * F_mT[2];
      //START INDEX i2=0, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[0] * p1ss_2 + WmP[0] * p1ss_3;
        scalar_type p1p2s_0 = PmB[0] * p1ss_0 + WmP[0] * p1ss_1;
        scalar_type p1p2s_1 = PmB[0] * p1ss_1 + WmP[0] * p1ss_2;
        scalar_type sp2s_2 = PmB[0] * F_mT[2] + WmP[0] * F_mT[3];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[0] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[6] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[6] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[6] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[6] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[6] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[6] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
      //START INDEX i2=1, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[1] * p1ss_2 + WmP[1] * p1ss_3;
        scalar_type p1p2s_0 = PmB[1] * p1ss_0 + WmP[1] * p1ss_1;
        scalar_type p1p2s_1 = PmB[1] * p1ss_1 + WmP[1] * p1ss_2;
        scalar_type sp2s_2 = PmB[1] * F_mT[2] + WmP[1] * F_mT[3];
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[7] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[1] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[7] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[7] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[7] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[7] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[7] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
      //START INDEX i2=2, CENTER 2
      {
        scalar_type p1p2s_2 = PmB[2] * p1ss_2 + WmP[2] * p1ss_3;
        scalar_type p1p2s_0 = PmB[2] * p1ss_0 + WmP[2] * p1ss_1;
        scalar_type p1p2s_1 = PmB[2] * p1ss_1 + WmP[2] * p1ss_2;
        scalar_type sp2s_2 = PmB[2] * F_mT[2] + WmP[2] * F_mT[3];
        p1p2s_2 += inv_two_zeta * (F_mT[2] - rho_zeta * F_mT[3]);
        p1p2s_0 += inv_two_zeta * (F_mT[0] - rho_zeta * F_mT[1]);
        p1p2s_1 += inv_two_zeta * (F_mT[1] - rho_zeta * F_mT[2]);
        //START INDEX i3=0, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[0] * p1p2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[8] += (double)( preterm * fit_dens_sh[j+0] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=1, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[1] * p1p2s_2;
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[8] += (double)( preterm * fit_dens_sh[j+1] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[8] += (double)( preterm * fit_dens_sh[j+2] * prefactor_dens *  p1p2d34_0 );
          }
        }
        //START INDEX i3=2, CENTER 3
        {
          scalar_type p1p2p3_1 = WmQ[2] * p1p2s_2;
          scalar_type sp2p3_1 = WmQ[2] * sp2s_2;
          scalar_type p1sp3_1 = WmQ[2] * p1ss_2;
          p1p2p3_1 += inv_two_zeta_eta * sp2s_2;
          p1sp3_1 += inv_two_zeta_eta * F_mT[2];
          p1p2p3_1 += inv_two_zeta_eta * p1ss_2;
          sp2p3_1 += inv_two_zeta_eta * F_mT[2];
          //START INDEX i4=0, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[0] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[8] += (double)( preterm * fit_dens_sh[j+3] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=1, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[1] * p1p2p3_1;
            scalar_type norm4 = 1.0f;
            scalar_type preterm = norm4;
            my_fock[8] += (double)( preterm * fit_dens_sh[j+4] * prefactor_dens *  p1p2d34_0 );
          }
          //START INDEX i4=2, CENTER 3
          {
            scalar_type p1p2d34_0 = WmQ[2] * p1p2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * sp2p3_1;
            p1p2d34_0 += inv_two_zeta_eta * p1sp3_1;
            scalar_type norm4 = 1.0f;
            p1p2d34_0 += inv_two_eta * (p1p2s_0 - rho_eta * p1p2s_1);
            norm4 = gpu_normalization_factor;
            scalar_type preterm = norm4;
            my_fock[8] += (double)( preterm * fit_dens_sh[j+5] * prefactor_dens *  p1p2d34_0 );
          }
        }
      }
    }
  }
}
