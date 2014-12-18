// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (S-S)-------------------------------------------
        scalar_type F_mU[4];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
          lio_gamma<scalar_type,3>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          scalar_type term;

          uint curr_ind = 0;
          // d_l1 == 0
          {
            scalar_type p_s0_1 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
            scalar_type p_s1_1 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
            scalar_type p_s2_1 = PmA[0] * F_mU[2] - PmC[0] * F_mU[3];
            // d_l2 == 0
            {
              scalar_type d_s0 = PmA[0] * p_s0_1 - PmC[0] * p_s1_1;
              scalar_type d_s1 = PmA[0] * p_s1_1 - PmC[0] * p_s2_1;
              d_s0            += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1            += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // p_l == 0
              {

                term  = p_s0_1 - p_s1_1;
                term += p_s0_1 - p_s1_1;
                term *= inv_two_zeta;
                term += PmB[0] * d_s0 - PmC[0] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 1
              {

                term  = PmB[1] * d_s0 - PmC[1] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 2
              {

                term  = PmB[2] * d_s0 - PmC[2] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
            }
          }
          // d_l1 == 1
          {
            scalar_type p_s0_1 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
            scalar_type p_s1_1 = PmA[1] * F_mU[1] - PmC[1] * F_mU[2];
            scalar_type p_s2_1 = PmA[1] * F_mU[2] - PmC[1] * F_mU[3];
            // d_l2 == 0
            {
              scalar_type d_s0 = PmA[0] * p_s0_1 - PmC[0] * p_s1_1;
              scalar_type d_s1 = PmA[0] * p_s1_1 - PmC[0] * p_s2_1;

              scalar_type p_s0_2 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
              scalar_type p_s1_2 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
              // p_l == 0
              {

                term  = inv_two_zeta * (p_s0_1 - p_s1_1);
                term += PmB[0] * d_s0 - PmC[0] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 1
              {

                term  = inv_two_zeta * (p_s0_2 - p_s1_2);
                term += PmB[1] * d_s0 - PmC[1] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 2
              {

                term  = PmB[2] * d_s0 - PmC[2] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
            }
            // d_l2 == 1
            {
              scalar_type d_s0 = PmA[1] * p_s0_1 - PmC[1] * p_s1_1;
              scalar_type d_s1 = PmA[1] * p_s1_1 - PmC[1] * p_s2_1;
              d_s0            += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1            += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // p_l == 0
              {

                term  = PmB[0] * d_s0 - PmC[0] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 1
              {

                term  = p_s0_1 - p_s1_1;
                term += p_s0_1 - p_s1_1;
                term *= inv_two_zeta;
                term += PmB[1] * d_s0 - PmC[1] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 2
              {

                term  = PmB[2] * d_s0 - PmC[2] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
            }
          }
          // d_l1 == 2
          {
            scalar_type p_s0_1 = PmA[2] * F_mU[0] - PmC[2] * F_mU[1];
            scalar_type p_s1_1 = PmA[2] * F_mU[1] - PmC[2] * F_mU[2];
            scalar_type p_s2_1 = PmA[2] * F_mU[2] - PmC[2] * F_mU[3];
            // d_l2 == 0
            {
              scalar_type d_s0 = PmA[0] * p_s0_1 - PmC[0] * p_s1_1;
              scalar_type d_s1 = PmA[0] * p_s1_1 - PmC[0] * p_s2_1;

              scalar_type p_s0_2 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
              scalar_type p_s1_2 = PmA[0] * F_mU[1] - PmC[0] * F_mU[2];
              // p_l == 0
              {

                term  = inv_two_zeta * (p_s0_1 - p_s1_1);
                term += PmB[0] * d_s0 - PmC[0] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 1
              {

                term  = PmB[1] * d_s0 - PmC[1] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 2
              {

                term  = inv_two_zeta * (p_s0_2 - p_s1_2);
                term += PmB[2] * d_s0 - PmC[2] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
            }
            // d_l2 == 1
            {
              scalar_type d_s0 = PmA[1] * p_s0_1 - PmC[1] * p_s1_1;
              scalar_type d_s1 = PmA[1] * p_s1_1 - PmC[1] * p_s2_1;

              scalar_type p_s0_2 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
              scalar_type p_s1_2 = PmA[1] * F_mU[1] - PmC[1] * F_mU[2];
              // p_l == 0
              {

                term  = PmB[0] * d_s0 - PmC[0] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 1
              {

                term  = inv_two_zeta * (p_s0_1 - p_s1_1);
                term += PmB[1] * d_s0 - PmC[1] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 2
              {

                term  = inv_two_zeta * (p_s0_2 - p_s1_2);
                term += PmB[2] * d_s0 - PmC[2] * d_s1;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
            }
            // d_l2 == 2
            {
              scalar_type d_s0 = PmA[2] * p_s0_1 - PmC[2] * p_s1_1;
              scalar_type d_s1 = PmA[2] * p_s1_1 - PmC[2] * p_s2_1;
              d_s0            += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1            += inv_two_zeta * (F_mU[1] - F_mU[2]);
              // p_l == 0
              {

                term  = PmB[0] * d_s0 - PmC[0] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 1
              {

                term  = PmB[1] * d_s0 - PmC[1] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
              // p_l == 2
              {

                term  = p_s0_1 - p_s1_1;
                term += p_s0_1 - p_s1_1;
                term *= inv_two_zeta;
                term += PmB[2] * d_s0 - PmC[2] * d_s1;
                term *= gpu_normalization_factor;

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
            }
          }
        }
        // END individual force terms 
//------------------------------------------END TERM-TYPE DEPENDENT PART (S-S)----------------------------------------------
