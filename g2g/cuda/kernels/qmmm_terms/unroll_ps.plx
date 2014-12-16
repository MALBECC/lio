#!/usr/bin/perl -w

print <<"END";
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (P-S)-------------------------------------------
        scalar_type F_mU[3];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          lio_gamma<scalar_type,2>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
END
for $p_l (0..2) {
  print <<"END";
          // p_l == $p_l
          {
            scalar_type pre_term = clatom_charge_sh[j] * dens[$p_l];
            scalar_type p_s0 = PmA[$p_l] * F_mU[0] - PmC[$p_l] * F_mU[1];
            scalar_type p_s1 = PmA[$p_l] * F_mU[1] - PmC[$p_l] * F_mU[2];
END
  for $grad_l (0..2) {
    print <<"END";
            // grad_l == $grad_l
            {
              C_force_term  = PmC[$grad_l] * p_s1;
END
    if ($p_l == $grad_l) {
      print <<"END";
              C_force_term += inv_two_zeta * F_mU[1];

END
    }
    if ($p_l == $grad_l) {
    print <<"END";
              A_force_term  = inv_two_zeta * F_mU[0] - C_force_term;
END
    } else {
    print <<"END";

              A_force_term  = -C_force_term;
END
    }
    print <<"END";
              B_force_term  = PmB[$grad_l] * p_s0 + A_force_term;
              A_force_term += PmA[$grad_l] * p_s0;

END
    if ($p_l == $grad_l) {
      print <<"END";
              A_force[$grad_l]     += pre_term * (2.0f * ai * A_force_term - F_mU[0]);
END
    } else {
      print <<"END";
              A_force[$grad_l]     += pre_term * (2.0f * ai * A_force_term);
END
    }
    print <<"END";
              B_force[$grad_l]     += pre_term * (2.0f * aj * B_force_term);
              C_force[$grad_l][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
            }
END
  }
  print <<"END";
          }
END
}
print <<"END";
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (P-S)-------------------------------------------
END
