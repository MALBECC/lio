#!/usr/bin/perl -w

print <<"END";
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-S)-------------------------------------------
        scalar_type F_mU[4];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          lio_gamma<scalar_type,3>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common, d_s0, d_s1;
          uint dens_ind = 0;

END
#          for (int d_l1 = 0; d_l1 < 3; d_l1++)
for $d_l1 (0..2) {
print <<"END";
          // d_l1 == $d_l1
          {
END
#            for (int d_l2 = 0; d_l2 <= d_l1; d_l2++)
for $d_l2 (0..$d_l1) {
print <<"END";
            // d_l2 == $d_l2
            {
END
if ($d_l1 == $d_l2) {
print <<"END";
              scalar_type pre_term = gpu_normalization_factor * clatom_charge_sh[j] * dens[dens_ind];
END
} else {
print <<"END";
              scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
END
}
print <<"END";
              dens_ind++;

              d_s0  = PmA[$d_l1] * (PmA[$d_l2] * F_mU[0] - PmC[$d_l2] * F_mU[1]); // p_s0 (d_l2)
              d_s0 -= PmC[$d_l1] * (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]); // p_s1 (d_l2)
              d_s1  = PmA[$d_l1] * (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]); // p_s1 (d_l2)
              d_s1 -= PmC[$d_l1] * (PmA[$d_l2] * F_mU[2] - PmC[$d_l2] * F_mU[3]); // p_s2 (d_l2)
END
if ($d_l1 == $d_l2) {
print <<"END";
              d_s0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
END
}
#              for (int grad_l = 0; grad_l < 3; grad_l++)
for $grad_l (0..2) {
print <<"END";
              // grad_l == $grad_l
              {
END
if ($d_l1 == $grad_l) {
print <<"END";
                C_force_term  = PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]; // p_s1 (d_l2)
                AB_common     = PmA[$d_l2] * F_mU[0] - PmC[$d_l2] * F_mU[1]; // p_s0 (d_l2)
END
} elsif ($d_l2 == $grad_l) {
print <<"END";
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
END
}
if ($d_l2 == $grad_l) {
print <<"END";
                C_force_term += PmA[$d_l1] * F_mU[1] - PmC[$d_l1] * F_mU[2]; // p_s1 (d_l1)
                AB_common    += PmA[$d_l1] * F_mU[0] - PmC[$d_l1] * F_mU[1]; // p_s0 (d_l1)
END
}
if ($d_l1 != $grad_l and $d_l2 != $grad_l) {
print <<"END";
                C_force_term  = PmC[$grad_l] * d_s1;

                A_force_term  = -C_force_term;
END
} else {
print <<"END";
                C_force_term  = PmC[$grad_l] * d_s1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
END
}
print <<"END";
                B_force_term  = PmB[$grad_l] * d_s0 + A_force_term;
                A_force_term += PmA[$grad_l] * d_s0;

                A_force_term *= 2.0f * ai;
END
if ($d_l1 == $grad_l) {
print <<"END";
                A_force_term -= PmA[$d_l2] * F_mU[0] - PmC[$d_l2] * F_mU[1]; // p_s0 (d_l2)
END
}
if ($d_l2 == $grad_l) {
print <<"END";
                A_force_term -= PmA[$d_l1] * F_mU[0] - PmC[$d_l1] * F_mU[1]; // p_s0 (d_l1)
END
}
print <<"END";
                A_force[$grad_l]     += pre_term * A_force_term;
                B_force[$grad_l]     += pre_term * 2.0f * aj * B_force_term;
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
END
}
print <<"END";
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-S)-------------------------------------------
END
