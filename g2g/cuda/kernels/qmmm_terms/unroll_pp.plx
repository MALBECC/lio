#!/usr/bin/perl -w

print <<"END";
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (P-P)-------------------------------------------
        scalar_type F_mU[4];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * zeta;
          //for (int m = 0; m <= 3; m++) 
          //{
          //  F_mU[m] = lio_gamma<scalar_type>(m,U);
          lio_gamma<scalar_type,3>(F_mU,U);
          //}
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common, p_p0, p_p1;
          scalar_type mm_charge = clatom_charge_sh[j];
          uint dens_ind = 0;

END
#          for (int p1_l = 0; p1_l < 3; p1_l++)
for $p1_l (0..2) {
print <<"END";
          // p1_l == $p1_l
          {
END
#            for (int p2_l = 0; p2_l < 3; p2_l++)
for $p2_l (0..2) {
print <<"END";
            // p2_l == $p2_l
            {
              scalar_type pre_term;
              {
END
if ($p2_l > $p1_l) {
#                bool skip = same_func && ($p2_l > $p1_l);
print <<"END";
                pre_term  = !same_func * mm_charge * dens[dens_ind];
                dens_ind += !same_func;
END
} else {
print <<"END";
                pre_term  = mm_charge * dens[dens_ind];
                dens_ind++;
END
}
print <<"END";
              }

              p_p0  = PmB[$p2_l] * (PmA[$p1_l] * F_mU[0] - PmC[$p1_l] * F_mU[1]); // p_s0
              p_p0 -= PmC[$p2_l] * (PmA[$p1_l] * F_mU[1] - PmC[$p1_l] * F_mU[2]); // p_s1
              p_p1  = PmB[$p2_l] * (PmA[$p1_l] * F_mU[1] - PmC[$p1_l] * F_mU[2]); // p_s1
              p_p1 -= PmC[$p2_l] * (PmA[$p1_l] * F_mU[2] - PmC[$p1_l] * F_mU[3]); // p_s2
END
if ($p2_l == $p1_l) {
print <<"END";
              p_p0 += inv_two_zeta * (F_mU[0] - F_mU[1]);
              p_p1 += inv_two_zeta * (F_mU[1] - F_mU[2]);
END
}
#              for (int grad_l = 0; grad_l < 3; grad_l++)
for $grad_l (0..2) {
print <<"END";
              // grad_l == $grad_l
              {
END
if ($p1_l == $grad_l) {
print <<"END";
                C_force_term  = PmB[$p2_l] * F_mU[1] - PmC[$p2_l] * F_mU[2]; // p_s1 (B)
                AB_common     = PmB[$p2_l] * F_mU[0] - PmC[$p2_l] * F_mU[1]; // p_s0 (B)
END
} elsif ($p2_l == $grad_l) {
print <<"END";
                C_force_term  = 0.0f;
                AB_common     = 0.0f;
END
}
if ($p2_l == $grad_l) {
print <<"END";
                C_force_term += PmA[$p1_l] * F_mU[1] - PmC[$p1_l] * F_mU[2]; // p_s1
                AB_common    += PmA[$p1_l] * F_mU[0] - PmC[$p1_l] * F_mU[1]; // p_s0
END
}
if ($p1_l != $grad_l and $p2_l != $grad_l) {
print <<"END";
                C_force_term  = PmC[$grad_l] * p_p1;

                A_force_term  = -C_force_term;
END
} else {
print <<"END";
                C_force_term  = PmC[$grad_l] * p_p1 + inv_two_zeta * C_force_term;

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
END
}
print <<"END";
                B_force_term  = PmB[$grad_l] * p_p0 + A_force_term;
                A_force_term += PmA[$grad_l] * p_p0;
                A_force_term *= 2.0f * ai;
                B_force_term *= 2.0f * aj;
END
if ($p1_l == $grad_l) {
print <<"END";
                A_force_term -= PmB[$p2_l] * F_mU[0] - PmC[$p2_l] * F_mU[1]; // p_s0 (B)
END
}
if ($p2_l == $grad_l) {
print <<"END";
                B_force_term -= PmA[$p1_l] * F_mU[0] - PmC[$p1_l] * F_mU[1]; // p_s0
END
}
print <<"END";
                A_force[$grad_l]     += pre_term * A_force_term;
                B_force[$grad_l]     += pre_term * B_force_term;
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
//-------------------------------------------END TERM-TYPE DEPENDENT PART (P-P)-------------------------------------------
END
