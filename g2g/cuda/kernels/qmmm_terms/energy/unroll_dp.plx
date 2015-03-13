#!/usr/bin/perl -w

print <<"END";
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
END
#          for (uint d_l1 = 0; d_l1 < 3; d_l1++) {
#$curr_ind = 0;
for $d_l1 (0..2) {
print <<"END";
          // d_l1 == $d_l1
          {
            scalar_type p_s0_1 = PmA[$d_l1] * F_mU[0] - PmC[$d_l1] * F_mU[1];
            scalar_type p_s1_1 = PmA[$d_l1] * F_mU[1] - PmC[$d_l1] * F_mU[2];
            scalar_type p_s2_1 = PmA[$d_l1] * F_mU[2] - PmC[$d_l1] * F_mU[3];
END
#            for (uint d_l2 = 0; d_l2 <= d_l1; d_l2++) {
for $d_l2 (0..$d_l1) {
print <<"END";
            // d_l2 == $d_l2
            {
              scalar_type d_s0 = PmA[$d_l2] * p_s0_1 - PmC[$d_l2] * p_s1_1;
              scalar_type d_s1 = PmA[$d_l2] * p_s1_1 - PmC[$d_l2] * p_s2_1;
END
if ($d_l1 == $d_l2) {
print <<"END";
              d_s0            += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1            += inv_two_zeta * (F_mU[1] - F_mU[2]);
END
}
if ($d_l1 != $d_l2) {
print <<"END";

              scalar_type p_s0_2 = PmA[$d_l2] * F_mU[0] - PmC[$d_l2] * F_mU[1];
              scalar_type p_s1_2 = PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2];
END
}
#              preterm  = (d_l1 != d_l2) * 1.0f + (d_l1 == d_l2) * gpu_normalization_factor;
#             for (uint p_l = 0; p_l < 3; p_l++) {
for $p_l (0..2) {
print <<"END";
              // p_l == $p_l
              {
END
if ($d_l1 == $p_l and $d_l2 == $p_l) {
if ($d_l1 != $d_l2) {
print <<"END";

                term  = p_s0_2 - p_s1_2;
END
} else {
print <<"END";

                term  = p_s0_1 - p_s1_1;
END
}
print <<"END";
                term += p_s0_1 - p_s1_1;
                term *= inv_two_zeta;
                term += PmB[$p_l] * d_s0 - PmC[$p_l] * d_s1;
END
} elsif ($d_l1 == $p_l) {
if ($d_l1 != $d_l2) {
print <<"END";

                term  = inv_two_zeta * (p_s0_2 - p_s1_2);
END
} else {
print <<"END";

                term  = inv_two_zeta * (p_s0_1 - p_s1_1);
END
}
print <<"END";
                term += PmB[$p_l] * d_s0 - PmC[$p_l] * d_s1;
END
} elsif ($d_l2 == $p_l) {
print <<"END";

                term  = inv_two_zeta * (p_s0_1 - p_s1_1);
                term += PmB[$p_l] * d_s0 - PmC[$p_l] * d_s1;
END
} else {
print <<"END";

                term  = PmB[$p_l] * d_s0 - PmC[$p_l] * d_s1;
END
}
if ($d_l1 == $d_l2) {
print <<"END";
                term *= gpu_normalization_factor;
END
}
print <<"END";

                my_fock[curr_ind] += clatom_charge_sh[j] * term;
                curr_ind++;
              }
END
#$curr_ind += 1;
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
//------------------------------------------END TERM-TYPE DEPENDENT PART (S-S)----------------------------------------------
END
