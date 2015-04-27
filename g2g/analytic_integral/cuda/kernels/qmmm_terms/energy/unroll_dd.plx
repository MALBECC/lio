#!/usr/bin/perl -w

print <<"END";
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//--------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-DS)-------------------------------------------
        scalar_type F_mU[5];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
          lio_gamma<scalar_type,4>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          scalar_type term;

          uint curr_ind = 0;
END
#          for (uint d1_l1 = 0; d1_l1 < 3; d1_l1++)
for $d1_l1 (0..2) {
print <<"END";
          // d1_l1 == $d1_l1
          {

            scalar_type p_s0_1 = PmA[$d1_l1] * F_mU[0] - PmC[$d1_l1] * F_mU[1];
            scalar_type p_s1_1 = PmA[$d1_l1] * F_mU[1] - PmC[$d1_l1] * F_mU[2];
            scalar_type p_s2_1 = PmA[$d1_l1] * F_mU[2] - PmC[$d1_l1] * F_mU[3];
            scalar_type p_s3_1 = PmA[$d1_l1] * F_mU[3] - PmC[$d1_l1] * F_mU[4];

END
#            for (uint d1_l2 = 0; d1_l2 <= d1_l1; d1_l2++)
for $d1_l2 (0..$d1_l1) {
print <<"END";
            // d1_l2 == $d1_l2
            {

END
if ($d1_l1 != $d1_l2) {
print <<"END";
              scalar_type p_s0_2 = PmA[$d1_l2] * F_mU[0] - PmC[$d1_l2] * F_mU[1];
              scalar_type p_s1_2 = PmA[$d1_l2] * F_mU[1] - PmC[$d1_l2] * F_mU[2];
              scalar_type p_s2_2 = PmA[$d1_l2] * F_mU[2] - PmC[$d1_l2] * F_mU[3];
END
}
print <<"END";

              scalar_type d_s0 = PmA[$d1_l2] * p_s0_1 - PmC[$d1_l2] * p_s1_1;
              scalar_type d_s1 = PmA[$d1_l2] * p_s1_1 - PmC[$d1_l2] * p_s2_1;
              scalar_type d_s2 = PmA[$d1_l2] * p_s2_1 - PmC[$d1_l2] * p_s3_1;
END
if ($d1_l1 == $d1_l2) {
print <<"END";
              d_s0            += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1            += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d_s2            += inv_two_zeta * (F_mU[2] - F_mU[3]);
END
}

#              preterm1  = (d1_l1 != d1_l2) * 1.0f + (d1_l1 == d1_l2) * G2G::gpu_normalization_factor;

#              for (uint d2_l1 = 0; d2_l1 < 3; d2_l1++)
for $d2_l1 (0..2) {
print <<"END";

              // d2_l1 == $d2_l1
              {

END
if ($d1_l1 == $d2_l1 and $d1_l2 == $d2_l1) {
if ($d1_l1 != $d1_l2) {
print <<"END";
                scalar_type d_p0 = p_s0_2 - p_s1_2;
                scalar_type d_p1 = p_s1_2 - p_s2_2;
END
} else {
print <<"END";
                scalar_type d_p0 = p_s0_1 - p_s1_1;
                scalar_type d_p1 = p_s1_1 - p_s2_1;
END
}
print <<"END";
                d_p0            += p_s0_1 - p_s1_1;
                d_p1            += p_s1_1 - p_s2_1;
                d_p0            *= inv_two_zeta;
                d_p1            *= inv_two_zeta;
                d_p0            += PmB[$d2_l1] * d_s0 - PmC[$d2_l1] * d_s1;
                d_p1            += PmB[$d2_l1] * d_s1 - PmC[$d2_l1] * d_s2;
END
} elsif ($d1_l1 == $d2_l1) {
if ($d1_l1 != $d1_l2) {
print <<"END";
                scalar_type d_p0 = inv_two_zeta * (p_s0_2 - p_s1_2);
                scalar_type d_p1 = inv_two_zeta * (p_s1_2 - p_s2_2);
END
} else {
print <<"END";
                scalar_type d_p0 = inv_two_zeta * (p_s0_1 - p_s1_1);
                scalar_type d_p1 = inv_two_zeta * (p_s1_1 - p_s2_1);
END
}
print <<"END";
                d_p0            += PmB[$d2_l1] * d_s0 - PmC[$d2_l1] * d_s1;
                d_p1            += PmB[$d2_l1] * d_s1 - PmC[$d2_l1] * d_s2;
END
} elsif ($d1_l2 == $d2_l1) {
print <<"END";
                scalar_type d_p0 = inv_two_zeta * (p_s0_1 - p_s1_1);
                scalar_type d_p1 = inv_two_zeta * (p_s1_1 - p_s2_1);
                d_p0            += PmB[$d2_l1] * d_s0 - PmC[$d2_l1] * d_s1;
                d_p1            += PmB[$d2_l1] * d_s1 - PmC[$d2_l1] * d_s2;
END
} else {
print <<"END";
                scalar_type d_p0 = PmB[$d2_l1] * d_s0 - PmC[$d2_l1] * d_s1;
                scalar_type d_p1 = PmB[$d2_l1] * d_s1 - PmC[$d2_l1] * d_s2;
END
}
if ($d2_l1 >= $d1_l2) {
print <<"END";

                scalar_type p_p0_1 = PmB[$d2_l1] * p_s0_1 - PmC[$d2_l1] * p_s1_1;
                scalar_type p_p1_1 = PmB[$d2_l1] * p_s1_1 - PmC[$d2_l1] * p_s2_1;
END
if ($d1_l1 == $d2_l1) {
print <<"END";
                p_p0_1            += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_1            += inv_two_zeta * (F_mU[1] - F_mU[2]);
END
}
}
if ($d1_l1 != $d1_l2 and $d2_l1 >= $d1_l1) {
print <<"END";

                scalar_type p_p0_2 = PmB[$d2_l1] * p_s0_2 - PmC[$d2_l1] * p_s1_2;
                scalar_type p_p1_2 = PmB[$d2_l1] * p_s1_2 - PmC[$d2_l1] * p_s2_2;
END
if ($d1_l2 == $d2_l1) {
print <<"END";
                p_p0_2            += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_2            += inv_two_zeta * (F_mU[1] - F_mU[2]);
END
}
}

#                for (uint d2_l2 = 0; d2_l2 <= d2_l1; d2_l2++)
for $d2_l2 (0..$d2_l1) {
print <<"END";
                // d2_l2 == $d2_l2
                {

END
if ($d1_l1 == $d2_l2 and $d1_l2 == $d2_l2 and $d2_l1 == $d2_l2) {
if ($d1_l1 != $d1_l2) {
print <<"END";
                  term  = p_p0_2 - p_p1_2;
END
} else {
print <<"END";
                  term  = p_p0_1 - p_p1_1;
END
}
print <<"END"
                  term += p_p0_1 - p_p1_1;
                  term += d_s0 - d_s1;
                  term *= inv_two_zeta;
                  term += PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
} elsif ($d1_l1 == $d2_l2 and $d1_l2 == $d2_l2) {
if ($d1_l1 != $d1_l2) {
print <<"END";
                  term  = p_p0_2 - p_p1_2;
END
} else {
print <<"END";
                  term  = p_p0_1 - p_p1_1;
END
}
print <<"END"
                  term += p_p0_1 - p_p1_1;
                  term *= inv_two_zeta;
                  term += PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
} elsif ($d1_l1 == $d2_l2 and $d2_l1 == $d2_l2) {
if ($d1_l1 != $d1_l2) {
print <<"END";
                  term  = p_p0_2 - p_p1_2;
END
} else {
print <<"END";
                  term  = p_p0_1 - p_p1_1;
END
}
print <<"END";
                  term += d_s0 - d_s1;
                  term *= inv_two_zeta;
                  term += PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
} elsif ($d1_l2 == $d2_l2 and $d2_l1 == $d2_l2) {
print <<"END";
                  term  = p_p0_1 - p_p1_1;
                  term += d_s0 - d_s1;
                  term *= inv_two_zeta;
                  term += PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
} elsif ($d1_l1 == $d2_l2) {
if ($d1_l1 != $d1_l2) {
print <<"END";
                  term  = inv_two_zeta * (p_p0_2 - p_p1_2);
END
} else {
print <<"END";
                  term  = inv_two_zeta * (p_p0_1 - p_p1_1);
END
}
print <<"END";
                  term += PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
} elsif ($d1_l2 == $d2_l2) {
print <<"END";
                  term  = inv_two_zeta * (p_p0_1 - p_p1_1);
                  term += PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
} elsif ($d2_l1 == $d2_l2) {
print <<"END";
                  term  = inv_two_zeta * (d_s0 - d_s1);
                  term += PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
} else {
print <<"END";
                  term  = PmB[$d2_l2] * d_p0 - PmC[$d2_l2] * d_p1;
END
}
if ($d1_l1 == $d1_l2) {
print <<"END";
                  term *= G2G::gpu_normalization_factor;
END
}
if ($d2_l1 == $d2_l2) {
print <<"END";
                  term *= G2G::gpu_normalization_factor;
END
}
if ($d2_l1 > $d1_l1 or ($d2_l1 == $d1_l1 and $d2_l2 > $d1_l2)) {
print <<"END";

                  my_fock[curr_ind] += !same_func * clatom_charge_sh[j] * term;
                  curr_ind += !same_func;
END
} else {
print <<"END";

                  my_fock[curr_ind] += clatom_charge_sh[j] * term;
                  curr_ind++;
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
