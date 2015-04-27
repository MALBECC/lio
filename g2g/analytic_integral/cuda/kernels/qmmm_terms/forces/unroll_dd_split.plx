#!/usr/bin/perl -w

print <<"END";
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-D)-------------------------------------------
        scalar_type F_mU[6];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
          lio_gamma<scalar_type,5>(F_mU,U);
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common;
          uint dens_ind = 0;

          bool del_d1 = d1_l1 == d1_l2;
          scalar_type d1_s0  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]); // p_s0 (d1_l2)
          d1_s0             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s1 (d1_l2)
          d1_s0             += del_d1 * inv_two_zeta * (F_mU[0] - F_mU[1]);
          scalar_type d1_s1  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s1 (d1_l2)
          d1_s1             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s2 (d1_l2)
          d1_s1             += del_d1 * inv_two_zeta * (F_mU[1] - F_mU[2]);
          scalar_type d1_s2  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s2 (d1_l2)
          d1_s2             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]); // p_s3 (d1_l2)
          d1_s2             += del_d1 * inv_two_zeta * (F_mU[2] - F_mU[3]);
          scalar_type d1_s3  = PmA[d1_l1] * (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]); // p_s3 (d1_l2)
          d1_s3             -= PmC[d1_l1] * (PmA[d1_l2] * F_mU[4] - PmC[d1_l2] * F_mU[5]); // p_s4 (d1_l2)
          d1_s3             += del_d1 * inv_two_zeta * (F_mU[3] - F_mU[4]);
END
#          for (int d2_l1 = 0; d2_l1 < 3; d2_l1++)
for $d2_l1 (0..2) {
print <<"END";
          // d2_l1 = $d2_l1
          {
            bool del_d1l1_d2l1 = d1_l1 == $d2_l1, del_d1l2_d2l1 = d1_l2 == $d2_l1;

            scalar_type p_p0_d1l1_d2l1  = PmB[$d2_l1] * (PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]); // p_s0 (d1_l1)
            p_p0_d1l1_d2l1             -= PmC[$d2_l1] * (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]); // p_s1 (d1_l1)
            p_p0_d1l1_d2l1             += del_d1l1_d2l1 * inv_two_zeta * (F_mU[0] - F_mU[1]);
            scalar_type p_p1_d1l1_d2l1  = PmB[$d2_l1] * (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]); // p_s1 (d1_l1)
            p_p1_d1l1_d2l1             -= PmC[$d2_l1] * (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]); // p_s2 (d1_l1)
            p_p1_d1l1_d2l1             += del_d1l1_d2l1 * inv_two_zeta * (F_mU[1] - F_mU[2]);
            scalar_type p_p2_d1l1_d2l1  = PmB[$d2_l1] * (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]); // p_s2 (d1_l1)
            p_p2_d1l1_d2l1             -= PmC[$d2_l1] * (PmA[d1_l1] * F_mU[3] - PmC[d1_l1] * F_mU[4]); // p_s3 (d1_l1)
            p_p2_d1l1_d2l1             += del_d1l1_d2l1 * inv_two_zeta * (F_mU[2] - F_mU[3]);
            scalar_type p_p0_d1l2_d2l1  = PmB[$d2_l1] * (PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]); // p_s0 (d1_l2)
            p_p0_d1l2_d2l1             -= PmC[$d2_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s1 (d1_l2)
            p_p0_d1l2_d2l1             += del_d1l2_d2l1 * inv_two_zeta * (F_mU[0] - F_mU[1]);
            scalar_type p_p1_d1l2_d2l1  = PmB[$d2_l1] * (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]); // p_s1 (d1_l2)
            p_p1_d1l2_d2l1             -= PmC[$d2_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s2 (d1_l2)
            p_p1_d1l2_d2l1             += del_d1l2_d2l1 * inv_two_zeta * (F_mU[1] - F_mU[2]);
            scalar_type p_p2_d1l2_d2l1  = PmB[$d2_l1] * (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]); // p_s2 (d1_l2)
            p_p2_d1l2_d2l1             -= PmC[$d2_l1] * (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]); // p_s3 (d1_l2)
            p_p2_d1l2_d2l1             += del_d1l2_d2l1 * inv_two_zeta * (F_mU[2] - F_mU[3]);

            scalar_type d1_p0_d2l1       = del_d1l1_d2l1 * ((PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]) - (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
            d1_p0_d2l1                  += del_d1l2_d2l1 * ((PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]) - (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
            d1_p0_d2l1                  *= inv_two_zeta;
            d1_p0_d2l1                  += PmB[$d2_l1] * d1_s0 - PmC[$d2_l1] * d1_s1;
            scalar_type d1_p1_d2l1       = del_d1l1_d2l1 * ((PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]) - (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]));  // p_s1 (d1_l2) - p_s2 (d1_l2)
            d1_p1_d2l1                  += del_d1l2_d2l1 * ((PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]) - (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]));  // p_s1 (d1_l1) - p_s2 (d1_l1)
            d1_p1_d2l1                  *= inv_two_zeta;
            d1_p1_d2l1                  += PmB[$d2_l1] * d1_s1 - PmC[$d2_l1] * d1_s2;
            scalar_type d1_p2_d2l1       = del_d1l1_d2l1 * ((PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]) - (PmA[d1_l2] * F_mU[3] - PmC[d1_l2] * F_mU[4]));  // p_s2 (d1_l2) - p_s3 (d1_l2)
            d1_p2_d2l1                  += del_d1l2_d2l1 * ((PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]) - (PmA[d1_l1] * F_mU[3] - PmC[d1_l1] * F_mU[4]));  // p_s2 (d1_l1) - p_s3 (d1_l1)
            d1_p2_d2l1                  *= inv_two_zeta;
            d1_p2_d2l1                  += PmB[$d2_l1] * d1_s2 - PmC[$d2_l1] * d1_s3;
END
#            for (int d2_l2 = 0; d2_l2 <= $d2_l1; d2_l2++)
for $d2_l2 (0..$d2_l1) {
print <<"END";
            // d2_l2 = $d2_l2
            {
              bool del_d1l1_d2l2 = d1_l1 == $d2_l2, del_d1l2_d2l2 = d1_l2 == $d2_l2;
              scalar_type pre_term;
              {
END
if ($d2_l1 == 0) {
print <<"END";
                pre_term  = clatom_charge_sh[j] * dens[dens_ind];
                dens_ind++;
END
} elsif ($d2_l2 == 0) {
print <<"END";
                bool skip = same_func && ($d2_l1 > d1_l1);
                pre_term  = !skip * clatom_charge_sh[j] * dens[dens_ind];
                dens_ind += !skip;
END
} else {
print <<"END";
                bool skip = same_func && ($d2_l1 > d1_l1 || ($d2_l1 == d1_l1 && $d2_l2 > d1_l2));
                pre_term  = !skip * clatom_charge_sh[j] * dens[dens_ind];
                dens_ind += !skip;
END
}
if ($d2_l1 == $d2_l2) {
print <<"END";
                pre_term *= !del_d1*G2G::gpu_normalization_factor + del_d1*G2G::gpu_normalization_factor*G2G::gpu_normalization_factor;
END
} else {
print <<"END";
                pre_term *= !del_d1*1.0f + del_d1*G2G::gpu_normalization_factor;
END
}
print <<"END";
              }

              scalar_type d1_p0_d2l2      = del_d1l1_d2l2 * ((PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]) - (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]));  // p_s0 (d1_l2) - p_s1 (d1_l2)
              d1_p0_d2l2                 += del_d1l2_d2l2 * ((PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]) - (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]));  // p_s0 (d1_l1) - p_s1 (d1_l1)
              d1_p0_d2l2                 *= inv_two_zeta;
              d1_p0_d2l2                 += PmB[$d2_l2] * d1_s0 - PmC[$d2_l2] * d1_s1;
              scalar_type d1_p1_d2l2      = del_d1l1_d2l2 * ((PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]) - (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]));  // p_s1 (d1_l2) - p_s2 (d1_l2)
              d1_p1_d2l2                 += del_d1l2_d2l2 * ((PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]) - (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]));  // p_s1 (d1_l1) - p_s2 (d1_l1)
              d1_p1_d2l2                 *= inv_two_zeta;
              d1_p1_d2l2                 += PmB[$d2_l2] * d1_s1 - PmC[$d2_l2] * d1_s2;
              scalar_type p_d2_0_d1l1     = del_d1l1_d2l2 * ((PmB[$d2_l1] * F_mU[0] - PmC[$d2_l1] * F_mU[1]) - (PmB[$d2_l1] * F_mU[1] - PmC[$d2_l1] * F_mU[2]));  // s_p0 (d2_l1) - s_p1 (d2_l1)
END
if ($d2_l1 == $d2_l2) {
print <<"END";
              p_d2_0_d1l1                += (PmA[d1_l1] * F_mU[0] - PmC[d1_l1] * F_mU[1]) - (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]);  // p_s0 (d1_l1) - p_s1 (d1_l1)
END
}
print <<"END";
              p_d2_0_d1l1                *= inv_two_zeta;
              p_d2_0_d1l1                += PmB[$d2_l2] * p_p0_d1l1_d2l1 - PmC[$d2_l2] * p_p1_d1l1_d2l1;
              scalar_type p_d2_1_d1l1     = del_d1l1_d2l2 * ((PmB[$d2_l1] * F_mU[1] - PmC[$d2_l1] * F_mU[2]) - (PmB[$d2_l1] * F_mU[2] - PmC[$d2_l1] * F_mU[3]));  // s_p1 (d2_l1) - s_p2 (d2_l1)
END
if ($d2_l1 == $d2_l2) {
print <<"END";
              p_d2_1_d1l1                += (PmA[d1_l1] * F_mU[1] - PmC[d1_l1] * F_mU[2]) - (PmA[d1_l1] * F_mU[2] - PmC[d1_l1] * F_mU[3]);  // p_s1 (d1_l1) - p_s2 (d1_l1)
END
}
print <<"END";
              p_d2_1_d1l1                *= inv_two_zeta;
              p_d2_1_d1l1                += PmB[$d2_l2] * p_p1_d1l1_d2l1 - PmC[$d2_l2] * p_p2_d1l1_d2l1;
              scalar_type p_d2_0_d1l2     = del_d1l2_d2l2 * ((PmB[$d2_l1] * F_mU[0] - PmC[$d2_l1] * F_mU[1]) - (PmB[$d2_l1] * F_mU[1] - PmC[$d2_l1] * F_mU[2]));  // s_p0 (d2_l1) - s_p1 (d2_l1)
END
if ($d2_l1 == $d2_l2) {
print <<"END";
              p_d2_0_d1l2                += (PmA[d1_l2] * F_mU[0] - PmC[d1_l2] * F_mU[1]) - (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]);  // p_s0 (d1_l2) - p_s1 (d1_l2)
END
}
print <<"END";
              p_d2_0_d1l2                *= inv_two_zeta;
              p_d2_0_d1l2                += PmB[$d2_l2] * p_p0_d1l2_d2l1 - PmC[$d2_l2] * p_p1_d1l2_d2l1;
              scalar_type p_d2_1_d1l2     = del_d1l2_d2l2 * ((PmB[$d2_l1] * F_mU[1] - PmC[$d2_l1] * F_mU[2]) - (PmB[$d2_l1] * F_mU[2] - PmC[$d2_l1] * F_mU[3]));  // s_p1 (d2_l1) - s_p2 (d2_l1)
END
if ($d2_l1 == $d2_l2) {
print <<"END";
              p_d2_1_d1l2                += (PmA[d1_l2] * F_mU[1] - PmC[d1_l2] * F_mU[2]) - (PmA[d1_l2] * F_mU[2] - PmC[d1_l2] * F_mU[3]);  // p_s1 (d1_l2) - p_s2 (d1_l2)
END
}
print <<"END";
              p_d2_1_d1l2                *= inv_two_zeta;
              p_d2_1_d1l2                += PmB[$d2_l2] * p_p1_d1l2_d2l1 - PmC[$d2_l2] * p_p2_d1l2_d2l1;

              scalar_type d_d0            = del_d1l1_d2l2 * (p_p0_d1l2_d2l1 - p_p1_d1l2_d2l1);
              d_d0                       += del_d1l2_d2l2 * (p_p0_d1l1_d2l1 - p_p1_d1l1_d2l1);
END
if ($d2_l1 == $d2_l2) {
print <<"END";
              d_d0                       += d1_s0 - d1_s1;
END
}
print <<"END";
              d_d0                       *= inv_two_zeta;
              d_d0                       += PmB[$d2_l2] * d1_p0_d2l1 - PmC[$d2_l2] * d1_p1_d2l1;
              scalar_type d_d1            = del_d1l1_d2l2 * (p_p1_d1l2_d2l1 - p_p2_d1l2_d2l1);
              d_d1                       += del_d1l2_d2l2 * (p_p1_d1l1_d2l1 - p_p2_d1l1_d2l1);
END
if ($d2_l1 == $d2_l2) {
print <<"END";
              d_d1                       += d1_s1 - d1_s2;
END
}
print <<"END";
              d_d1                       *= inv_two_zeta;
              d_d1                       += PmB[$d2_l2] * d1_p1_d2l1 - PmC[$d2_l2] * d1_p2_d2l1;
END
#              for (int grad_l = 0; grad_l < 3; grad_l++)
for $grad_l (0..2) {
print <<"END";
              // grad_l = $grad_l
              {
                bool del_d1l1g = d1_l1 == $grad_l, del_d1l2g = d1_l2 == $grad_l;
                C_force_term  = del_d1l1g * p_d2_1_d1l2;
                C_force_term += del_d1l2g * p_d2_1_d1l1;
END
if ($d2_l1 == $grad_l) {
print <<"END";
                C_force_term += d1_p1_d2l2;
END
}
if ($d2_l2 == $grad_l) {
print <<"END";
                C_force_term += d1_p1_d2l1;
END
}
print <<"END";
                C_force_term  = PmC[$grad_l] * d_d1 + inv_two_zeta * C_force_term;
  
                AB_common     = del_d1l1g * p_d2_0_d1l2;
                AB_common    += del_d1l2g * p_d2_0_d1l1;
END
if ($d2_l1 == $grad_l) {
print <<"END";
                AB_common    += d1_p0_d2l2;
END
}
if ($d2_l2 == $grad_l) {
print <<"END";
                AB_common    += d1_p0_d2l1;
END
}
print <<"END";

                A_force_term  = inv_two_zeta * AB_common - C_force_term;
                B_force_term  = PmB[$grad_l] * d_d0 + A_force_term;
                A_force_term += PmA[$grad_l] * d_d0;
                A_force_term *= 2.0f * ai;
                A_force_term -= del_d1l1g * p_d2_0_d1l2;
                A_force_term -= del_d1l2g * p_d2_0_d1l1;
                B_force_term *= 2.0f * aj;
END
if ($d2_l1 == $grad_l) {
print <<"END";
                B_force_term -= d1_p0_d2l2;
END
}
if ($d2_l2 == $grad_l) {
print <<"END";
                B_force_term -= d1_p0_d2l1;
END
}
print <<"END";

                A_force[$grad_l]     += pre_term * A_force_term;
                B_force[$grad_l]     += pre_term * B_force_term;
                C_force[$grad_l][tid]+= valid_thread * prefactor_mm * pre_term * C_force_term;
              }
END
}
print <<"END"
            }
END
}
print <<"END"
          }
END
}
print <<"END";
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-D)-------------------------------------------
END
