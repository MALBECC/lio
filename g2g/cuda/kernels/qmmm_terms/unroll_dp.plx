#!/usr/bin/perl -w

print <<"END";
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NOTE: THIS FILE WAS GENERATED AUTOMATICALLY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//-------------------------------------------BEGIN TERM-TYPE DEPENDENT PART (D-P)-------------------------------------------
        scalar_type F_mU[5];
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * zeta;
          //for (int m = 0; m <= 4; m++) 
          //{
          //  F_mU[m] = lio_gamma<scalar_type>(m,U);
          lio_gamma<scalar_type,4>(F_mU,U);
          //}
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        {
          C_force[0][tid] = 0.0f; C_force[1][tid] = 0.0f; C_force[2][tid] = 0.0f;
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type AB_common;
          //scalar_type mm_charge = clatom_charge_sh[j];
          uint dens_ind = 0;

END
#          for (int d_l1 = 0; d_l1 < 3; d_l1++)
for $d_l1 (0..2) {
print <<"END";
          {
END
#            for (int d_l2 = 0; d_l2 <= d_l1; d_l2++)
for $d_l2 (0..$d_l1) {
print <<"END";
            {

              scalar_type d_s0  = PmA[$d_l1] * (PmA[$d_l2] * F_mU[0] - PmC[$d_l2] * F_mU[1]); // p_s0 (d_l2)
              d_s0             -= PmC[$d_l1] * (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]); // p_s1 (d_l2)
              scalar_type d_s1  = PmA[$d_l1] * (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]); // p_s1 (d_l2)
              d_s1             -= PmC[$d_l1] * (PmA[$d_l2] * F_mU[2] - PmC[$d_l2] * F_mU[3]); // p_s2 (d_l2)
              scalar_type d_s2  = PmA[$d_l1] * (PmA[$d_l2] * F_mU[2] - PmC[$d_l2] * F_mU[3]); // p_s2 (d_l2)
              d_s2             -= PmC[$d_l1] * (PmA[$d_l2] * F_mU[3] - PmC[$d_l2] * F_mU[4]); // p_s3 (d_l2)
END
if ($d_l1 == $d_l2) {
print <<"END";
              d_s0             += inv_two_zeta * (F_mU[0] - F_mU[1]);
              d_s1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
              d_s2             += inv_two_zeta * (F_mU[2] - F_mU[3]);
END
}
#              for (int p_l = 0; p_l < 3; p_l++)
for $p_l (0..2) {
print <<"END";
              {
END
if ($d_l1 == $d_l2) {
print <<"END";
                scalar_type pre_term = gpu_normalization_factor*clatom_charge_sh[j] * dens[dens_ind];
END
} else {
print <<"END";
                scalar_type pre_term = clatom_charge_sh[j] * dens[dens_ind];
END
}
print <<"END";
                dens_ind++;

                scalar_type p_p0_d1  = PmB[$p_l] * (PmA[$d_l1] * F_mU[0] - PmC[$d_l1] * F_mU[1]); // p_s0 (d_l1)
                p_p0_d1             -= PmC[$p_l] * (PmA[$d_l1] * F_mU[1] - PmC[$d_l1] * F_mU[2]); // p_s1 (d_l1)
                scalar_type p_p1_d1  = PmB[$p_l] * (PmA[$d_l1] * F_mU[1] - PmC[$d_l1] * F_mU[2]); // p_s1 (d_l1)
                p_p1_d1             -= PmC[$p_l] * (PmA[$d_l1] * F_mU[2] - PmC[$d_l1] * F_mU[3]); // p_s2 (d_l1)
END
if ($d_l1 == $p_l) {
print <<"END";
                p_p0_d1             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d1             += inv_two_zeta * (F_mU[1] - F_mU[2]);
END
}
print <<"END";
                scalar_type p_p0_d2  = PmB[$p_l] * (PmA[$d_l2] * F_mU[0] - PmC[$d_l2] * F_mU[1]); // p_s0 (d_l2)
                p_p0_d2             -= PmC[$p_l] * (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]); // p_s1 (d_l2)
                scalar_type p_p1_d2  = PmB[$p_l] * (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]); // p_s1 (d_l2)
                p_p1_d2             -= PmC[$p_l] * (PmA[$d_l2] * F_mU[2] - PmC[$d_l2] * F_mU[3]); // p_s2 (d_l2)
END
if ($d_l2 == $p_l) {
print <<"END";
                p_p0_d2             += inv_two_zeta * (F_mU[0] - F_mU[1]);
                p_p1_d2             += inv_two_zeta * (F_mU[1] - F_mU[2]);
END
}
print <<"END";
                scalar_type d_p0 = 0.0f,d_p1 = 0.0f;
END
if ($d_l1 == $p_l) {
print <<"END";
                d_p0  = (PmA[$d_l2] * F_mU[0] - PmC[$d_l2] * F_mU[1]) - (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]);  // p_s0 (d_l2) - p_s1 (d_l2)
END
}
if ($d_l2 == $p_l) {
print <<"END";
                d_p0 += (PmA[$d_l1] * F_mU[0] - PmC[$d_l1] * F_mU[1]) - (PmA[$d_l1] * F_mU[1] - PmC[$d_l1] * F_mU[2]);  // p_s0 (d_l1) - p_s1 (d_l1)
END
}
print <<"END";
                d_p0 *= inv_two_zeta;
                d_p0 += PmB[$p_l] * d_s0 - PmC[$p_l] * d_s1;
END
if ($d_l1 == $p_l) {
print <<"END";
                d_p1  = (PmA[$d_l2] * F_mU[1] - PmC[$d_l2] * F_mU[2]) - (PmA[$d_l2] * F_mU[2] - PmC[$d_l2] * F_mU[3]);  // p_s1 (d_l2) - p_s2 (d_l2)
END
}
if ($d_l2 == $p_l) {
print <<"END";
                d_p1 += (PmA[$d_l1] * F_mU[1] - PmC[$d_l1] * F_mU[2]) - (PmA[$d_l1] * F_mU[2] - PmC[$d_l1] * F_mU[3]);  // p_s1 (d_l1) - p_s2 (d_l1)
END
}
print <<"END";
                d_p1 *= inv_two_zeta;
                d_p1 += PmB[$p_l] * d_s1 - PmC[$p_l] * d_s2;
END
#                for (int grad_l = 0; grad_l < 3; grad_l++)
for $grad_l (0..2) {
print <<"END";
                {
                  C_force_term  = 0.0f;
                  AB_common = 0.0f;
END
if ($d_l1 == $grad_l) {
print <<"END";
                  C_force_term += p_p1_d2;
                  AB_common    += p_p0_d2;
END
}
if ($d_l2 == $grad_l) {
print <<"END";
                  C_force_term += p_p1_d1;
                  AB_common    += p_p0_d1;
END
}
if ($p_l == $grad_l) {
print <<"END";
                  C_force_term += d_s1;
                  AB_common    += d_s0;
END
}
print <<"END";
                  C_force_term  = PmC[$grad_l] * d_p1 + inv_two_zeta * C_force_term;
  
                  A_force_term  = inv_two_zeta * AB_common - C_force_term;
                  B_force_term  = PmB[$grad_l] * d_p0 + A_force_term;
                  A_force_term += PmA[$grad_l] * d_p0;
                  A_force_term *= 2.0f * ai;
END
if ($d_l1 == $grad_l) {
print <<"END";
                  A_force_term -= p_p0_d2;
END
}
if ($d_l2 == $grad_l) {
print <<"END";
                  A_force_term -= p_p0_d1;
END
}
print <<"END";
                  B_force_term *= 2.0f * aj;
END
if ($p_l == $grad_l) {
print <<"END";
                  B_force_term -= d_s0;
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
END
}
print <<"END";
        }
        // END individual force terms
//-------------------------------------------END TERM-TYPE DEPENDENT PART (D-P)-------------------------------------------
END
