//
// QM/MM forces kernel - calculate gradients for QM/MM 1-e operator over significant basis primitives
// Each thread maps to a pair of primitives, so each thread contributes partial forces on 1 or 2 QM nuclei
// Each thread iterates over every MM atom, so each thread contributes partial forces on every MM atom
// The partial forces are calculated using the Obara-Saika recursion relations, and then reduced per-block
//
// The template parameter term_type defines which type of functions are being calculated
// 0 = s-s , 1 = p-s , 2 = p-p, 3 = d-s , 4 = d-p , 5 = d-d
//
// TODO: currently, one thread maps to one primitive-primitive overlap force term; is there a better mapping? (thread to function, thread to sub-shell, etc)
// TODO: should the loop over MM atoms be broken up to be done by multiple blocks rather than a block looping over every MM atom?
//
template<class scalar_type, uint term_type>
__global__ void gpu_coulomb_fit1( uint num_terms, G2G::vec_type<scalar_type,2>* ac_values, uint* func2nuc, scalar_type* dens_values, uint* func_code, uint* local_dens,
                                    double* rc_partial, uint global_stride, G2G::vec_type<scalar_type,2>* ac_values_dens, G2G::vec_type<scalar_type,3>* nuc_pos_dens,
                                    uint s_end, uint p_end, uint d_end, uint p_offset, uint d_offset )
{

  uint ffnum = index_x(blockDim, blockIdx, threadIdx);
  int tid = threadIdx.x;
  bool valid_thread = (ffnum < num_terms);// && term_type == 1;

  scalar_type prefactor_mo;
  // Each thread accumulates its own energy terms
  // TODO: are these staying on registers or going into local memory? might need to rethink it...
  /*double my_fock[term_type==0? 1 : (term_type==1? 3 : (term_type==2? 9 : (term_type==3? 6 : (term_type==4? 18 : 36))))];
  for (uint i = 0; i < TERM_TYPE_GAUSSIANS[term_type]; i++) {
    my_fock[i] = 0.0f;
  }
  uint fock_ind;*/
  bool same_func;

  {
    __shared__ G2G::vec_type<scalar_type,3> nuc_pos_dens_sh[QMMM_BLOCK_SIZE];
    __shared__ G2G::vec_type<scalar_type,2> ac_val_dens_sh[QMMM_BLOCK_SIZE];
    //__shared__ scalar_type fit_dens_sh[QMMM_BLOCK_SIZE];
    __shared__ double rc_sh[6][QMMM_BLOCK_SIZE];

    scalar_type ai, aj, inv_two_zeta;
    scalar_type dens[term_type==0? 1 : (term_type==1? 3 : (term_type==2? 9 : (term_type==3? 6 : (term_type==4? 18 : 36))))];
    scalar_type P[3], PmA[3], PmB[3];

    // TODO: each thread calculates its own zeta, overlap, etc here; should these be precalculated and saved (for use here and in Coulomb calculation)?
    {
      //fock_ind = local_fock_ind[ffnum];

      scalar_type cc;
      uint nuc1, nuc2;
      {
        //
        // Decode the function code to figure out which two functions and two primitives this thread maps to
        //
        uint my_func_code = func_code[ffnum];

        uint div = MAX_CONTRACTIONS;
        uint cont2 = my_func_code % div;
        my_func_code /= div;
        uint cont1 = my_func_code % div;
        my_func_code /= div;

        div = gpu_m;
        uint f2 = my_func_code % div;
        my_func_code /= div;
        uint f1 = my_func_code;

        same_func = f1 == f2;
        //valid_thread = valid_thread && f1==62 && f2<=1;
        uint dens_ind = local_dens[ffnum];
        if (term_type == 2 && same_func) {
          uint true_ind = 0, false_ind = 0;
          for (uint p1 = 0; p1 < 3; p1++) {
            for (uint p2 = 0; p2 < 3; p2++) {
              if (p2 <= p1) {
                dens[false_ind] = dens_values[dens_ind+true_ind];
                true_ind++;
              } else {
                dens[false_ind] = 0.0f;
              }
              false_ind++;
            }
          }
        } else if (term_type == 5 && same_func) {
          uint true_ind = 0, false_ind = 0;
          for (uint d1_1 = 0; d1_1 < 3; d1_1++) {
            for (uint d1_2 = 0; d1_2 <= d1_1; d1_2++) {
              for (uint d2_1 = 0; d2_1 < 3; d2_1++) {
                for (uint d2_2 = 0; d2_2 <= d2_1; d2_2++) {
                  if (!(d2_1 > d1_1 || (d2_1 == d1_1 && d2_2 > d1_2))) {
                    dens[false_ind] = dens_values[dens_ind+true_ind];
                    true_ind++;
                  } else {
                    dens[false_ind] = 0.0f;
                  }
                  false_ind++;
                }
              }
            }
          }
        } else {
          for (uint i = 0; i < TERM_TYPE_GAUSSIANS[term_type]; i++) {
            dens[i] = dens_values[dens_ind+i];
          }
        }


        //
        // Get the function values and nuclei for this thread
        //
        G2G::vec_type<scalar_type,2> ac1 = ac_values[f1 + cont1 * COALESCED_DIMENSION(gpu_m)];//total_funcs)];
        G2G::vec_type<scalar_type,2> ac2 = ac_values[f2 + cont2 * COALESCED_DIMENSION(gpu_m)];//total_funcs)];
        ai = ac1.x;
        aj = ac2.x;
        cc = ac1.y * ac2.y;

        nuc1 = func2nuc[f1];
        nuc2 = func2nuc[f2];
      }

      //
      // Precalulate the terms and prefactors that will show up in the forces calculation
      //
      scalar_type ovlap;
  
      G2G::vec_type<scalar_type,3> A, B;
      //A = G2G::gpu_atom_positions[nuc1];
      A.x = G2G::gpu_atom_positions[nuc1].x;
      A.y = G2G::gpu_atom_positions[nuc1].y;
      A.z = G2G::gpu_atom_positions[nuc1].z;
      //B = G2G::gpu_atom_positions[nuc2];
      B.x = G2G::gpu_atom_positions[nuc2].x;
      B.y = G2G::gpu_atom_positions[nuc2].y;
      B.z = G2G::gpu_atom_positions[nuc2].z;
  
      //
      // ai and aj can differ by several orders of magnitude
      // They're involved in two additions here, with the results involved in a division
      // Using double precision here is important to maintain precision in the final results
      //
      double zeta = (double)ai + (double)aj;
      inv_two_zeta = 1.0 / (2.0 * zeta);
      P[0] = (A.x*(double)ai + B.x*(double)aj) / zeta;
      P[1] = (A.y*(double)ai + B.y*(double)aj) / zeta;
      P[2] = (A.z*(double)ai + B.z*(double)aj) / zeta;

      PmA[0] = P[0] - A.x;
      PmA[1] = P[1] - A.y;
      PmA[2] = P[2] - A.z;
      PmB[0] = P[0] - B.x;
      PmB[1] = P[1] - B.y;
      PmB[2] = P[2] - B.z;

      G2G::vec_type<scalar_type,3> AmB = A - B;
      scalar_type ds2 = length2(AmB);
      scalar_type ksi = ((double)ai*(double)aj)/zeta;
      ovlap = exp(-ds2*ksi);

      prefactor_mo = (double)(cc * 2.0f * PI52 * ovlap) / zeta;
    }
    __shared__ uint term_start[3];
    term_start[0] = 0; term_start[1] = p_offset; term_start[2] = d_offset;
    __shared__ uint term_end[3];
    term_end[0] = s_end; term_end[1] = p_end; term_end[2] = d_end;
    __shared__ uint inner_stop[3];
    inner_stop[0] = QMMM_BLOCK_SIZE; inner_stop[1] = 126; inner_stop[2] = 126;
    __shared__ uint inner_step[3];
    inner_step[0] = 1; inner_step[1] = 3; inner_step[2] = 6;

    uint rc_ind = 0;
    #pragma unroll 3
    for (int func_type = 0; func_type < 3; func_type++) {
      //
      // Outer loop: read in block of MM atom information into shared memory
      //
      for (int i = term_start[func_type]; i < term_end[func_type]; i += QMMM_BLOCK_SIZE)
      {
        if (i + tid < term_end[func_type]) {
          nuc_pos_dens_sh[tid] = nuc_pos_dens[i+tid];
          ac_val_dens_sh[tid] = ac_values_dens[i+tid];
          //fit_dens_sh[tid] = fit_dens[i+tid];
        }
        __syncthreads();
        //
        // Inner loop: process block of MM atoms; each thread calculates a single primitive/primitive overlap force term
        //
        for (int j = 0; j < inner_stop[func_type] && i+j < term_end[func_type]; j += inner_step[func_type])
        {
          {
            scalar_type WmP[3], WmQ[3], inv_two_zeta_eta, rho, rho_zeta, prefactor_dens;
            scalar_type inv_two_eta, rho_eta;
            {
              double zeta = (double)ai + (double)aj;
              double zeta_eta = zeta + (double)ac_val_dens_sh[j].x;
              scalar_type W[3];
              W[0] = (P[0]*zeta + nuc_pos_dens_sh[j].x*(double)ac_val_dens_sh[j].x) / zeta_eta;
              W[1] = (P[1]*zeta + nuc_pos_dens_sh[j].y*(double)ac_val_dens_sh[j].x) / zeta_eta;
              W[2] = (P[2]*zeta + nuc_pos_dens_sh[j].z*(double)ac_val_dens_sh[j].x) / zeta_eta;
              WmP[0] = W[0] - P[0];
              WmP[1] = W[1] - P[1];
              WmP[2] = W[2] - P[2];
              WmQ[0] = W[0] - nuc_pos_dens_sh[j].x;
              WmQ[1] = W[1] - nuc_pos_dens_sh[j].y;
              WmQ[2] = W[2] - nuc_pos_dens_sh[j].z;
              rho = ((double)ac_val_dens_sh[j].x * zeta) / zeta_eta;
              inv_two_zeta_eta = 1.0 / (2.0 * zeta_eta);
              rho_zeta = (double)ac_val_dens_sh[j].x / zeta_eta;
              inv_two_eta = 1.0 / (2.0 * (double)ac_val_dens_sh[j].x);
              rho_eta = zeta / zeta_eta;
              prefactor_dens = (double)ac_val_dens_sh[j].y / ((double)ac_val_dens_sh[j].x * sqrt(zeta_eta));
            }
            //
            // Do the core part of the forces calculation - the evaluation of the Obara-Saika recursion equations
            // This is where the different term types differ the most, so these are moved into separate files in the qmmm_terms directory
            // Current version: p-s through d-d are manually unrolled, and d-d is split up over six threads per primitive pair
            //
            // BEGIN TERM-TYPE DEPENDENT PART
            rc_sh[0][tid] = 0.0; rc_sh[1][tid] = 0.0; rc_sh[2][tid] = 0.0;
            rc_sh[3][tid] = 0.0; rc_sh[4][tid] = 0.0; rc_sh[5][tid] = 0.0;
#ifdef FOCK_CALC
#undef FOCK_CALC
#endif
            switch (term_type)
            {
              case 0:
                switch (func_type)
                {
                  case 0:
                    #include "coulomb_terms/energy/ss_s.h"
                    break;
                  case 1:
                    #include "coulomb_terms/energy/ss_p.h"
                    break;
                  case 2:
                    #include "coulomb_terms/energy/ss_d.h"
                    break;
                }
                break;
              case 1:
                switch (func_type)
                {
                  case 0:
                    #include "coulomb_terms/energy/ps_s.h"
                    break;
                  case 1:
                    #include "coulomb_terms/energy/ps_p.h"
                    break;
                  case 2:
                    #include "coulomb_terms/energy/ps_d.h"
                    break;
                }
                break;
              case 2:
                switch (func_type)
                {
                  case 0:
                    #include "coulomb_terms/energy/pp_s.h"
                    break;
                  case 1:
                    #include "coulomb_terms/energy/pp_p.h"
                    break;
                  case 2:
                    #include "coulomb_terms/energy/pp_d.h"
                    break;
                }
                break;
              case 3:
                switch (func_type)
                {
                  case 0:
                    #include "coulomb_terms/energy/ds_s.h"
                    break;
                  case 1:
                    #include "coulomb_terms/energy/ds_p.h"
                    break;
                  case 2:
                    #include "coulomb_terms/energy/ds_d.h"
                    break;
                }
                break;
              case 4:
                switch (func_type)
                {
                  case 0:
                    #include "coulomb_terms/energy/dp_s.h"
                    break;
                  case 1:
                    #include "coulomb_terms/energy/dp_p.h"
                    break;
                  case 2:
                    #include "coulomb_terms/energy/dp_d.h"
                    break;
                }
                break;
              case 5:
                switch (func_type)
                {
                  case 0:
                    #include "coulomb_terms/energy/dd_s.h"
                    break;
                  case 1:
                    #include "coulomb_terms/energy/dd_p.h"
                    break;
                  case 2:
                    #include "coulomb_terms/energy/dd_d.h"
                    break;
                }
                break;
            }
            // END TERM-TYPE DEPENDENT PART
          }
          switch(func_type)
          {
            case 0:
              rc_sh[0][tid] = valid_thread * prefactor_mo * rc_sh[0][tid];

              __syncthreads();
              if (tid < QMMM_FORCES_HALF_BLOCK)
              {
                rc_sh[0][tid] += rc_sh[0][tid+QMMM_FORCES_HALF_BLOCK];
              }
              __syncthreads();
              if (tid < WARP_SIZE) { warpReduce<double>(rc_sh[0], tid); }
              if (tid == 0)        { rc_partial[global_stride*blockIdx.x+rc_ind] += rc_sh[0][0]; }
              break;
            case 1:
              rc_sh[0][tid] = valid_thread * prefactor_mo * rc_sh[0][tid];
              rc_sh[1][tid] = valid_thread * prefactor_mo * rc_sh[1][tid];
              rc_sh[2][tid] = valid_thread * prefactor_mo * rc_sh[2][tid];
              __syncthreads();
              if (tid < QMMM_FORCES_HALF_BLOCK)
              {
                rc_sh[0][tid] += rc_sh[0][tid+QMMM_FORCES_HALF_BLOCK];
                rc_sh[1][tid] += rc_sh[1][tid+QMMM_FORCES_HALF_BLOCK];
              }
              else
              {
                rc_sh[2][tid-QMMM_FORCES_HALF_BLOCK] += rc_sh[2][tid];
              }
              __syncthreads();
              if (tid < WARP_SIZE)        { warpReduce<double>(rc_sh[0], tid); }
              else if (tid < WARP_SIZE2)  { warpReduce<double>(rc_sh[1], tid-WARP_SIZE); }
              else if (tid < WARP_SIZE3)  { warpReduce<double>(rc_sh[2], tid-WARP_SIZE2); }
              if (tid == 0)               { rc_partial[global_stride*blockIdx.x+rc_ind] += rc_sh[0][0]; }
              else if (tid == WARP_SIZE)  { rc_partial[global_stride*blockIdx.x+(rc_ind+1)] += rc_sh[1][0]; }
              else if (tid == WARP_SIZE2) { rc_partial[global_stride*blockIdx.x+(rc_ind+2)] += rc_sh[2][0]; }
              break;
            case 2:
              rc_sh[0][tid] = valid_thread * prefactor_mo * rc_sh[0][tid];
              rc_sh[1][tid] = valid_thread * prefactor_mo * rc_sh[1][tid];
              rc_sh[2][tid] = valid_thread * prefactor_mo * rc_sh[2][tid];
              rc_sh[3][tid] = valid_thread * prefactor_mo * rc_sh[3][tid];
              rc_sh[4][tid] = valid_thread * prefactor_mo * rc_sh[4][tid];
              rc_sh[5][tid] = valid_thread * prefactor_mo * rc_sh[5][tid];
              __syncthreads();
              if (tid < QMMM_FORCES_HALF_BLOCK)
              {
                rc_sh[0][tid] += rc_sh[0][tid+QMMM_FORCES_HALF_BLOCK];
                rc_sh[1][tid] += rc_sh[1][tid+QMMM_FORCES_HALF_BLOCK];
                rc_sh[2][tid] += rc_sh[2][tid+QMMM_FORCES_HALF_BLOCK];
              }
              else
              {
                rc_sh[3][tid-QMMM_FORCES_HALF_BLOCK] += rc_sh[3][tid];
                rc_sh[4][tid-QMMM_FORCES_HALF_BLOCK] += rc_sh[4][tid];
                rc_sh[5][tid-QMMM_FORCES_HALF_BLOCK] += rc_sh[5][tid];
              }
              __syncthreads();
              if (tid < WARP_SIZE)        { warpReduce<double>(rc_sh[0], tid); 
                                            warpReduce<double>(rc_sh[1], tid); }
              else if (tid < WARP_SIZE2)  { warpReduce<double>(rc_sh[2], tid-WARP_SIZE);
                                            warpReduce<double>(rc_sh[3], tid-WARP_SIZE); }
              else if (tid < WARP_SIZE3)  { warpReduce<double>(rc_sh[4], tid-WARP_SIZE2);
                                            warpReduce<double>(rc_sh[5], tid-WARP_SIZE2); }
              if (tid == 0)               { rc_partial[global_stride*blockIdx.x+(rc_ind)] += rc_sh[0][0];
                                            rc_partial[global_stride*blockIdx.x+(rc_ind+1)] += rc_sh[1][0]; }
              else if (tid == WARP_SIZE)  { rc_partial[global_stride*blockIdx.x+(rc_ind+2)] += rc_sh[2][0]; 
                                            rc_partial[global_stride*blockIdx.x+(rc_ind+3)] += rc_sh[3][0]; }
              else if (tid == WARP_SIZE2) { rc_partial[global_stride*blockIdx.x+(rc_ind+4)] += rc_sh[4][0]; 
                                            rc_partial[global_stride*blockIdx.x+(rc_ind+5)] += rc_sh[5][0]; }
              break;
          }
          rc_ind += inner_step[func_type];
          __syncthreads();
        }
      }
    }
  }
}

//
// Reduce the partial Fock array into the first row
// Also, reduce corresponding partial energies over the block
//
template<class scalar_type>
__global__ void gpu_coulomb_rc_term_reduce( double* rc_partial,uint stride, uint start, uint end, uint width )
{

  uint my_rc_ind = index_x(blockDim, blockIdx, threadIdx);
  double my_partial_rc = 0.0;

  if (my_rc_ind < width) {
    for (uint i = start; i < end; i++) {
      my_partial_rc += rc_partial[stride*i + my_rc_ind];
    }
    rc_partial[stride*start+my_rc_ind] = my_partial_rc;
  }

}

template<class scalar_type>
__global__ void gpu_coulomb_rc_reduce( double* rc_partial,uint stride,uint width,uint num_terms )
{

  uint my_rc_ind = index_x(blockDim, blockIdx, threadIdx);
  double my_partial_rc = 0.0;

  if (my_rc_ind < width) {
    for (uint i = 0; i < num_terms; i++) {
      my_partial_rc += rc_partial[stride*gpu_out_offsets[i] + my_rc_ind];
    }
    rc_partial[my_rc_ind] = my_partial_rc;
  }

}

template<class scalar_type>
__global__ void gpu_coulomb_fit2( double* rc, scalar_type* fit_dens, double* G_inv, uint* input_ind, uint m_dens )
{

  uint ffnum = index_x(blockDim, blockIdx, threadIdx);
  bool valid_thread = ffnum < m_dens;
  uint af_ind = input_ind[ffnum*valid_thread];
  int tid = threadIdx.x;

  __shared__ double rc_sh[QMMM_BLOCK_SIZE];
  double my_fit_dens = 0.0;


  for (int i = 0; i < m_dens; i += QMMM_BLOCK_SIZE)
  {
    if (i + tid < m_dens)
    {
      rc_sh[tid] = rc[i+tid];
    }
    __syncthreads();
    for (int j = 0; j < QMMM_BLOCK_SIZE && i+j < m_dens; j++)
    {
      my_fit_dens += rc_sh[j] * G_inv[(i+j)*COALESCED_DIMENSION(m_dens)+ffnum*valid_thread];
    }
    __syncthreads();
  }
  fit_dens[af_ind] = my_fit_dens*valid_thread;
}

