//
// QM/MM energy kernel - calculate Fock elements for QM/MM 1-e operator over
// significant basis primitives
// Each thread maps to a pair of primitives, and calculates part of each element
// in one block (s-s,p-s,d-d,etc) of the Fock matrix
// Each thread iterates over every MM atom
// The Fock integrals are calculated using the Obara-Saika recursion relations,
// and then reduced per-block
//
// The template parameter term_type defines which type of functions are being
// calculated
// 0 = s-s , 1 = p-s , 2 = p-p, 3 = d-s , 4 = d-p , 5 = d-d
//
// TODO: currently, one thread maps to one primitive-primitive overlap force
// term; is there a better mapping? (thread to function, thread to sub-shell,
// etc)
// TODO: should the loop over MM atoms be broken up to be done by multiple
// blocks rather than a block looping over every MM atom?
//
template <class scalar_type, uint term_type, bool do_cl, bool do_qm>
__global__ void gpu_qmmm_fock(
    uint num_terms, G2G::vec_type<scalar_type, 2>* ac_values, uint* func2nuc,
    uint* func_code, uint* local_fock_ind, double* fock, uint global_stride,
    G2G::vec_type<scalar_type, 3>* clatom_pos, scalar_type* clatom_chg) {
  uint ffnum = index_x(blockDim, blockIdx, threadIdx);
  int tid = threadIdx.x;
  bool valid_thread = (ffnum < num_terms);

  scalar_type prefactor;
  // Each thread accumulates its own energy terms
  // TODO: are these staying on registers or going into local memory? might need
  // to rethink it...
  double my_fock
      [term_type == 0
           ? 1
           : (term_type == 1
                  ? 3
                  : (term_type == 2
                         ? 9
                         : (term_type == 3 ? 6 : (term_type == 4 ? 18 : 36))))];
  for (uint i = 0; i < TERM_TYPE_GAUSSIANS[term_type]; i++) {
    my_fock[i] = 0.0;
  }
  uint fock_ind;
  bool same_func;

  {
    __shared__ G2G::vec_type<scalar_type, 3>
        clatom_position_sh[QMMM_BLOCK_SIZE];
    __shared__ scalar_type clatom_charge_sh[QMMM_BLOCK_SIZE];
    //////////////////////////////

    scalar_type ai, aj, inv_two_zeta;
    scalar_type P[3], PmA[3], PmB[3];

    // TODO: each thread calculates its own zeta, overlap, etc here; should
    // these be precalculated and saved (for use here and in Coulomb
    // calculation)?
    {
      fock_ind = local_fock_ind[ffnum];

      scalar_type cc;
      uint nuc1, nuc2;
      {
        //
        // Decode the function code to figure out which two functions and two
        // primitives this thread maps to
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

        //
        // Get the function values and nuclei for this thread
        //
        G2G::vec_type<scalar_type, 2> ac1 =
            ac_values[f1 + cont1 * COALESCED_DIMENSION(gpu_m)];
        G2G::vec_type<scalar_type, 2> ac2 =
            ac_values[f2 + cont2 * COALESCED_DIMENSION(gpu_m)];
        ai = ac1.x;
        aj = ac2.x;
        cc = ac1.y * ac2.y;

        nuc1 = func2nuc[f1];
        nuc2 = func2nuc[f2];
      }

      //
      // Precalulate the terms and prefactors that will show up in the energy
      // calculation
      //
      scalar_type ovlap;

      G2G::vec_type<scalar_type, 3> A, B;
      // A = G2G::gpu_atom_positions[nuc1];
      A.x = G2G::gpu_atom_positions[nuc1].x;
      A.y = G2G::gpu_atom_positions[nuc1].y;
      A.z = G2G::gpu_atom_positions[nuc1].z;
      // B = G2G::gpu_atom_positions[nuc2];
      B.x = G2G::gpu_atom_positions[nuc2].x;
      B.y = G2G::gpu_atom_positions[nuc2].y;
      B.z = G2G::gpu_atom_positions[nuc2].z;

      //
      // ai and aj can differ by several orders of magnitude
      // They're involved in two additions here, with the results involved in a
      // division
      // Using double precision here is important to maintain precision in the
      // final results
      //
      double zeta = (double)ai + (double)aj;
      inv_two_zeta = 1.0 / (2.0 * zeta);
      P[0] = (A.x * (double)ai + B.x * (double)aj) / zeta;
      P[1] = (A.y * (double)ai + B.y * (double)aj) / zeta;
      P[2] = (A.z * (double)ai + B.z * (double)aj) / zeta;

      PmA[0] = P[0] - A.x;
      PmA[1] = P[1] - A.y;
      PmA[2] = P[2] - A.z;
      PmB[0] = P[0] - B.x;
      PmB[1] = P[1] - B.y;
      PmB[2] = P[2] - B.z;

      G2G::vec_type<scalar_type, 3> AmB = A - B;
      scalar_type ds2 = length2(AmB);
      scalar_type ksi = ((double)ai * (double)aj) / zeta;
      ovlap = exp(-ds2 * ksi);

      prefactor = -2.0 * PI * cc * ovlap / zeta;
    }

    if (do_qm) {
      //
      // Outer loop: read in block of MM atom information into shared memory
      //
      for (int i = 0; i < G2G::gpu_atoms; i += QMMM_BLOCK_SIZE) {
        if (i + tid < G2G::gpu_atoms) {
          clatom_position_sh[tid].x = G2G::gpu_atom_positions[i + tid].x;
          clatom_position_sh[tid].y = G2G::gpu_atom_positions[i + tid].y;
          clatom_position_sh[tid].z = G2G::gpu_atom_positions[i + tid].z;
          clatom_charge_sh[tid] = (scalar_type)(gpu_atom_Z[i + tid]);
        }
        __syncthreads();
        //
        // Inner loop: process block of MM atoms; each thread calculates a
        // single primitive/primitive overlap energy term
        //
        for (int j = 0; j < QMMM_BLOCK_SIZE && i + j < G2G::gpu_atoms; j++) {
          scalar_type PmC[3];
          PmC[0] = P[0] - clatom_position_sh[j].x;
          PmC[1] = P[1] - clatom_position_sh[j].y;
          PmC[2] = P[2] - clatom_position_sh[j].z;
          //
          // Do the core part of the Fock element calculation - the evaluation
          // of the Obara-Saika recursion equations
          // This is where the different term types differ the most, so these
          // are moved into separate files in the qmmm_terms directory
          // Current version: d-p and d-d are manually unrolled
          //
          // BEGIN TERM-TYPE DEPENDENT PART
          switch (term_type) {
            case 0: {
#include "qmmm_terms/energy/ss.h"
              break;
            }
            case 1: {
#include "qmmm_terms/energy/ps.h"
              break;
            }
            case 2: {
#include "qmmm_terms/energy/pp.h"
              break;
            }
            case 3: {
#include "qmmm_terms/energy/ds.h"
              break;
            }
            case 4: {
#include "qmmm_terms/energy/dp.h"
              break;
            }
            case 5: {
#include "qmmm_terms/energy/dd.h"
              break;
            }
          }
          // END TERM-TYPE DEPENDENT PART
        }
        __syncthreads();
      }
    }
    if (do_cl) {
      //
      // Outer loop: read in block of MM atom information into shared memory
      //
      for (int i = 0; i < gpu_clatoms; i += QMMM_BLOCK_SIZE) {
        if (i + tid < gpu_clatoms) {
          clatom_position_sh[tid] = clatom_pos[i + tid];
          clatom_charge_sh[tid] = clatom_chg[i + tid];
        }
        __syncthreads();
        //
        // Inner loop: process block of MM atoms; each thread calculates a
        // single primitive/primitive overlap energy term
        //
        for (int j = 0; j < QMMM_BLOCK_SIZE && i + j < gpu_clatoms; j++) {
          scalar_type PmC[3];
          {
            G2G::vec_type<scalar_type, 3> clatom_pos = clatom_position_sh[j];
            PmC[0] = P[0] - clatom_pos.x;
            PmC[1] = P[1] - clatom_pos.y;
            PmC[2] = P[2] - clatom_pos.z;
          }
          //
          // Do the core part of the Fock element calculation - the evaluation
          // of the Obara-Saika recursion equations
          // This is where the different term types differ the most, so these
          // are moved into separate files in the qmmm_terms directory
          // Current version: d-p and d-d are manually unrolled
          //
          // BEGIN TERM-TYPE DEPENDENT PART
          switch (term_type) {
            case 0: {
#include "qmmm_terms/energy/ss.h"
              break;
            }
            case 1: {
#include "qmmm_terms/energy/ps.h"
              break;
            }
            case 2: {
#include "qmmm_terms/energy/pp.h"
              break;
            }
            case 3: {
#include "qmmm_terms/energy/ds.h"
              break;
            }
            case 4: {
#include "qmmm_terms/energy/dp.h"
              break;
            }
            case 5: {
#include "qmmm_terms/energy/dd.h"
              break;
            }
          }
          // END TERM-TYPE DEPENDENT PART
        }
        __syncthreads();
      }
    }
  }

  if (term_type == 2) {
    if (same_func) {
      uint true_ind = 0, false_ind = 0;
      for (uint p1 = 0; p1 < 3; p1++) {
        for (uint p2 = 0; p2 < 3; p2++) {
          if (p2 <= p1) {
            my_fock[true_ind] = my_fock[false_ind];
            true_ind++;
          }
          false_ind++;
        }
      }
    }
  } else if (term_type == 5) {
    if (same_func) {
      uint true_ind = 0, false_ind = 0;
      for (uint d1_1 = 0; d1_1 < 3; d1_1++) {
        for (uint d1_2 = 0; d1_2 <= d1_1; d1_2++) {
          for (uint d2_1 = 0; d2_1 < 3; d2_1++) {
            for (uint d2_2 = 0; d2_2 <= d2_1; d2_2++) {
              if (!(d2_1 > d1_1 || (d2_1 == d1_1 && d2_2 > d1_2))) {
                my_fock[true_ind] = my_fock[false_ind];
                true_ind++;
              }
              false_ind++;
            }
          }
        }
      }
    }
  }

  //
  // Reduce the partial fock elements
  //
  {
    uint bl_fock_ind = 0;

    __shared__ uint fock_ind_sh[QMMM_BLOCK_SIZE];
    __shared__ bool same_func_sh[QMMM_BLOCK_SIZE];
    __shared__ double fock_sh[QMMM_BLOCK_SIZE];

    //
    // First figure out which fock elements this block contains
    // TODO: there's probably a better way to do this, but this part is tiny
    // compared to the main loop, so...
    //
    fock_ind_sh[tid] = fock_ind;
    __syncthreads();
    uint curr_ind = fock_ind_sh[0], curr_bl_ind = 0;
    //
    // Each thread loops through the block's Fock indices, keeps track of the
    // number of unique indices, and finds its own
    // index relative to the block (bl_fock_ind)
    //
    for (int i = 1; i < QMMM_BLOCK_SIZE; i++) {
      curr_bl_ind += curr_ind != fock_ind_sh[i];
      curr_ind = fock_ind_sh[i];
      bl_fock_ind = (curr_ind == fock_ind) * curr_bl_ind +
                    (curr_ind != fock_ind) * bl_fock_ind;
    }

    //
    // Each thread tells the block which global Fock index corresponds to its
    // own block-local index, and if its Fock index
    // corresponds to an element where function i = function j
    //
    same_func_sh[tid] = false;
    __syncthreads();
    same_func_sh[bl_fock_ind] = same_func;
    fock_ind_sh[bl_fock_ind] = fock_ind;
    __syncthreads();
    //
    // Loop over the Fock elements in this block
    //
    for (int i = 0; i <= curr_bl_ind; i++) {
      bool use_fock = bl_fock_ind == i;
      //
      // For the symmetric cases (p-p and d-d) we need to check if we're looking
      // at a Fock element where function i = function j
      // For these diagonal Fock blocks, only the lower (upper?) triangular part
      // of the block is needed
      //
      uint last_j = TERM_TYPE_GAUSSIANS[term_type];
      if (term_type == 2 && same_func_sh[i]) {
        last_j = PP_SAME_FUNC_SIZE;
      } else if (term_type == 5 && same_func_sh[i]) {
        last_j = DD_SAME_FUNC_SIZE;
      }
      for (int j = 0; j < last_j; j++) {
        //
        // Load the individual thread's Fock terms into the appropriate shared
        // location
        //
        fock_sh[tid] = valid_thread * use_fock * prefactor * my_fock[j];
        __syncthreads();

        //
        // Reduce the Fock terms
        //
        if (tid < QMMM_FORCES_HALF_BLOCK) {
          fock_sh[tid] += fock_sh[tid + QMMM_FORCES_HALF_BLOCK];
        }
        __syncthreads();

        if (tid < WARP_SIZE) {
          warpReduce<double>(fock_sh, tid);
        }
        if (tid == 0) {
          atomicAdd(fock + fock_ind_sh[i] + j, fock_sh[0]);
        }
        __syncthreads();
      }
    }
  }
}
