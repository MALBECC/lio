//
// QM/MM forces kernel - calculate gradients for QM/MM 1-e operator over
// significant basis primitives
// Each thread maps to a pair of primitives, so each thread contributes partial
// forces on 1 or 2 QM nuclei
// Each thread iterates over every MM atom, so each thread contributes partial
// forces on every MM atom
// The partial forces are calculated using the Obara-Saika recursion relations,
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
__global__ void gpu_qmmm_forces(
    uint num_terms, G2G::vec_type<scalar_type, 2>* ac_values, uint* func2nuc,
    scalar_type* dens_values, uint* func_code, uint* local_dens,
    G2G::vec_type<scalar_type, 3>* mm_forces,
    G2G::vec_type<scalar_type, 3>* qm_forces, uint global_stride,
    G2G::vec_type<scalar_type, 3>* clatom_pos, scalar_type* clatom_chg) {
  uint ffnum = index_x(blockDim, blockIdx, threadIdx);
  int tid = threadIdx.x;
  bool valid_thread = (ffnum < num_terms);

  // Each thread maps to a single pair of QM nuclei, so these forces are
  // computed locally and accumulated at the end
  scalar_type A_force[3] = {0.0, 0.0, 0.0}, B_force[3] = {0.0, 0.0, 0.0};
  uint nuc1, nuc2;
  scalar_type prefactor_qm;

  {
    __shared__ G2G::vec_type<scalar_type, 3>
        clatom_position_sh[QMMM_BLOCK_SIZE];
    __shared__ scalar_type clatom_charge_sh[QMMM_BLOCK_SIZE];
    // Shared memory space for reduction of MM atom force terms
    __shared__ scalar_type C_force[3][QMMM_BLOCK_SIZE];

    scalar_type ai, aj, prefactor_mm, inv_two_zeta;
    scalar_type dens[term_type == 0
                         ? 1
                         : (term_type == 1
                                ? 3
                                : (term_type == 2
                                       ? 9
                                       : (term_type == 3
                                              ? 6
                                              : (term_type == 4 ? 18 : 36))))];
    scalar_type P[3], PmA[3], PmB[3];
    bool same_func = false;

    // uint d1_l1, d1_l2; // These are only used by the d-d kernel; currently
    // this is the only type where each thread maps to a different orbital

    // TODO: each thread calculates its own zeta, overlap, etc here; should
    // these be precalculated and saved (for use here and in Coulomb
    // calculation)?
    {
      scalar_type cc;
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
        // Currently, d-d threads map to a single dxx,dyx,etc for the first
        // function
        // Here we figure out which function (dxx,dyx,dyy,dzx,dzy,or dzz) this
        // thread is doing
        //
        /*uint orb1;
        if (term_type == 5) {
          orb1 = (f1 - gpu_d_offset) % 6;
          switch (orb1) {
            case 0: d1_l1 = 0; d1_l2 = 0; break;
            case 1: d1_l1 = 1; d1_l2 = 0; break;
            case 2: d1_l1 = 1; d1_l2 = 1; break;
            case 3: d1_l1 = 2; d1_l2 = 0; break;
            case 4: d1_l1 = 2; d1_l2 = 1; break;
            case 5: d1_l1 = 2; d1_l2 = 2; break;
          }
          same_func = (f1-orb1) == f2;
        } else {
          same_func = f1 == f2;
        }*/

        //
        // Get the density matrix elements this thread will need
        // The two symmetric cases (p-p and d-d) need to check if function 1 and
        // 2 are the same, and only take the lower triangle
        // of the density block if they are
        //
        uint dens_ind = local_dens[ffnum];
        if (term_type == 2 && same_func) {
          uint true_ind = 0, false_ind = 0;
          for (uint p1 = 0; p1 < 3; p1++) {
            for (uint p2 = 0; p2 < 3; p2++) {
              if (p2 <= p1) {
                dens[false_ind] = dens_values[dens_ind + true_ind];
                true_ind++;
              } else {
                dens[false_ind] = 0.0;
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
                    dens[false_ind] = dens_values[dens_ind + true_ind];
                    true_ind++;
                  } else {
                    dens[false_ind] = 0.0;
                  }
                  false_ind++;
                }
              }
            }
          }
        } else {
          for (uint i = 0; i < TERM_TYPE_GAUSSIANS[term_type]; i++) {
            dens[i] = dens_values[dens_ind + i];
          }
        }

        //
        // Get the function values and nuclei for this thread
        //
        G2G::vec_type<scalar_type, 2> ac1 =
            ac_values[f1 +
                      cont1 * COALESCED_DIMENSION(gpu_m)];  // total_funcs)];
        G2G::vec_type<scalar_type, 2> ac2 =
            ac_values[f2 +
                      cont2 * COALESCED_DIMENSION(gpu_m)];  // total_funcs)];
        ai = ac1.x;
        aj = ac2.x;
        cc = ac1.y * ac2.y;

        nuc1 = func2nuc[f1];
        nuc2 = func2nuc[f2];
      }

      //
      // Precalulate the terms and prefactors that will show up in the forces
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

      prefactor_mm = -cc * 4.0 * PI * ovlap;
      prefactor_qm = prefactor_mm * inv_two_zeta;
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
        // single primitive/primitive overlap force term
        //
        for (int j = 0; j < QMMM_BLOCK_SIZE && i + j < G2G::gpu_atoms; j++) {
          {
            scalar_type PmC[3];
            PmC[0] = P[0] - clatom_position_sh[j].x;
            PmC[1] = P[1] - clatom_position_sh[j].y;
            PmC[2] = P[2] - clatom_position_sh[j].z;
            //
            // Do the core part of the forces calculation - the evaluation of
            // the Obara-Saika recursion equations
            // This is where the different term types differ the most, so these
            // are moved into separate files in the qmmm_terms directory
            // Current version: p-s through d-d are manually unrolled, and d-d
            // is split up over six threads per primitive pair
            //
            // BEGIN TERM-TYPE DEPENDENT PART
            C_force[0][tid] = 0.0;
            C_force[1][tid] = 0.0;
            C_force[2][tid] = 0.0;
            switch (term_type) {
              case 0: {
#include "qmmm_terms/forces/ss.h"
                break;
              }
              case 1: {
#include "qmmm_terms/forces/ps.h"
                break;
              }
              case 2: {
#include "qmmm_terms/forces/pp.h"
                break;
              }
              case 3: {
#include "qmmm_terms/forces/ds.h"
                break;
              }
              case 4: {
#include "qmmm_terms/forces/dp.h"
                break;
              }
              case 5: {
#include "qmmm_terms/forces/dd.h"
                break;
              }
            }
            C_force[0][tid] *= valid_thread * prefactor_mm;
            C_force[1][tid] *= valid_thread * prefactor_mm;
            C_force[2][tid] *= valid_thread * prefactor_mm;
            // END TERM-TYPE DEPENDENT PART
          }

          __syncthreads();

          //
          // BEGIN reduction of MM atom force terms
          //
          // TODO: should we do the per-block reduction here in this loop? or
          // should each thread save its value to global memory for later
          // accumulation?
          //
          // IMPORTANT: ASSUMING BLOCK SIZE OF 128
          //
          // First half of block does x,y
          if (tid < QMMM_FORCES_HALF_BLOCK) {
            C_force[0][tid] += C_force[0][tid + QMMM_FORCES_HALF_BLOCK];
            C_force[1][tid] += C_force[1][tid + QMMM_FORCES_HALF_BLOCK];
          }
          // Second half does z (probably doesn't make much of a difference)
          else {
            C_force[2][tid - QMMM_FORCES_HALF_BLOCK] += C_force[2][tid];
          }
          __syncthreads();
          // first warp does x
          if (tid < WARP_SIZE) {
            warpReduce<scalar_type>(C_force[0], tid);
          }
          // second warp does y
          else if (tid < WARP_SIZE2) {
            warpReduce<scalar_type>(C_force[1], tid - WARP_SIZE);
          }
          // third warp does z
          else if (tid < WARP_SIZE3) {
            warpReduce<scalar_type>(C_force[2], tid - WARP_SIZE2);
          }

          // TODO: tried turning this into one global read to get the force
          // vector object, but didn't seem to improve performance, maybe
          // there's a better way?
          if (tid == 0) {
            qm_forces[global_stride * (i + j) + blockIdx.x].x = C_force[0][0];
          } else if (tid == WARP_SIZE) {
            qm_forces[global_stride * (i + j) + blockIdx.x].y = C_force[1][0];
          } else if (tid == WARP_SIZE2) {
            qm_forces[global_stride * (i + j) + blockIdx.x].z = C_force[2][0];
          }
          //
          // END reduction
          //

          __syncthreads();
        }
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
        // single primitive/primitive overlap force term
        //
        for (int j = 0; j < QMMM_BLOCK_SIZE && i + j < gpu_clatoms; j++) {
          {
            scalar_type PmC[3];
            PmC[0] = P[0] - clatom_position_sh[j].x;
            PmC[1] = P[1] - clatom_position_sh[j].y;
            PmC[2] = P[2] - clatom_position_sh[j].z;
            //
            // Do the core part of the forces calculation - the evaluation of
            // the Obara-Saika recursion equations
            // This is where the different term types differ the most, so these
            // are moved into separate files in the qmmm_terms directory
            // Current version: p-s through d-d are manually unrolled, and d-d
            // is split up over six threads per primitive pair
            //
            // BEGIN TERM-TYPE DEPENDENT PART
            C_force[0][tid] = 0.0;
            C_force[1][tid] = 0.0;
            C_force[2][tid] = 0.0;
            switch (term_type) {
              case 0: {
#include "qmmm_terms/forces/ss.h"
                break;
              }
              case 1: {
#include "qmmm_terms/forces/ps.h"
                break;
              }
              case 2: {
#include "qmmm_terms/forces/pp.h"
                break;
              }
              case 3: {
#include "qmmm_terms/forces/ds.h"
                break;
              }
              case 4: {
#include "qmmm_terms/forces/dp.h"
                break;
              }
              case 5: {
#include "qmmm_terms/forces/dd.h"
                break;
              }
            }
            C_force[0][tid] *= valid_thread * prefactor_mm;
            C_force[1][tid] *= valid_thread * prefactor_mm;
            C_force[2][tid] *= valid_thread * prefactor_mm;
            // END TERM-TYPE DEPENDENT PART
          }

          __syncthreads();

          //
          // BEGIN reduction of MM atom force terms
          //
          // TODO: should we do the per-block reduction here in this loop? or
          // should each thread save its value to global memory for later
          // accumulation?
          //
          // IMPORTANT: ASSUMING BLOCK SIZE OF 128
          //
          // First half of block does x,y
          if (tid < QMMM_FORCES_HALF_BLOCK) {
            C_force[0][tid] += C_force[0][tid + QMMM_FORCES_HALF_BLOCK];
            C_force[1][tid] += C_force[1][tid + QMMM_FORCES_HALF_BLOCK];
          }
          // Second half does z (probably doesn't make much of a difference)
          else {
            C_force[2][tid - QMMM_FORCES_HALF_BLOCK] += C_force[2][tid];
          }
          __syncthreads();
          // first warp does x
          if (tid < WARP_SIZE) {
            warpReduce<scalar_type>(C_force[0], tid);
          }
          // second warp does y
          else if (tid < WARP_SIZE2) {
            warpReduce<scalar_type>(C_force[1], tid - WARP_SIZE);
          }
          // third warp does z
          else if (tid < WARP_SIZE3) {
            warpReduce<scalar_type>(C_force[2], tid - WARP_SIZE2);
          }

          // TODO: tried turning this into one global read to get the force
          // vector object, but didn't seem to improve performance, maybe
          // there's a better way?
          if (tid == 0) {
            mm_forces[global_stride * (i + j) + blockIdx.x].x = C_force[0][0];
          } else if (tid == WARP_SIZE) {
            mm_forces[global_stride * (i + j) + blockIdx.x].y = C_force[1][0];
          } else if (tid == WARP_SIZE2) {
            mm_forces[global_stride * (i + j) + blockIdx.x].z = C_force[2][0];
          }
          //
          // END reduction
          //

          __syncthreads();
        }
      }
    }
  }

  //
  // Reduce the QM force terms
  //
  // TODO: (same question as for the MM forces) - should we do the per-block
  // reduction in this kernel?
  {
    __shared__ bool nuc_flags[MAX_ATOMS];
    __shared__ scalar_type QM_force[3][QMMM_BLOCK_SIZE];

    //
    // First figure out which nuclei are present in this block
    //
    for (int i = 0; i < G2G::gpu_atoms; i += QMMM_BLOCK_SIZE) {
      if (i + tid < G2G::gpu_atoms) nuc_flags[i + tid] = false;
    }
    __syncthreads();
    nuc_flags[nuc1] = true;
    nuc_flags[nuc2] = true;
    __syncthreads();
    for (int i = 0; i < G2G::gpu_atoms; i++) {
      // Only for this block's nuclei
      if (nuc_flags[i] == true) {
        //
        // Load the individual thread's force terms into the appropriate shared
        // location
        //
        bool useA = nuc1 == i, useB = nuc2 == i;
        QM_force[0][tid] = valid_thread * prefactor_qm *
                           (useA * A_force[0] + useB * B_force[0]);
        QM_force[1][tid] = valid_thread * prefactor_qm *
                           (useA * A_force[1] + useB * B_force[1]);
        QM_force[2][tid] = valid_thread * prefactor_qm *
                           (useA * A_force[2] + useB * B_force[2]);
        __syncthreads();

        //
        // Reduce the force terms
        //
        // First half of block does x,y
        if (tid < QMMM_FORCES_HALF_BLOCK) {
          QM_force[0][tid] += QM_force[0][tid + QMMM_FORCES_HALF_BLOCK];
          QM_force[1][tid] += QM_force[1][tid + QMMM_FORCES_HALF_BLOCK];
        }
        // Second half does z
        else {
          QM_force[2][tid - QMMM_FORCES_HALF_BLOCK] += QM_force[2][tid];
        }
        __syncthreads();
        // first warp does x
        if (tid < WARP_SIZE) {
          warpReduce<scalar_type>(QM_force[0], tid);
        }
        // second warp does y
        else if (tid < WARP_SIZE2) {
          warpReduce<scalar_type>(QM_force[1], tid - WARP_SIZE);
        }
        // third warp does z
        else if (tid < WARP_SIZE3) {
          warpReduce<scalar_type>(QM_force[2], tid - WARP_SIZE2);
        }

        if (tid == 0) {
          qm_forces[global_stride * i + blockIdx.x].x += QM_force[0][0];
        } else if (tid == WARP_SIZE) {
          qm_forces[global_stride * i + blockIdx.x].y += QM_force[1][0];
        } else if (tid == WARP_SIZE2) {
          qm_forces[global_stride * i + blockIdx.x].z += QM_force[2][0];
        }
        __syncthreads();
      }
      // At this point, the global QM array is uninitialized; since we'll
      // accumulate all entries, it needs to be zeroed
      // TODO: Zeroing out the partial qm array before this kernel might be
      // better
      // else
      //{
      //  if (tid == 0)       { qm_forces[global_stride*i+blockIdx.x].x = 0.0; }
      //  else if (tid == 32) { qm_forces[global_stride*i+blockIdx.x].y = 0.0; }
      //  else if (tid == 64) { qm_forces[global_stride*i+blockIdx.x].z = 0.0; }
      //}
    }
  }
}
