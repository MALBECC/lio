#include <fstream>
#include <iostream>
#include <vector>

#include "../common.h"
#include "../init.h"
#include "../matrix.h"
#include "../timer.h"
#include "../scalar_vector_types.h"

#include "aint_init.h"
#include "aint_common.h"
#include "os_integral.h"

using std::cout;
using std::vector;
using std::endl;

namespace AINT {

uint NUM_TERM_TYPES;

//
// Do check between all basis primitives to find those with significant overlap
// Check the resulting Gaussian argument from two primitives to the rmax
// parameter; only use primitives within that cut-off
//
// A single thread gets mapped to a pair of significant primitives
// We set up here arrays that tell which two functions/two primitives a thread
// is calculating
// We also pick out the density matrix elements for significant functions here
//
template <class scalar_type>
void OSIntegral<scalar_type>::new_cutoff(void) {
#if GPU_KERNELS
  clear();
#endif

  uint i, j, ni, nj;
  uint i_orbitals, j_orbitals;
  uint nuc_i, nuc_j;
  G2G::vec_type<double, 3> A, B, AmB;
  double ai, aj;
  double dsq, ksi, zeta;
  uint num_terms = 0, total_num_terms = 0;

  this->term_type_offsets[0] = 0;  // s-s starts at 0
  uint i_begin, i_end, j_begin, j_end;
  uint tmp_ind = 0;
  uint s_start = 0, p_start = G2G::fortran_vars.s_funcs,
       d_start = G2G::fortran_vars.s_funcs + G2G::fortran_vars.p_funcs * 3,
       m = G2G::fortran_vars.m;

  //                                       s-s        p-s        p-p        d-s
  //                                       d-p        d-d
  uint i_begin_vals[MAX_TERM_TYPES] = {s_start, p_start, p_start,
                                       d_start, d_start, d_start};
  uint i_end_vals[MAX_TERM_TYPES] = {p_start, d_start, d_start, m, m, m};
  uint j_begin_vals[MAX_TERM_TYPES] = {s_start, s_start, p_start,
                                       s_start, p_start, d_start};
  uint j_end_vals[MAX_TERM_TYPES] = {p_start - 1, p_start - 1, d_start - 1,
                                     p_start - 1, d_start - 1, m - 1};
  uint i_orbital_vals[MAX_TERM_TYPES] = {1, 3, 3, 6, 6, 6};
  uint j_orbital_vals[MAX_TERM_TYPES] = {1, 1, 3, 1, 3, 6};

  uint local_dens_ind, num_dens_terms = 0, total_dens_terms = 0;
  this->dens_offsets[0] = 0;
  uint tmp_dens_ind = 0;

  if (p_start == m) {         // When there are only s functions.
    NUM_TERM_TYPES = 1;
  } else if (d_start == m) {  // When there are only s and p functions.
    NUM_TERM_TYPES = 3;
  } else {
    NUM_TERM_TYPES = 6;
  }

  for (uint current_term_type = 0; current_term_type < NUM_TERM_TYPES;
       current_term_type++) {
    this->term_type_counts[current_term_type] = 0;
    i_begin = i_begin_vals[current_term_type];
    i_end = i_end_vals[current_term_type];
    j_begin = j_begin_vals[current_term_type];
    j_end = j_end_vals[current_term_type];
    i_orbitals = i_orbital_vals[current_term_type];
    j_orbitals = j_orbital_vals[current_term_type];

    this->dens_counts[current_term_type] = 0;
    local_dens_ind = 0;

    // We pad the input arrays between term types, so the offsets for each term
    // type need to be tracked
    if (current_term_type > 0) {
      tmp_ind += term_type_counts[current_term_type - 1];
      tmp_ind += QMMM_BLOCK_SIZE -
                 (term_type_counts[current_term_type - 1] % QMMM_BLOCK_SIZE);
      this->term_type_offsets[current_term_type] = tmp_ind;
      tmp_dens_ind += dens_counts[current_term_type - 1];
      tmp_dens_ind += QMMM_BLOCK_SIZE -
                      (dens_counts[current_term_type - 1] % QMMM_BLOCK_SIZE);
      this->dens_offsets[current_term_type] = tmp_dens_ind;
    }

    // function i, center A
    for (i = i_begin; i < i_end; i += i_orbitals) {
      nuc_i = G2G::fortran_vars.nucleii(i) - 1;
      A = G2G::fortran_vars.atom_positions(nuc_i);
      // function j, center B
      for (j = j_begin; j <= ((i > j_end) ? j_end : i); j += j_orbitals) {
        nuc_j = G2G::fortran_vars.nucleii(j) - 1;
        B = G2G::fortran_vars.atom_positions(nuc_j);
        AmB = A - B;
        dsq = length2(AmB);
        bool use_funcs = false;  // Do these two functions have any significant
                                 // primitive pairs?
        // primitive ni, function i
        for (ni = 0; ni < G2G::fortran_vars.contractions(i); ni++) {
          // primitive nj, function j
          for (nj = 0; nj < G2G::fortran_vars.contractions(j); nj++) {
            ai = G2G::fortran_vars.a_values(i, ni);
            aj = G2G::fortran_vars.a_values(j, nj);
            zeta = ai + aj;
            ksi = ai * aj / zeta;
            total_num_terms += (i == j) ? i_orbitals * (i_orbitals + 1) / 2
                                        : i_orbitals * j_orbitals;

            if (!this->apply_cutoff || dsq * ksi < integral_vars.rmax) {
              use_funcs = true;
              num_terms++;
              this->term_type_counts[current_term_type]++;

              // Encoding which two primitives this thread will calculate a
              // force term for in one number
              // NOTE: with a long integer, we can use this scheme up to an m of
              // about 9000
              //       if m needs to ever be larger than that, we need to break
              //       this up into multiple arrays
              uint this_func_code =
                  nj;  // First, primitive # nj in the lowest part
              this_func_code +=
                  ni *
                  MAX_CONTRACTIONS;  // Primitive # ni after the space for nj
              this_func_code +=
                  j * MAX_CONTRACTIONS *
                  MAX_CONTRACTIONS;  // Function # j after space for primitives
              this_func_code +=
                  i * MAX_CONTRACTIONS * MAX_CONTRACTIONS *
                  G2G::fortran_vars
                      .m;  // Finally, function # i in the highest part

              this->func_code.push_back(
                  this_func_code);  // Which primitives the thread represents
              this->local_dens.push_back(local_dens_ind);  // Which part of the
                                                           // (reduced) density
                                                           // matrix the thread
                                                           // needs
            }
          }
        }

        total_dens_terms += (i == j) ? i_orbitals * (i_orbitals + 1) / 2
                                     : i_orbitals * j_orbitals;
        // dens_values is a reduced density matrix that only keeps the elements
        // of functions with significant primitive pairs
        // local2globaldens maps from a reduced density/Fock index back to the
        // full matrix
        if (use_funcs) {
          for (uint i_orbital = 0; i_orbital < i_orbitals; i_orbital++) {
            uint j_orbital_finish = (i == j) ? i_orbital + 1 : j_orbitals;
            for (uint j_orbital = 0; j_orbital < j_orbital_finish;
                 j_orbital++) {
              num_dens_terms++;
              this->dens_counts[current_term_type]++;

              uint dens_ind =
                  (i + i_orbital) +
                  (2 * G2G::fortran_vars.m - ((j + j_orbital) + 1)) *
                      (j + j_orbital) / 2;
              this->dens_values.push_back(
                  G2G::fortran_vars.rmm_input_ndens1.data[dens_ind]);
              this->local2globaldens.push_back(dens_ind);
              local_dens_ind++;
            }
          }
        }
      }
    }
    // Pad the input arrays so the next term type has an aligned offset
    for (j = 0; j < QMMM_BLOCK_SIZE -
                        (term_type_counts[current_term_type] % QMMM_BLOCK_SIZE);
         j++) {
      this->func_code.push_back(
          func_code[term_type_offsets[current_term_type]]);  // Use the first
                                                             // code from this
                                                             // term type
      this->local_dens.push_back(
          local_dens[term_type_offsets[current_term_type]]);
    }
    for (j = 0; j < QMMM_BLOCK_SIZE -
                        (dens_counts[current_term_type] % QMMM_BLOCK_SIZE);
         j++) {
      this->dens_values.push_back(dens_values[dens_offsets[current_term_type]]);
      this->local2globaldens.push_back(
          local2globaldens[dens_offsets[current_term_type]]);
    }
  }

  // Pad the input so that out-of-range threads do a dummy calculation (same as
  // the first thread), rather than branching and idling
  /*  for (i = 0; i < QMMM_BLOCK_SIZE -
    (COALESCED_DIMENSION(term_type_counts[NUM_TERM_TYPES-1]) % QMMM_BLOCK_SIZE);
    i++) {
      this->func_code.push_back(func_code[term_type_offsets[NUM_TERM_TYPES-1]]);
      this->local_dens.push_back(local_dens[term_type_offsets[NUM_TERM_TYPES-1]]);
    }
    for (i = 0; i < QMMM_BLOCK_SIZE -
    (COALESCED_DIMENSION(dens_counts[NUM_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++)
    {
      this->dens_values.push_back(dens_values[dens_offsets[NUM_TERM_TYPES-1]]);
        this->local2globaldens.push_back(local2globaldens[dens_offsets[NUM_TERM_TYPES-1]]);
    }*/

  //  cout << "AINT NUMBER OF THREADS: " << num_terms << " (" <<
  //  this->func_code.size() << ")" << endl;
  //  cout << "AINT DENSITY TERMS: " << num_dens_terms << " (" <<
  //  this->dens_values.size() << ")" << endl;
}

#if AINT_MP && !FULL_DOUBLE
template class OSIntegral<float>;
#else
template class OSIntegral<double>;
#endif
}
