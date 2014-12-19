__device__ __constant__ uint TERM_TYPE_GAUSSIANS_EN[6] = { 1, 3, 9, 6, 18, 36 }; // How many individual partial fock terms each type (s-s,etc) is calculating
#define PP_SAME_FUNC_SIZE 6
#define DD_SAME_FUNC_SIZE 21

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
__global__ void gpu_qmmm_fock( uint num_terms, vec_type<scalar_type,2>* ac_values, uint* func2nuc, uint* func_code, uint* local_fock_ind,// uint fock_size,
                                 scalar_type* fock, uint global_stride, vec_type<scalar_type,3>* clatom_pos, scalar_type *clatom_chg )//, uint fock_out )
{

  assert(QMMM_FORCES_BLOCK_SIZE == 128);
  uint ffnum = index_x(blockDim, blockIdx, threadIdx);
  int tid = threadIdx.x;
  bool valid_thread = (ffnum < num_terms);

  // Each thread maps to a single pair of QM nuclei, so these forces are computed locally and accumulated at the end
  scalar_type prefactor;
  scalar_type my_fock[term_type==0? 1 : (term_type==1? 3 : (term_type==2? 9 : (term_type==3? 6 : (term_type==4? 18 : 36))))];
  for (uint i = 0; i < TERM_TYPE_GAUSSIANS_EN[term_type]; i++) {
    my_fock[i] = 0.0f;
  }
  uint fock_ind;
  bool same_func;

  {
    __shared__ vec_type<scalar_type,3> clatom_position_sh[QMMM_FORCES_BLOCK_SIZE];
    __shared__ scalar_type clatom_charge_sh[QMMM_FORCES_BLOCK_SIZE];

    scalar_type ai, aj, inv_two_zeta;
    scalar_type P[3], PmA[3], PmB[3];

    // TODO: each thread calculates its own zeta, overlap, etc here; should these be precalculated and saved (for use here and in Coulomb calculation)?
    {
      fock_ind = local_fock_ind[ffnum];
 
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

        //
        // Get the function values and nuclei for this thread
        //
        vec_type<scalar_type,2> ac1 = ac_values[f1 + cont1 * COALESCED_DIMENSION(gpu_m)];
        vec_type<scalar_type,2> ac2 = ac_values[f2 + cont2 * COALESCED_DIMENSION(gpu_m)];
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
  
      vec_type<scalar_type,3> A, B;
      A = gpu_atom_positions[nuc1];
      B = gpu_atom_positions[nuc2];
  
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

      vec_type<scalar_type,3> AmB = A - B;
      scalar_type ds2 = length2(AmB);
      scalar_type ksi = ((double)ai*(double)aj)/zeta;
      ovlap = exp(-ds2*ksi);

      prefactor = -2.0f * PI * cc * ovlap / zeta;
    }

    //
    // Outer loop: read in block of MM atom information into shared memory
    //
    for (int i = 0; i < gpu_clatoms; i += QMMM_FORCES_BLOCK_SIZE)
    {
      if (i + tid < gpu_clatoms) {
        clatom_position_sh[tid] = clatom_pos[i+tid];
        clatom_charge_sh[tid] = clatom_chg[i+tid];
      }
      __syncthreads();
      //
      // Inner loop: process block of MM atoms; each thread calculates a single primitive/primitive overlap force term
      //
      for (int j = 0; j < QMMM_FORCES_BLOCK_SIZE && i+j < gpu_clatoms; j++)
      {
        scalar_type PmC[3];
        {
          vec_type<scalar_type, 3> clatom_pos = clatom_position_sh[j];
          PmC[0] = P[0] - clatom_pos.x;
          PmC[1] = P[1] - clatom_pos.y;
          PmC[2] = P[2] - clatom_pos.z;
        }
        //
        // Do the core part of the forces calculation - the evaluation of the Obara-Saika recursion equations
        // This is where the different term types differ the most, so these are moved into separate files in the qmmm_terms directory
        // Current version: p-s through d-d are manually unrolled, and d-d is split up over six threads per primitive pair
        //
        // BEGIN TERM-TYPE DEPENDENT PART
        switch (term_type)
        {
          case 0:
          {
            #include "qmmm_terms/energy/ss.h"
            break;
          }
          case 1:
          {
            #include "qmmm_terms/energy/ps.h"
            break;
          }
          case 2:
          {
            #include "qmmm_terms/energy/pp.h"
            break;
          }
          case 3:
          {
            #include "qmmm_terms/energy/ds.h"
            break;
          }
          case 4:
          {
            #include "qmmm_terms/energy/dp_unrolled.h"
            break;
          }
          case 5:
          {
            #include "qmmm_terms/energy/dd_unrolled.h"
            break;
          }
        }
        // END TERM-TYPE DEPENDENT PART
      }
      __syncthreads();
    }
  }

  {

    uint bl_fock_ind = 0;//(fock_ind - min_ind) % TERM_TYPE_GAUSSIANS_EN[term_type];

    __shared__ uint fock_ind_sh[QMMM_FORCES_BLOCK_SIZE];
    //__shared__ bool fock_flags[QMMM_FORCES_BLOCK_SIZE];
    __shared__ bool same_func_sh[QMMM_FORCES_BLOCK_SIZE];
    __shared__ scalar_type fock_sh[QMMM_FORCES_BLOCK_SIZE];

    fock_ind_sh[tid] = fock_ind;
    __syncthreads();
    uint curr_ind = fock_ind_sh[0], curr_bl_ind = 0;
    for (int i = 1; i < QMMM_FORCES_BLOCK_SIZE; i++) {
      curr_bl_ind += curr_ind != fock_ind_sh[i];
      curr_ind = fock_ind_sh[i];
      bl_fock_ind = (curr_ind == fock_ind) * curr_bl_ind + (curr_ind != fock_ind) * bl_fock_ind;
    }

    //
    // Reduce the partial fock elements
    //
    //
    // First figure out which fock elements this block contains
    //
    //fock_flags[tid] = false;
    same_func_sh[tid] = false;
    __syncthreads();
    //fock_flags[bl_fock_ind] = true;
    same_func_sh[bl_fock_ind] = same_func;
    fock_ind_sh[bl_fock_ind] = fock_ind;
    __syncthreads();
    for (int i = 0; i <= curr_bl_ind/*QMMM_FORCES_BLOCK_SIZE*/; i++)
    {
      // Only for this block's fock
      //if (fock_flags[i] == true)
      //{
        //
        // Load the individual thread's fock terms into the appropriate shared location
        //
        bool use_fock = bl_fock_ind == i;
        uint last_j = TERM_TYPE_GAUSSIANS_EN[term_type];
        if (term_type == 2 && same_func_sh[i]) {
          last_j = PP_SAME_FUNC_SIZE;
        } else if (term_type == 5 && same_func_sh[i]) {
          last_j = DD_SAME_FUNC_SIZE;
        }
        for (int j = 0; j < last_j/*TERM_TYPE_GAUSSIANS_EN[term_type]*/; j++) {
          fock_sh[tid] = valid_thread * use_fock * prefactor * my_fock[j];
          __syncthreads();

          //
          // Reduce the fock terms
          //
          if (tid < QMMM_FORCES_HALF_BLOCK)
          {
            fock_sh[tid] += fock_sh[tid+QMMM_FORCES_HALF_BLOCK];
          }
          __syncthreads();

          if (tid < WARP_SIZE) { warpReduce<scalar_type>(fock_sh, tid); }
          if (tid == 0) {
            //printf("%d %d %d %d %d\n",fock_ind_sh[i],j,blockIdx.x,fock_out,global_stride*(fock_ind_sh[i]+j)+(blockIdx.x-fock_out));
            fock[fock_ind_sh[i]+j+global_stride*blockIdx.x] = fock_sh[0];
          }
          __syncthreads();
        }
      //}
    }
  }
}

template<class scalar_type>
__global__ void zero_fock( scalar_type* fock, uint global_stride, uint fock_length )
{
  uint my_x = index_x(blockDim, blockIdx, threadIdx);
  uint my_y = blockDim.y * blockIdx.y + threadIdx.y;

  if (my_x < global_stride && my_y < fock_length) {
    fock[my_x + global_stride * my_y] = 0.0f;
  }
}

template<class scalar_type>
__global__ void gpu_qmmm_fock_reduce( scalar_type* fock, scalar_type* dens, scalar_type* energies, uint stride, uint depth, uint width )
{
  uint my_fock_ind = index_x(blockDim, blockIdx, threadIdx);
  uint tid = threadIdx.x;
  scalar_type my_partial_fock = 0.0f;

  if (my_fock_ind < width) {
    for (uint i = 0; i < depth; i++) {
      my_partial_fock += fock[stride*i + my_fock_ind];
    }
    fock[my_fock_ind] = my_partial_fock;
  }

  __shared__ scalar_type energies_sh[QMMM_REDUCE_BLOCK_SIZE];
  scalar_type my_dens = 0.0f;
  if (my_fock_ind < width) {
    my_dens = dens[my_fock_ind];
  }
  energies_sh[tid] = my_partial_fock * my_dens;
  __syncthreads();

  if (tid < 64)
  {
    energies_sh[tid] += energies_sh[tid+64];
  }
  __syncthreads();

  if (tid < WARP_SIZE) { warpReduce<scalar_type>(energies_sh, tid); }
  if (tid == 0) { energies[blockIdx.x] = energies_sh[0]; }
}





