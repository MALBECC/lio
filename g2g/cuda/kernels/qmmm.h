#define WARP_SIZE 32
#define WARP_SIZE2 2*WARP_SIZE
#define WARP_SIZE3 3*WARP_SIZE
#define QMMM_FORCES_HALF_BLOCK QMMM_FORCES_BLOCK_SIZE/2

#define PI 3.141592653589793238462643383f

template<class scalar_type,int m_max>
__device__ void lio_gamma(scalar_type* __restrict__ F_mU, scalar_type U)
{
  int it;
  scalar_type ti,delt,delt2,delt3,delt4,delt5;

  // Calculate small-U branch value of F(m,U)
  // TODO: need to rethink how this branch (Taylor series expansion) is calculated
  // There's 6 reads to gpu_str, and currently U is not ordered wrt thread order, so the access pattern is terrible
  // Ideas: -reorder threads wrt U
  //        -place gpu_str in texture memory
  if (U <= 43.975)
  {
    it = 20.0 * (U + 0.025);
    //it = s * it;
    ti = it;
    delt = U - 0.05 * ti;
    delt3 = delt * 0.333333333333333;
    delt4 = 0.25 * delt;
    delt2 = delt4 + delt4;
    delt5 = 0.20 * delt;

    scalar_type tf0,tf1,tf2,tf3,tf4,tf5;
    tf0 = fetch(qmmm_str_tex,(float)it,0.0f);//qmmm_str[it];
    tf1 = fetch(qmmm_str_tex,(float)it,1.0f);//qmmm_str[it+880];
    tf2 = fetch(qmmm_str_tex,(float)it,2.0f);//qmmm_str[it+1760];
    tf3 = fetch(qmmm_str_tex,(float)it,3.0f);//qmmm_str[it+2640];
    tf4 = fetch(qmmm_str_tex,(float)it,4.0f);//qmmm_str[it+3520];
    tf5 = fetch(qmmm_str_tex,(float)it,5.0f);//qmmm_str[it+4400];
    F_mU[0] = tf0-delt * (tf1-delt2 * (tf2-delt3 * (tf3-delt4 * (tf4-delt5 * tf5))));
    for (uint m = 1; m <= m_max; m++) {
      tf0 = tf1;
      tf1 = tf2;
      tf2 = tf3;
      tf3 = tf4;
      tf4 = tf5;
      tf5 = fetch(qmmm_str_tex,(float)it,(float)(m+5.0f));//qmmm_str[it+(m+5)*880];

      F_mU[m] = tf0-delt * (tf1-delt2 * (tf2-delt3 * (tf3-delt4 * (tf4-delt5 * tf5))));
    }
  }
  // Calculate large-U branch value of F(m,U)
  else
  {
    scalar_type sqrtU = sqrtf(U);
    for (uint m = 0; m <= m_max; m++) {
      F_mU[m] = gpu_fac[m]/(powf(U,m)*sqrtU);//sqrtf(U));
    }
  }
}

/*template<class scalar_type>
__global__ void precompute_gamma( uint gamma_length, scalar_type gamma_inc, scalar_type* gpu_gamma )
{
  uint x = index_x(blockDim, blockIdx, threadIdx);
  uint m = threadIdx.y;
  scalar_type U = (x+1) * gamma_inc;

  if (x < gamma_length) {
    gpu_gamma[x + gamma_length*m] = lio_gamma<scalar_type>(m, U);
  }
}

template<class scalar_type>
__global__ void gpu_test_fmu_tex( scalar_type gamma_inc_test, scalar_type gamma_inc )
{

  int x = index_x(blockDim,blockIdx,threadIdx);
  int m = threadIdx.y;

  scalar_type U = (x+1)*gamma_inc_test;
  printf("F(%d,%f) = %.6e (lio) = %.6e (tex) \n", m,U,lio_gamma(m,U),fetch(qmmm_F_values_tex,(float)(U/gamma_inc-0.5f),(float)(m+0.5f)));

}*/

// Modified from presentation "Optimizing Parallel Reduction in CUDA" by Mark Harris (Nvidia)
template<class scalar_type>
__device__ void warpReduce(volatile scalar_type *sdata, unsigned int tid)
{
  sdata[tid] += sdata[tid + 32];
  sdata[tid] += sdata[tid + 16];
  sdata[tid] += sdata[tid + 8];
  sdata[tid] += sdata[tid + 4];
  sdata[tid] += sdata[tid + 2];
  sdata[tid] += sdata[tid + 1];
}

__device__ __constant__ uint TERM_TYPE_GAUSSIANS[6] = { 1, 3, 9, 6, 18, 6 };

// TODO: currently, one thread maps to one primitive-primitive overlap force term; is there a better mapping? (thread to function, thread to sub-shell, etc)
// Also, should the loop over MM atoms be broken up to be done by multiple blocks rather than a block looping over every MM atom?
template<class scalar_type, uint term_type>
__global__ void gpu_qmmm_forces( uint num_terms, vec_type<scalar_type,2>* ac_values, uint* func2nuc,/*scalar_type* a_values1, scalar_type* a_values2, scalar_type* cc_values,*/ scalar_type* dens_values, uint* func_code, uint* local_dens,//uint* orbital1, uint* orbital2, //uint* func1, uint* func2,
                                 /*uint* nuclei1, uint* nuclei2,*/ vec_type<scalar_type,3>* mm_forces, vec_type<scalar_type,3>* qm_forces, uint global_stride, vec_type<scalar_type,3>* clatom_pos, scalar_type *clatom_chg )//, uint s_func_end, uint p_func_end )
{

  assert(QMMM_FORCES_BLOCK_SIZE == 128);
  uint ffnum = index_x(blockDim, blockIdx, threadIdx);
  int tid = threadIdx.x;
  bool valid_thread = (ffnum < num_terms);

  // Each thread maps to a single pair of QM nuclei, so these forces are computed locally and accumulated at the end
  scalar_type A_force[3] = { 0.0f,0.0f,0.0f }, B_force[3] = { 0.0f,0.0f,0.0f };
  uint nuc1, nuc2;
  //const uint num_gauss = TERM_TYPE_GAUSSIANS[term_type];
  scalar_type prefactor_qm;

  {
    __shared__ vec_type<scalar_type,3> clatom_position_sh[QMMM_FORCES_BLOCK_SIZE];
    __shared__ scalar_type clatom_charge_sh[QMMM_FORCES_BLOCK_SIZE];
    // Shared memory space for reduction of MM atom force terms
    __shared__ scalar_type C_force[3][QMMM_FORCES_BLOCK_SIZE];

    scalar_type ai, aj, prefactor_mm, inv_two_zeta;
    scalar_type dens[term_type==0? 1 : (term_type==1? 3 : (term_type==2? 9 : (term_type==3? 6 : (term_type==4? 18 : 6))))];
    scalar_type P[3], PmA[3], PmB[3];
    bool same_func = false;
    //uint orb1 = 0, orb2 = 0;

    uint d1_l1, d1_l2; // These are only used by the d-d kernel

    //int max_m = 1;
    // TODO: each thread calculates its own zeta, overlap, etc here; should these be precalculated and saved (for use here and in Coulomb calculation)?
    {
      scalar_type cc;
      {
        uint my_func_code = func_code[ffnum];

        uint div = MAX_CONTRACTIONS;
        uint cont2 = my_func_code % div;
        my_func_code /= div;
        uint cont1 = my_func_code % div;
        my_func_code /= div;

        div = gpu_m;
        uint f2 = my_func_code % div;
        my_func_code /= div;
        uint f1 = my_func_code;// % div;
        //my_func_code /= div;

        uint orb1;
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
        }

        uint dens_ind = local_dens[ffnum];
        if (term_type == 2 && same_func) {
          for (uint i = 0; i < 6; i++) {
            dens[i] = dens_values[dens_ind+i];
          }
        } else if (term_type == 5 && same_func) {
          for (uint i = 0; i <= orb1; i++) {
            dens[i] = dens_values[dens_ind+i];
          }
        } else {
          for (uint i = 0; i < TERM_TYPE_GAUSSIANS[term_type]; i++) {
            dens[i] = dens_values[dens_ind+i];
          }
        }

        // find1, find2 needed if nuc / function value arrays contracted to total # of functions rather than m (p/d expanded out)
        // The extra math here seems to slow down the kernel more than the memory savings, but behavior could be different in bigger systems
        //uint find1 = f1, find2 = f2;
        /*if (term_type == 1 || term_type == 2) {
          orb1  = (f1 - gpu_s_funcs) % 3;
          //find1 = (f1 - gpu_s_funcs) / 3 + gpu_s_funcs;
        } else if (term_type >= 3) {
          orb1  = (f1 - (gpu_s_funcs+gpu_p_funcs*3)) % 6;
          //find1 = (f1 - (gpu_s_funcs+gpu_p_funcs*3)) / 6 + gpu_s_funcs+gpu_p_funcs*3;
        }
        if (term_type == 2 || term_type == 4) {
          orb2  = (f2 - gpu_s_funcs) % 3;
          //find2 = (f2 - gpu_s_funcs) / 3 + gpu_s_funcs;
        } else if (term_type == 5) {
          orb2  = (f2 - (gpu_s_funcs+gpu_p_funcs*3)) % 6;
          //find2 = (f2 - (gpu_s_funcs+gpu_p_funcs*3)) / 6 + gpu_s_funcs+gpu_p_funcs*3;
        }*/

        //printf("%d %d %d\n",f1,find1,term_type);
        vec_type<scalar_type,2> ac1 = ac_values[f1 + cont1 * COALESCED_DIMENSION(gpu_m)];//total_funcs)];
        vec_type<scalar_type,2> ac2 = ac_values[f2 + cont2 * COALESCED_DIMENSION(gpu_m)];//total_funcs)];
        ai = ac1.x;
        aj = ac2.x;
        //printf("%d: %d %d %d %d %.4f %.4f (%d)\n",term_type,f1,cont1,f2,cont2,ai,aj,valid_thread);
        cc = ac1.y * ac2.y;

        nuc1 = func2nuc[f1];
        nuc2 = func2nuc[f2];
      }

      scalar_type ovlap;
      //ai = a_values1[ffnum]; aj = a_values2[ffnum];
      //cc = cc_values[ffnum];
      //dens = dens_values[ffnum];
  
      vec_type<scalar_type,3> A, B;
      //nuc1 = nuclei1[ffnum]; nuc2 = nuclei2[ffnum];
      A = gpu_atom_positions[nuc1];
      B = gpu_atom_positions[nuc2];
  
      scalar_type zeta = ai + aj;
      //if (valid_thread && term_type==0) printf("THREADNUM: %d %f\n",ffnum,zeta);
      //scalar_type inv_zeta = 1.0f / zeta;
      inv_two_zeta = 1.0f / (2.0f * zeta);
      // Noticed that precomputing the inverse and multiplying rather than dividing here led to a 3x worse round-off error in the QM forces (compared to intsolG)
      // However, using inv_two_zeta for prefactor_qm further down doesn't have an effect...maybe because P goes into PmA/PmB/PmC/U which are then used in more
      // operations than prefactor_qm, which is just multiplied once at the end?
      P[0] = (A.x*ai + B.x*aj) / zeta;//* inv_zeta;
      P[1] = (A.y*ai + B.y*aj) / zeta;//* inv_zeta;
      P[2] = (A.z*ai + B.z*aj) / zeta;//* inv_zeta;
      PmA[0] = P[0] - A.x;
      PmA[1] = P[1] - A.y;
      PmA[2] = P[2] - A.z;
      PmB[0] = P[0] - B.x;
      PmB[1] = P[1] - B.y;
      PmB[2] = P[2] - B.z;

      vec_type<scalar_type,3> AmB = A - B;
      scalar_type ds2 = length2(AmB);
      scalar_type ksi = ai*aj/zeta;
      ovlap = exp(-ds2*ksi);

      if (term_type == 0) {
        prefactor_mm = -dens[0] * cc * 4.0f * PI * ovlap;
      } else {
        prefactor_mm = -cc * 4.0f * PI * ovlap;
      }
      prefactor_qm = prefactor_mm * inv_two_zeta;
    }

    // Outer loop: read in block of MM atom information into shared memory
    for (int i = 0; i < gpu_clatoms; i += QMMM_FORCES_BLOCK_SIZE)
    {
      if (i + tid < gpu_clatoms) {
        clatom_position_sh[tid] = clatom_pos[i+tid];
        clatom_charge_sh[tid] = clatom_chg[i+tid];
      }
      __syncthreads();
      // Inner loop: process block of MM atoms; each thread calculates a single primitive/primitive overlap force term
      for (int j = 0; j < QMMM_FORCES_BLOCK_SIZE && i+j < gpu_clatoms; j++)
      {
        {
          scalar_type PmC[3];
          {
            vec_type<scalar_type, 3> clatom_pos = clatom_position_sh[j];
            PmC[0] = P[0] - clatom_pos.x;
            PmC[1] = P[1] - clatom_pos.y;
            PmC[2] = P[2] - clatom_pos.z;
          }
          // BEGIN TERM-TYPE DEPENDENT PART
          switch (term_type)
          {
            case 0:
            {
              #include "qmmm_terms/ss.h"
              break;
            }
            case 1:
            {
              #include "qmmm_terms/ps_unrolled.h"
              break;
            }
            case 2:
            {
              #include "qmmm_terms/pp_unrolled.h"
              break;
            }
            case 3:
            {
              #include "qmmm_terms/ds_unrolled.h"
              break;
            }
            case 4:
            {
              #include "qmmm_terms/dp_unrolled.h"
              break;
            }
            case 5:
            {
              #include "qmmm_terms/dd_split_unrolled.h"//_unrolled.h"
              break;
            }
          }
          // END TERM-TYPE DEPENDENT PART
        }

        __syncthreads();

        // TODO: should we do the per-block reduction here in this loop? or should each thread save its value to global memory for later accumulation?
        // BEGIN reduction of MM atom force terms
        // IMPORTANT: ASSUMING WARP SIZE OF 32 (or maybe assuming memory access granularity = warp size?...assuming something here anyway)
        // ALSO ASSUMING BLOCK SIZE OF 128
        // First half of block does x,y
        if (tid < QMMM_FORCES_HALF_BLOCK)
        {
          C_force[0][tid] += C_force[0][tid+QMMM_FORCES_HALF_BLOCK];
          C_force[1][tid] += C_force[1][tid+QMMM_FORCES_HALF_BLOCK];
          //C_force[2][tid] += C_force[2][tid+64];
        }
        // Second half does z (probably doesn't make much of a difference)
        else
        {
          C_force[2][tid-QMMM_FORCES_HALF_BLOCK] += C_force[2][tid];
        }
        __syncthreads();
        // first warp does x
        if (tid < WARP_SIZE)       { warpReduce<scalar_type>(C_force[0], tid); }
        // second warp does y
        else if (tid < WARP_SIZE2) { warpReduce<scalar_type>(C_force[1], tid-WARP_SIZE); }
        // third warp does z
        else if (tid < WARP_SIZE3) { warpReduce<scalar_type>(C_force[2], tid-WARP_SIZE2); }

        //uint global_stride = COALESCED_DIMENSION(gridDim.x);
        // TODO: tried turning this into one global read to get the force vector object, but didn't seem to improve performance, maybe there's a better way?
        if (tid == 0)               { mm_forces[global_stride*(i+j)+blockIdx.x].x = C_force[0][0]; }
        else if (tid == WARP_SIZE)  { mm_forces[global_stride*(i+j)+blockIdx.x].y = C_force[1][0]; }
        else if (tid == WARP_SIZE2) { mm_forces[global_stride*(i+j)+blockIdx.x].z = C_force[2][0]; }
        // END reduction

        __syncthreads();
      }
    }
  }

  // TODO: (same question as for the MM forces) - should we do the per-block reduction in this kernel?
  // Reduce the QM force terms (not sure this is the best way to do things yet)
  {
    __shared__ bool nuc_flags[MAX_ATOMS];
    __shared__ scalar_type QM_force[3][QMMM_FORCES_BLOCK_SIZE];

    // First figure out which nuclei are present in this block
    for (int i = 0; i < gpu_atoms; i += QMMM_FORCES_BLOCK_SIZE) {
      if (i+tid<gpu_atoms) nuc_flags[i+tid] = false;
    }
    __syncthreads();
    nuc_flags[nuc1] = true;
    nuc_flags[nuc2] = true;
    __syncthreads();
    for (int i = 0; i < gpu_atoms; i++)
    {
      // Only for this block's nuclei
      if (nuc_flags[i] == true)
      {
        // Load the individual thread's force terms into the appropriate shared location
        bool useA = nuc1 == i, useB = nuc2 == i; 
        QM_force[0][tid] = valid_thread * prefactor_qm * (useA * A_force[0] + useB * B_force[0]);
        QM_force[1][tid] = valid_thread * prefactor_qm * (useA * A_force[1] + useB * B_force[1]);
        QM_force[2][tid] = valid_thread * prefactor_qm * (useA * A_force[2] + useB * B_force[2]);
        __syncthreads();

        // Reduce the force terms
        // First half of block does x,y
        if (tid < QMMM_FORCES_HALF_BLOCK)
        {
          QM_force[0][tid] += QM_force[0][tid+QMMM_FORCES_HALF_BLOCK];
          QM_force[1][tid] += QM_force[1][tid+QMMM_FORCES_HALF_BLOCK];
        }
        // Second half does z
        else
        {
          QM_force[2][tid-QMMM_FORCES_HALF_BLOCK] += QM_force[2][tid];
        }
        __syncthreads();
        // first warp does x
        if (tid < WARP_SIZE)       { warpReduce<scalar_type>(QM_force[0], tid); }
        // second warp does y
        else if (tid < WARP_SIZE2) { warpReduce<scalar_type>(QM_force[1], tid-WARP_SIZE); }
        // third warp does z
        else if (tid < WARP_SIZE3) { warpReduce<scalar_type>(QM_force[2], tid-WARP_SIZE2); }

        //uint global_stride = COALESCED_DIMENSION(gridDim.x);
        if (tid == 0)               { qm_forces[global_stride*i+blockIdx.x].x = QM_force[0][0]; }
        else if (tid == WARP_SIZE)  { qm_forces[global_stride*i+blockIdx.x].y = QM_force[1][0]; }
        else if (tid == WARP_SIZE2) { qm_forces[global_stride*i+blockIdx.x].z = QM_force[2][0]; }
        __syncthreads();
      }
      // At this point, the global QM array is uninitialized; since we'll accumulate all entries, it needs to be zeroed
      // TODO: Zeroing out the partial qm array before this kernel might be better
      else
      {
        //uint global_stride = COALESCED_DIMENSION(gridDim.x);
        if (tid == 0)       { qm_forces[global_stride*i+blockIdx.x].x = 0.0f; }
        else if (tid == 32) { qm_forces[global_stride*i+blockIdx.x].y = 0.0f; }
        else if (tid == 64) { qm_forces[global_stride*i+blockIdx.x].z = 0.0f; }
      }
    }
  }
}

