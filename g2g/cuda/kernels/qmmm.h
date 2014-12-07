//#define SQRT_PI 1.772453851
#define PI 3.141592653589793238462643383f

template<class scalar_type>
__device__ scalar_type lio_gamma(uint m, scalar_type U)
{
  scalar_type funct1,funct2;
  int it;
  scalar_type ti,delt,delt2,delt3,delt4,delt5;

  int s = (U<=43.975);
  // Calculate small-U branch value of F(m,U)
  //if (U <= 43.975) {
  it = 20.0 * (U + 0.025);
  it = s * it;
  ti = it;
  delt = U - 0.05 * ti;
  delt3 = delt * 0.333333333333333;
  delt4 = 0.25 * delt;
  delt2 = delt4 + delt4;
  delt5 = 0.20 * delt;

  scalar_type tf0,tf1,tf2,tf3,tf4,tf5;
  tf0 = gpu_str[it+m*880];
  tf1 = gpu_str[it+(m+1)*880];
  tf2 = gpu_str[it+(m+2)*880];
  tf3 = gpu_str[it+(m+3)*880];
  tf4 = gpu_str[it+(m+4)*880];
  tf5 = gpu_str[it+(m+5)*880];

  funct1 = tf0-delt * (tf1-delt2 * (tf2-delt3 * (tf3-delt4 * (tf4-delt5 * tf5))));
  //} else {

  // Calculate large-U branch value of F(m,U)
  funct2 = gpu_fac[m]/(powf(U,m)*sqrtf(U));
  //}

  // Return the appropriate branch while avoiding warp divergence
  // TODO: need to verify calculating both branches is better for performance; it works better on small systems but need to test on big test cases
  return s*funct1+(1-s)*funct2;
  //return funct1;
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

// TODO: currently, one thread maps to one primitive-primitive overlap force term; is there a better mapping? (thread to function, thread to sub-shell, etc)
// Also, should the loop over MM atoms be broken up to be done by multiple blocks rather than a block looping over every MM atom?
template<class scalar_type>
__global__ void gpu_qmmm_forces( uint num_terms, scalar_type* a_values1, scalar_type* a_values2, scalar_type* cc_values, scalar_type* dens_values,// uint* func1, uint* func2,
                                 uint* nuclei1, uint* nuclei2, vec_type<scalar_type,3>* mm_forces, vec_type<scalar_type,3>* qm_forces )//, uint s_func_end, uint p_func_end )
{

  uint ffnum = index_x(blockDim, blockIdx, threadIdx);
  int tid = threadIdx.x;
  bool valid_thread = (ffnum < num_terms);

  // Each thread maps to a single pair of QM nuclei, so these forces are computed locally and accumulated at the end
  scalar_type A_force[3], B_force[3];
  uint nuc1, nuc2;
  scalar_type prefactor_qm;

  {
    __shared__ vec_type<scalar_type,3> clatom_position_sh[QMMM_FORCES_BLOCK_SIZE];
    __shared__ scalar_type clatom_charge_sh[QMMM_FORCES_BLOCK_SIZE];
    // Shared memory space for reduction of MM atom force terms
    __shared__ scalar_type C_force[3][QMMM_FORCES_BLOCK_SIZE];

    scalar_type zeta;
    scalar_type ai, aj, prefactor_mm;
    scalar_type P[3], PmA[3], PmB[3];
    //int max_m = 1;
    // TODO: each thread calculates its own zeta, overlap, etc here; should these be precalculated and saved (for use here and in Coulomb calculation)?
    {
      scalar_type ovlap, cc, dens;
      ai = a_values1[ffnum]; aj = a_values2[ffnum];
      cc = cc_values[ffnum];
      dens = dens_values[ffnum];
  
      vec_type<scalar_type,3> A, B;
      nuc1 = nuclei1[ffnum]; nuc2 = nuclei2[ffnum];
      A = gpu_atom_positions[nuc1];
      B = gpu_atom_positions[nuc2];
  
      zeta = ai + aj;
      P[0] = (A.x*ai + B.x*aj) / zeta;
      P[1] = (A.y*ai + B.y*aj) / zeta;
      P[2] = (A.z*ai + B.z*aj) / zeta;
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

      prefactor_mm = -dens * cc * 4.0f * PI * ovlap;
      prefactor_qm = prefactor_mm / (2.0f * zeta);
    
      //uint f1 = func1[ffnum], f2 = func2[ffnum];
      //max_m += (f1 >= s_func_end) + (f1 >= p_func_end);
      //max_m += (f2 >= s_func_end) + (f2 >= p_func_end);
    }

    // Outer loop: read in block of MM atom information into shared memory
    for (int i = 0; i < gpu_clatoms; i += QMMM_FORCES_BLOCK_SIZE)
    {
      if (i + tid < gpu_clatoms) {
        clatom_position_sh[tid] = gpu_clatom_positions[i+tid];
        clatom_charge_sh[tid] = gpu_clatom_charges[i+tid];
      }
      __syncthreads();
      // Inner loop: process block of MM atoms; each thread calculates a single primitive/primitive overlap force term
      for (int j = 0; j < QMMM_FORCES_BLOCK_SIZE && i+j < gpu_clatoms; j++)
      {
        scalar_type PmC[3], F_mU[2]; // ONLY FOR S-S
        {
          vec_type<scalar_type, 3> clatom_pos = clatom_position_sh[j];
          PmC[0] = P[0] - clatom_pos.x;
          PmC[1] = P[1] - clatom_pos.y;
          PmC[2] = P[2] - clatom_pos.z;
        }
        {
          scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * zeta;
          //F_mU[0] = (SQRT_PI / (2*sqrtU)) * erff(sqrtU);
          for (int m = 0; m <= 1; m++) // ONLY FOR S-S
          {
            // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
            F_mU[m] = lio_gamma<scalar_type>(m,U);
            //F_mU[m] = fetch(qmmm_F_values_tex,(float)(U/gamma_inc-0.5f),(float)(m+0.5f));
          }
        }

        // BEGIN calculation of individual (single primitive-primitive overlap) force terms
        // TODO: This block can probably be put into separate headers for s-s, p-s, etc then included in a single main kernel file
        {
          scalar_type A_force_term, B_force_term, C_force_term;
          scalar_type mm_charge = clatom_charge_sh[j];
  
          for (int grad_l = 0; grad_l < 3; grad_l++)
          {
            C_force_term = PmC[grad_l] * F_mU[1];
            A_force_term = PmA[grad_l] * F_mU[0] - C_force_term;
            B_force_term = PmB[grad_l] * F_mU[0] - C_force_term;
  
            A_force[grad_l] += 2.0f * ai * mm_charge * A_force_term;
            B_force[grad_l] += 2.0f * aj * mm_charge * B_force_term;
            // Out-of-range threads contribute 0 to the force
            C_force[grad_l][tid] = valid_thread * prefactor_mm * mm_charge * C_force_term;
          }
        }
        // END individual force terms
        __syncthreads();

        // TODO: should we do the per-block reduction here in this loop? or should each thread save its value to global memory for later accumulation?
        // BEGIN reduction of MM atom force terms
        // IMPORTANT: ASSUMING WARP SIZE OF 32 (or maybe assuming memory access granularity = warp size?...assuming something here anyway)
        // First half of block does x,y
        if (tid < 64)
        {
          C_force[0][tid] += C_force[0][tid+64];
          C_force[1][tid] += C_force[1][tid+64];
          //C_force[2][tid] += C_force[2][tid+64];
        }
        // Second half does z (probably doesn't make much of a difference)
        else
        {
          C_force[2][tid-64] += C_force[2][tid];
        }
        __syncthreads();
        // first warp does x
        if (tid < 32)      { warpReduce<scalar_type>(C_force[0], tid); }
        // second warp does y
        else if (tid < 64) { warpReduce<scalar_type>(C_force[1], tid-32); }
        // third warp does z
        else if (tid < 96) { warpReduce<scalar_type>(C_force[2], tid-64); }

        {
          uint global_stride = COALESCED_DIMENSION(gridDim.x);
          // TODO: tried turning this into one global read to get the force vector object, but didn't seem to improve performance, maybe there's a better way?
          if (tid == 0)       { mm_forces[global_stride*(i+j)+blockIdx.x].x = C_force[0][0]; }
          else if (tid == 32) { mm_forces[global_stride*(i+j)+blockIdx.x].y = C_force[1][0]; }
          else if (tid == 64) { mm_forces[global_stride*(i+j)+blockIdx.x].z = C_force[2][0]; }
        }
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
        if (tid < 64)
        {
          QM_force[0][tid] += QM_force[0][tid+64];
          QM_force[1][tid] += QM_force[1][tid+64];
        }
        // Second half does z
        else
        {
          QM_force[2][tid-64] += QM_force[2][tid];
        }
        __syncthreads();
        // first warp does x
        if (tid < 32)      { warpReduce<scalar_type>(QM_force[0], tid); }
        // second warp does y
        else if (tid < 64) { warpReduce<scalar_type>(QM_force[1], tid-32); }
        // third warp does z
        else if (tid < 96) { warpReduce<scalar_type>(QM_force[2], tid-64); }
        {
          uint global_stride = COALESCED_DIMENSION(gridDim.x);
          if (tid == 0)       { qm_forces[global_stride*i+blockIdx.x].x = QM_force[0][0]; }
          else if (tid == 32) { qm_forces[global_stride*i+blockIdx.x].y = QM_force[1][0]; }
          else if (tid == 64) { qm_forces[global_stride*i+blockIdx.x].z = QM_force[2][0]; }
        }
        __syncthreads();
      }
      // At this point, the global QM array is uninitialized; since we'll accumulate all entries, it needs to be zeroed
      // TODO: Zeroing out the partial qm array before this kernel might be better
      else
      {
        uint global_stride = COALESCED_DIMENSION(gridDim.x);
        if (tid == 0)       { qm_forces[global_stride*i+blockIdx.x].x = 0.0f; }
        else if (tid == 32) { qm_forces[global_stride*i+blockIdx.x].y = 0.0f; }
        else if (tid == 64) { qm_forces[global_stride*i+blockIdx.x].z = 0.0f; }
      }
    }
  }
}

