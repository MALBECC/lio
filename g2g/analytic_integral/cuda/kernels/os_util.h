#define WARP_SIZE 32
#define WARP_SIZE2 2*WARP_SIZE
#define WARP_SIZE3 3*WARP_SIZE

#define QMMM_FORCES_HALF_BLOCK QMMM_BLOCK_SIZE/2

#define PP_SAME_FUNC_SIZE 6
#define DD_SAME_FUNC_SIZE 21

#define PI 3.14159265358979312
#define PI52 34.9868366552497108
//#define PI 3.141592653589793238462643383
//#define PI52 17.49341832762486284626282167987

#if FULL_DOUBLE || !AINT_MP
static __inline__ __device__ double os_fetch_double(texture<int2, 2> t, float x, float y)
{
    int2 v = tex2D(t,x,y);
    return __hiloint2double(v.y, v.x);
}
#define os_fetch(t,x,y) os_fetch_double(t,x,y)
#else
#define os_fetch(t,x,y) tex2D(t,x,y)
#endif

//
// Calculates F(m,U) values for m = 0 to max_m (F(m,U) is used in the Obara-Saika recursion relations)
// F(m,U) is calculated by one of two methods based on the value of U
// For small U (<=43.975), a truncated Taylor series approximation is used, using the precalculated table STR (qmmm_str_tex)
// For large U (>43.975) , the asymptotic value of the incomplete gamma (basically a regular gamma) is used, with prefactors precalculated in FAC (gpu_fac)
//
template<class scalar_type,int m_max>
__device__ void lio_gamma(scalar_type* __restrict__ F_mU, scalar_type U)
{
  int it;
  scalar_type ti,delt,delt2,delt3,delt4,delt5;

  // Calculate small-U branch value of F(m,U)
  // TODO: need to rethink how this branch (Taylor series expansion) is calculated
  // There's 6 reads to gpu_str, and currently U is not ordered wrt thread order, so the access pattern is terrible
  // Ideas: -reorder threads wrt U (seems to make things worse)
  //        -place gpu_str in texture memory (current implementation, doesn't seem to improve over global memory)
  if (U <= 43.975)
  {
    it = 20.0 * (U + 0.025);
    ti = it;
    delt = U - 0.05 * ti;
    delt3 = delt * 0.333333333333333;
    delt4 = 0.25 * delt;
    delt2 = delt4 + delt4;
    delt5 = 0.20 * delt;

    scalar_type tf0,tf1,tf2,tf3,tf4,tf5;
    tf0 = os_fetch(str_tex,(float)it,0.0);//qmmm_str[it];
    tf1 = os_fetch(str_tex,(float)it,1.0);//qmmm_str[it+880];
    tf2 = os_fetch(str_tex,(float)it,2.0);//qmmm_str[it+1760];
    tf3 = os_fetch(str_tex,(float)it,3.0);//qmmm_str[it+2640];
    tf4 = os_fetch(str_tex,(float)it,4.0);//qmmm_str[it+3520];
    tf5 = os_fetch(str_tex,(float)it,5.0);//qmmm_str[it+4400];

    F_mU[0] = tf0-delt * (tf1-delt2 * (tf2-delt3 * (tf3-delt4 * (tf4-delt5 * tf5))));
    for (uint m = 1; m <= m_max; m++) {
      tf0 = tf1;
      tf1 = tf2;
      tf2 = tf3;
      tf3 = tf4;
      tf4 = tf5;
      tf5 = os_fetch(str_tex,(float)it,(float)(m+5.0));//qmmm_str[it+(m+5)*880];

      F_mU[m] = tf0-delt * (tf1-delt2 * (tf2-delt3 * (tf3-delt4 * (tf4-delt5 * tf5))));
    }
  }
  // Calculate large-U branch value of F(m,U)
  else
  {
    scalar_type sqrtU = sqrt(U);
    scalar_type powU = 1.0;
    for (uint m = 0; m <= m_max; m++) {
      F_mU[m] = gpu_fac[m]/(powU*sqrtU);//sqrtf(U));
      powU *= U;
    }
  }
}

// 
// Reduce an array within a single warp - no need for synchronization between steps (but the "volatile" keyword is needed to avoid register-caching data race issues)
//
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

//
// Zero the partial Fock array
//
template<class scalar_type>
__global__ void zero_fock( scalar_type* fock, uint global_stride, uint fock_length )
{
  uint my_x = index_x(blockDim, blockIdx, threadIdx);
  uint my_y = blockDim.y * blockIdx.y + threadIdx.y;

  if (my_x < global_stride && my_y < fock_length) {
    fock[my_x + global_stride * my_y] = 0.0;
  }
}
//
// Zero the partial forces array
//
template<class scalar_type>
__global__ void zero_forces(G2G::vec_type<scalar_type,3>* forces,uint width,uint height)
{
  uint my_x = index_x(blockDim,blockIdx,threadIdx);
  uint my_y = blockIdx.y * blockDim.y + threadIdx.y;

  if (my_x < width && my_y < height) {
    forces[width*my_y + my_x].x = 0.0;
    forces[width*my_y + my_x].y = 0.0;
    forces[width*my_y + my_x].z = 0.0;
  }
}

//
// Reduce the partial Fock array into the first row
// Also, reduce corresponding partial energies over the block
//
template<class scalar_type>
__global__ void gpu_fock_reduce( double* fock, scalar_type* dens, double* energies, uint stride, uint depth, uint width )
{

  uint my_fock_ind = index_x(blockDim, blockIdx, threadIdx);
  uint tid = threadIdx.x;
  double my_fock = 0.0;

  //if (my_fock_ind < width) {
  //  for (uint i = 0; i < depth; i++) {
  //    my_partial_fock += fock[stride*i + my_fock_ind];
  //  }
  //  fock[my_fock_ind] = my_partial_fock;
  //}

  __shared__ double energies_sh[QMMM_REDUCE_BLOCK_SIZE];
  scalar_type my_dens = 0.0;
  if (my_fock_ind < width) {
    my_dens = dens[my_fock_ind];
    my_fock = fock[my_fock_ind];
  }
  energies_sh[tid] = my_fock * (double)my_dens;
  __syncthreads();

  if (tid < 64)
  {
    energies_sh[tid] += energies_sh[tid+64];
  }
  __syncthreads();

  if (tid < WARP_SIZE) { warpReduce<double>(energies_sh, tid); }
  if (tid == 0) { energies[blockIdx.x] = energies_sh[0]; }
}

// Double precision addition as an atomic operation
// Taken from the CUDA Toolkit Documentation
static __device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

