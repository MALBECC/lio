
// This function is called for each (i,j) in RMM, using a thread for each point.
// TODO: This wastes half the threads, may be we can build a grid without
// ignored blocks.
template <class scalar_type, bool check_pos>
__global__ void gpu_update_rmm(const scalar_type* __restrict__ factors,
                               int points, scalar_type* rmm,
                               const scalar_type* __restrict__ function_values,
                               int m) {
  //    if (blockIdx.x * blockDim.x > blockIdx.y * blockDim.y) return;
  int i, j, first_fi, first_fj;
  // There's more than one block; do the math to get position
  if (check_pos) {
    // Figure out where our block is in the lower triangle
    // We get 1D index k = blockIdx.x; column 1 value of k in row j is always
    // k_1 = j*(j+1); solving for j gives determinant 1+8*k_1
    // Thus, 1+8*k_1 must be square of (odd) integer; take sqrt(1+8*k) and get
    // first odd integer below it - get row and column from there
    int n = sqrtf(1.0f + 8.0f * blockIdx.x);
    n -= (1 - n % 2);
    int block_j = (n - 1) / 2;
    int block_i = blockIdx.x - (block_j + 1) * block_j / 2;

    first_fi = block_i * blockDim.x;
    first_fj = block_j * blockDim.y;
    i = first_fi + threadIdx.x;
    j = first_fj + threadIdx.y;
    // There's only one block; we don't need to take the square root
  } else {
    uint3 pos = index(blockDim, blockIdx, threadIdx);

    i = pos.x;  // Column
    j = pos.y;  // Row
    first_fi = blockIdx.x * blockDim.x;
    first_fj = blockIdx.y * blockDim.y;
  }

  // Keep only the lower triangle.
  bool valid_thread = ( (i < m) && (j < m) && (i <= j) );

  // This stores the RMM section to be calculated.
  scalar_type rmm_local = 0.0f;

  __shared__ scalar_type // Fi[point][i]
      functions_i_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY + 1];
  __shared__ scalar_type // Fj[point][j]
      functions_j_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY + 1];
  __shared__ scalar_type // factor[point] May have shared bank conflics
      factor_local[RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY];

  int inc = RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY;
  // absolute threadId inside block
  int abs_threadIdx = threadIdx.y * blockDim.x + threadIdx.x;
  bool valid_fi_thread = (first_fi + threadIdx.y) < m;
  bool valid_fj_thread = (first_fj + threadIdx.y) < m;
  for (int point_base = 0; point_base < points; point_base += inc) {
    __syncthreads();

    /* all threads load a point */
    if (point_base + abs_threadIdx < points)
      factor_local[abs_threadIdx] = factors[point_base + abs_threadIdx];

    int last_point = point_base + inc;

    #pragma unroll 16
    for (int point = point_base; point < last_point;
         point += RMM_BLOCK_SIZE_XY) {
      if (point < points) {
        /* every RMM_BLOCK_SIZE_X iterations, Fi and Fj get filled with
         * RMM_BLOCK_SIZE_Y functions,
         * for RMM_BLOCK_SIZE_X different points */
        __syncthreads();
        bool valid_point = point + threadIdx.x < points;
        bool validFi = valid_point * valid_fi_thread;
        bool validFj = valid_point * valid_fj_thread;
        int function_values_fi_index =
            COALESCED_DIMENSION(points) * (first_fi + threadIdx.y) +
            (point + threadIdx.x);
        int function_values_fj_index =
            COALESCED_DIMENSION(points) * (first_fj + threadIdx.y) +
            (point + threadIdx.x);
        int factor_local_fi_index = (point - point_base) + threadIdx.x;

        /* on blocks with some invalid threads, the computation contributes 0 to
         * rmm_local (
         * this avoids instruction serialization).
         * We avoid if-elses using the boolean predicates of those if to
         * multiply by 0 or 1 if needed.
         * We also do it on the array subindeces to avoid segfaulting.
         */
        if (validFi) {
           scalar_type fi_times_factor =
               function_values[function_values_fi_index] *
               factor_local[factor_local_fi_index];

           functions_i_local[threadIdx.x][threadIdx.y] = fi_times_factor;
        }

        if (validFj) {
        functions_j_local[threadIdx.x][threadIdx.y] =
            function_values[function_values_fj_index];
        }

        __syncthreads();
        if (validFj && validFi) {
        for (int point_sub = 0; point_sub < RMM_BLOCK_SIZE_XY; point_sub++) {
          rmm_local += functions_i_local[point_sub][threadIdx.x] *
                       functions_j_local[point_sub][threadIdx.y];
        }
        }
      }
    }
  }

  if (valid_thread) {
     rmm[COALESCED_DIMENSION(m) * j + i] = rmm_local;
  }
}
