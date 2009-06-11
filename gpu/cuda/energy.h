template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* const energy, float* const factor, const float* const point_weights,
  uint points, const float* rdm, const float* const function_values, uint m, float* w)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);

  float partial_density = 0.0f;

  bool valid_thread = (point < points);
  float point_weight;
  if (valid_thread) point_weight = point_weights[point];

  __shared__ float rdm_sh[NCO_BATCH_SIZE];
  __shared__ float w_sh[NCO_BATCH_SIZE + 1][DENSITY_BLOCK_SIZE];

  /***** compute density ******/
  for (uint i = 0; i < gpu_nco; i += NCO_BATCH_SIZE) {
    for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) w_sh[ii][threadIdx.x] = 0.0f;

    for (uint j = 0; j < m; j++) {
      __syncthreads();
      if ((i + threadIdx.x < NCO_BATCH_SIZE) && (i + threadIdx.x < gpu_nco)) rdm_sh[threadIdx.x] = rdm[COALESCED_DIMENSION(gpu_nco) * j + (i + threadIdx.x)];
      __syncthreads();

      if (valid_thread) {
        float f = function_values[COALESCED_DIMENSION(points) * j + point];
        for (uint ii = 0; ii < NCO_BATCH_SIZE && (i + ii < gpu_nco); ii++) {
          w_sh[ii][threadIdx.x] += f * rdm_sh[ii];
        }
      }
    }

    for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) { partial_density += w_sh[ii][threadIdx.x] * w_sh[ii][threadIdx.x]; }
  }
  partial_density *= 2.0f;

  /***** compute energy / factor *****/
  float y2a, exc_corr;
  if (compute_energy) {
    if (compute_derivs) {
			gpu_pot<true, true>(partial_density, exc_corr, y2a);
      if (valid_thread) factor[point] = point_weight * y2a;
		}
		else gpu_pot<true, false>(partial_density, exc_corr, y2a);

    if (valid_thread) energy[point] = (partial_density * point_weight) * exc_corr;
	}
	else {
		gpu_pot<false, true>(partial_density, exc_corr, y2a);
    if (valid_thread) factor[point] = point_weight * y2a;
  } 
}
