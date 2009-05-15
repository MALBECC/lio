template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* const energy, float* const factor, const float* const point_weights,
  uint points, const float* rdm, const float* const function_values, uint m, float* w)
{
  dim3 pos = index(blockDim, blockIdx, threadIdx);
  uint point = pos.x;

  float partial_density = 0.0f;

  __shared__ float rdm_sh[DENSITY_BLOCK_SIZE];

  bool valid_thread = (point < points);

  float point_weight;
  if (valid_thread) point_weight = point_weights[point];

  /***** compute density ******/
  for (uint i = 0; i < gpu_nco; i++) {
    float w_local = 0.0f;
    for (uint j = 0; j < m; j += DENSITY_BLOCK_SIZE) {
      if (j + threadIdx.x < m) rdm_sh[threadIdx.x] = rdm[j + threadIdx.x];

      __syncthreads();

      if (valid_thread) {
        for (uint jj = 0; jj < DENSITY_BLOCK_SIZE && (j + jj < m); jj++) {
          w_local += function_values[COALESCED_DIMENSION(points) * (j + jj) + point] * rdm_sh[jj];
        }
      }
      __syncthreads();
    }

    if (compute_derivs && valid_thread) {
      w[point] = w_local;
      w += COALESCED_DIMENSION(points);
    }
    partial_density += w_local * w_local;

    rdm += COALESCED_DIMENSION(m);
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
