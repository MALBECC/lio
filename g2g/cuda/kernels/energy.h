template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* const energy, float* const factor, const float* const point_weights,
  uint points, const float* rdm, const float* const function_values, uint m, float* w)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);

  float partial_density = 0.0f;

  bool valid_thread = (point < points);
  float point_weight;
  if (valid_thread) point_weight = point_weights[point];

  __shared__ float rdm_sh[ENERGY_BATCH_SIZE];

  /***** compute density ******/
  for (uint i = 0; i < m; i++) {
    float fi = function_values[COALESCED_DIMENSION(points) * i + point];
    float w = 0.0f;

    for (uint bj = i; bj < m; bj += ENERGY_BATCH_SIZE) {
      __syncthreads();
      if (threadIdx.x < ENERGY_BATCH_SIZE) {
        if (bj + threadIdx.x < m) rdm_sh[threadIdx.x] = rdm[COALESCED_DIMENSION(m) * i + (bj + threadIdx.x)];
        else rdm_sh[threadIdx.x] = 0.0f;
      }
      __syncthreads();

      if (valid_thread) {
        for (uint j = 0; j < ENERGY_BATCH_SIZE && (bj + j) < m; j++) {
          float fj = function_values[COALESCED_DIMENSION(points) * (bj + j) + point];
          w += rdm_sh[j] * fj;
        }
      }
    }

    partial_density += fi * w;
  }

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
