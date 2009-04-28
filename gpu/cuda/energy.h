#include "pot.h"

/**
 * gpu_compute_functions
 * 	grid: points 
 */
template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* energy, float* factor, float* point_weights, uint points, float* rdm,
  float* function_values, uint m, float* w)
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
    uint rdm_idx = i * COALESCED_DIMENSION(m);

    float w_local = 0.0f;
    for (uint j = 0; j < m; j += DENSITY_BLOCK_SIZE) {
      if (j + threadIdx.x < m) rdm_sh[threadIdx.x] = rdm[rdm_idx + j + threadIdx.x];

      __syncthreads();

      if (valid_thread) {
        for (uint jj = 0; jj < DENSITY_BLOCK_SIZE && (j + jj < m) ; jj++) {
          w_local += function_values[COALESCED_DIMENSION(points) * (j + jj) + point] * rdm_sh[jj];
        }
      }

      __syncthreads();
    }

    if (compute_derivs && valid_thread) w[COALESCED_DIMENSION(points) * i + point] = w_local;
    partial_density += w_local * w_local;
  }
  partial_density *= 2.0f;

  /***** compute energy / factor *****/
  float exc, corr, y2a;
  if (compute_energy) {
    if (compute_derivs) {
			gpu_pot<true, true>(partial_density, exc, corr, y2a);
      if (valid_thread) factor[point] = point_weight * y2a;
		}
		else gpu_pot<true, false>(partial_density, exc, corr, y2a);

    if (valid_thread) energy[point] = (partial_density * point_weight) * (exc + corr);
	}
	else {
		gpu_pot<false, true>(partial_density, exc, corr, y2a);
    if (valid_thread) factor[point] = point_weight * y2a;
  } 
}

__global__ void gpu_compute_density_derivs(uint points, float* rdmt, float4* gradient_values, float4* density_deriv, uint* nuc,
                                           uint nucleii_count, uint m, float* w)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);
  bool valid_thread = (point < points);
  
  __shared__ uint nuc_sh[DENSITY_DERIV_BLOCK_SIZE];
  __shared__ float rdm_sh[DENSITY_DERIV_BLOCK_SIZE];

  if (valid_thread) { for (uint i = 0; i < nucleii_count; i++) density_deriv[COALESCED_DIMENSION(points) * i + point] = make_float4(0.0f,0.0f,0.0f,0.0f); }

  for (uint j = 0; j < m; j += DENSITY_DERIV_BLOCK_SIZE) {
    if (j + threadIdx.x < m) nuc_sh[threadIdx.x] = nuc[j + threadIdx.x];

    for (uint jj = 0; jj < DENSITY_DERIV_BLOCK_SIZE && (j + jj < m); jj++) {
      float wrdm = 0.0f;

      for (uint i = 0; i < gpu_nco; i += DENSITY_DERIV_BLOCK_SIZE) {
        if (i + threadIdx.x < gpu_nco) rdm_sh[threadIdx.x] = rdmt[COALESCED_DIMENSION(gpu_nco) * (j + jj) + i + threadIdx.x];

        __syncthreads();
        if (valid_thread) {
          for (uint ii = 0; ii < DENSITY_DERIV_BLOCK_SIZE  && (i + ii < gpu_nco); ii++) {
            wrdm += rdm_sh[ii] * w[COALESCED_DIMENSION(points) * (i + ii) + point];
          }
        }
        __syncthreads();
      }

      if (valid_thread) {
        uint this_nuc = nuc_sh[jj]; // Parece ser necesario para que el compilador coalescee la escritura de abajo
        density_deriv[COALESCED_DIMENSION(points) * this_nuc + point] += gradient_values[COALESCED_DIMENSION(points) * (j + jj) + point] * wrdm;
      }
    }

    __syncthreads();
  }
}
