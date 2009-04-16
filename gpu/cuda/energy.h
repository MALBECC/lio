#include "pot.h"

/**
 * gpu_compute_functions
 * 	grid: points 
 */

template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* energy, float* factor, float* point_weights, uint points, float* rdm,
																		float* function_values, float4* gradient_values, float4* density_deriv, uint* nuc,
																		uint nucleii_count, uint4 functions)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;

	uint m = functions.w;
	
	float partial_density = 0.0f;

	__shared__ float rdm_sh[DENSITY_BLOCK_SIZE];
  __shared__ uint nuc_sh[DENSITY_BLOCK_SIZE];
	
	bool valid_thread = (point < points);
	
	float point_weight;
	if (valid_thread) point_weight = point_weights[point];

  if (compute_derivs) {
    if (valid_thread) { for (uint i = 0; i < nucleii_count; i++) density_deriv[COALESCED_DIMENSION(points) * i + point] = make_float4(0.0f,0.0f,0.0f,0.0f); }
  }

	/* density */	
	for (uint i = 0; i < gpu_nco; i++) {
		uint rdm_idx = i * COALESCED_DIMENSION(m);

		float w = 0.0f;
		/* coalesced access to RDM */
		for (uint j = 0; j < m; j += DENSITY_BLOCK_SIZE) {
			if (j + threadIdx.x < m) rdm_sh[threadIdx.x] = rdm[rdm_idx + j + threadIdx.x];

			__syncthreads();
			
			if (valid_thread) {
				for (uint jj = 0; jj < DENSITY_BLOCK_SIZE && (j + jj < m); jj++) {
          uint base = COALESCED_DIMENSION(points) * (j + jj) + point;
					w += function_values[base] * rdm_sh[jj];
				}
			}

			__syncthreads();
		}

		if (compute_derivs) {
      for (uint j = 0; j < m; j += DENSITY_BLOCK_SIZE) {
        if (j + threadIdx.x < m) rdm_sh[threadIdx.x] = rdm[rdm_idx + j + threadIdx.x];
        if (j + threadIdx.x < functions.x + functions.y + functions.z) nuc_sh[threadIdx.x] = nuc[j + threadIdx.x];
        __syncthreads();

        if (valid_thread) {
          for (uint jj = 0; jj < DENSITY_BLOCK_SIZE && (j + jj < m); jj++) {
            uint k;
            if ((j + jj) < functions.x) k = (j+jj);
            else if ((j + jj) < functions.x + functions.y * 3) k = ((j + jj) - functions.x) / 3 + functions.x;
            else k = ((j + jj) - (functions.x + functions.y * 3)) / 6 + (functions.x + functions.y);
            density_deriv[COALESCED_DIMENSION(points) * nuc_sh[k] + point] +=
                gradient_values[COALESCED_DIMENSION(points) * (j + jj) + point] * (rdm_sh[jj] * w);
          }
        }

        __syncthreads();
      }
		}
		partial_density += w * w;
	}
	partial_density *= 2.0f;

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
