#include "pot.h"

/**
 * gpu_compute_functions
 * 	grid: points 
 *  TODO: paralelizar por puntos (pueden cargar cosas comunes)
 */
template<bool compute_energy> __global__ void gpu_compute_density(float* result, float4* point_weight_positions, uint points, float* rdm, float* function_values, uint m)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;
	
	float partial_density = 0.0f;
	
	/**** Load Point Information ****/
	if (point >= points) return; 	// esto se puede evitar computando basura
	float4 point_weight_position = point_weight_positions[point];
	float3 point_position = make_float3(point_weight_position.x, point_weight_position.y, point_weight_position.z);
	float point_weight = point_weight_position.w;
	
	uint base = point * m;

	/* density */	
	for (uint i = 0; i < gpu_nco; i++) {
		float w = 0.0f;
		for (uint j = 0; j < m; j++) {
			w += rdm[i * m + j] * function_values[base + j];	
		}
		partial_density += w * w;
	}
	partial_density *= 2.0f;

	float exc, corr, y2a;
	if (compute_energy) {
		gpu_pot<false>(partial_density, exc, corr, y2a);
		result[point] = (partial_density * point_weight) * (exc + corr);
	}
	else {
		gpu_pot<true>(partial_density, exc, corr, y2a);
		result[point] = point_weight * y2a;
	}
}	
