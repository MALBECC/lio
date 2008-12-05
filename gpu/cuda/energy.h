#include "pot.h"

/**
 * gpu_compute_functions
 * 	grid: points 
 *  TODO: paralelizar por puntos (pueden cargar cosas comunes)
 */
template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* energy, float* factor, float4* point_weight_positions, uint points, float* rdm,
																		float* function_values, float3* gradient_values, float3* density_deriv, uint* nuc,
																		uint nucleii_count, uint4 functions)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;

	uint m = functions.w;
	
	float partial_density = 0.0f;
	
	/**** Load Point Information ****/
	if (point >= points) return;
	float4 point_weight_position = point_weight_positions[point];
	float3 point_position = make_float3(point_weight_position.x, point_weight_position.y, point_weight_position.z);
	float point_weight = point_weight_position.w;
	
	uint base = point * m;
	uint deriv_base = point * nucleii_count;

	if (compute_derivs) {
		for (uint atom = 0; atom < nucleii_count; atom++) density_deriv[deriv_base + atom] = make_float3(0.0f,0.0f,0.0f);
	}

	/* density */	
	for (uint i = 0; i < gpu_nco; i++) {
		float w = 0.0f;
		for (uint j = 0; j < m; j++) {
			w += function_values[base + j] * rdm[i * m + j];	
		}

		if (compute_derivs) {
			uint jj = 0;
			for (uint j = 0; j < functions.x; j++, jj++) {
				density_deriv[deriv_base + nuc[j]] += gradient_values[base + j] * rdm[i * m + j] * w;
			}
			for (uint j = functions.x; j < functions.x + functions.y; j++, jj+=3) {
				uint nuc_j = nuc[j];
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 0] * rdm[i * m + jj + 0] * w;	
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 1] * rdm[i * m + jj + 1] * w;	
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 2] * rdm[i * m + jj + 2] * w;
			}
			for (uint j = functions.x + functions.y; j < functions.x + functions.y + functions.z; j++, jj+=6) {
				uint nuc_j = nuc[j];
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 0] * rdm[i * m + jj + 0] * w;	
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 1] * rdm[i * m + jj + 1] * w;	
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 2] * rdm[i * m + jj + 2] * w;	
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 3] * rdm[i * m + jj + 3] * w;	
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 4] * rdm[i * m + jj + 4] * w;	
				density_deriv[deriv_base + nuc_j] += gradient_values[base + jj + 5] * rdm[i * m + jj + 5] * w;
			}
		}

		partial_density += w * w;
	}
	partial_density *= 2.0f;

	float exc, corr, y2a;
	if (compute_energy) {
		if (compute_derivs) {
			gpu_pot<true, true>(partial_density, exc, corr, y2a);
			factor[point] = point_weight * y2a;
		}
		else {
			gpu_pot<true, false>(partial_density, exc, corr, y2a);
		}
		energy[point] = (partial_density * point_weight) * (exc + corr);
	}
	else {
		gpu_pot<false, true>(partial_density, exc, corr, y2a);
		factor[point] = point_weight * y2a;
	}
}	
