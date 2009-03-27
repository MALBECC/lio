#include "pot.h"

/**
 * gpu_compute_functions
 * 	grid: points 
 *  TODO: paralelizar por puntos (pueden cargar cosas comunes)
 */

#define NEW_DERIVS 0

template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* energy, float* factor, float4* point_weight_positions, uint points, float* rdm,
																		float* function_values, float3* gradient_values, float3* density_deriv, uint* nuc,
																		uint nucleii_count, uint4 functions)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;

	uint m = functions.w;
	
	float partial_density = 0.0f;

	__shared__ float rdm_sh[DENSITY_BLOCK_SIZE];
	
	bool valid_thread = (point < points);
	
	float4 point_weight_position;
	if (valid_thread) point_weight_position = point_weight_positions[point];
	//float3 point_position = make_float3(point_weight_position.x, point_weight_position.y, point_weight_position.z);
	float point_weight = point_weight_position.w;
	
	uint base = point * COALESCED_DIMENSION(m);
  uint grad_base = point * m;
	uint deriv_base = point * COALESCED_DIMENSION(nucleii_count);

  if (compute_derivs) {
    if (valid_thread) { for (uint i = 0; i < nucleii_count; i++) density_deriv[deriv_base + i] = make_float3(0.0f,0.0f,0.0f); }
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
					w += function_values[base + j + jj] * rdm_sh[jj];
				}
			}

			__syncthreads();
		}

		// TODO: coalescear RDM aca tambien
		if (compute_derivs) {
      if (valid_thread) {
        // TODO: esto se tiene que poder hacer mejor
        #if NEW_DERIVS
        uint j = 0;
        for (uint jj = 0; jj < m; jj++) {
          if (jj < functions.x) j = jj;
          else if (jj < functions.x + 3 * functions.y) j = (jj - functions.x)/3;
          else j = (jj - functions.x - functions.y * 3) / 6;

          /* TODO: nuc podria ser uno por cada jj o se podria cargar cada 3, 6 */
          density_deriv[deriv_base + nuc[j]] += gradient_values[grad_base + jj] * rdm[rdm_idx + jj] * w;
        }
        #else
        uint j = 0, jj = 0;
        for (; j < functions.x; j++, jj++) {
          density_deriv[deriv_base + nuc[j]] += gradient_values[grad_base + j] * rdm[rdm_idx + j] * w;
        }
        for (; j < functions.x + functions.y; j++, jj+=3) {
          uint deriv_idx = deriv_base + nuc[j];
          uint grad_idx = grad_base + jj;
          float3 partial_result = make_float3(0.0f,0.0f,0.0f);
          partial_result += gradient_values[grad_idx + 0] * rdm[rdm_idx + jj + 0];
          partial_result += gradient_values[grad_idx + 1] * rdm[rdm_idx + jj + 1];
          partial_result += gradient_values[grad_idx + 2] * rdm[rdm_idx + jj + 2];
          density_deriv[deriv_idx] += partial_result * w;
        }
        for (; j < functions.x + functions.y + functions.z; j++, jj+=6) {
          uint deriv_idx = deriv_base + nuc[j];
          uint grad_idx = grad_base + jj;
          float3 partial_result = make_float3(0.0f,0.0f,0.0f);
          partial_result += gradient_values[grad_idx + 0] * rdm[rdm_idx + jj + 0];
          partial_result += gradient_values[grad_idx + 1] * rdm[rdm_idx + jj + 1];
          partial_result += gradient_values[grad_idx + 2] * rdm[rdm_idx + jj + 2];
          partial_result += gradient_values[grad_idx + 3] * rdm[rdm_idx + jj + 3];
          partial_result += gradient_values[grad_idx + 4] * rdm[rdm_idx + jj + 4];
          partial_result += gradient_values[grad_idx + 5] * rdm[rdm_idx + jj + 5];
          density_deriv[deriv_idx] += partial_result * w;
        }
        #endif
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
		else {
			gpu_pot<true, false>(partial_density, exc, corr, y2a);
		}
		if (valid_thread) energy[point] = (partial_density * point_weight) * (exc + corr);
	}
	else {
		gpu_pot<false, true>(partial_density, exc, corr, y2a);
		if (valid_thread) factor[point] = point_weight * y2a;
	}
}	

/*********** test **************/
#define NCO_PER_BATCH 4

#define FUNCS_S_PER_BATCH 4
#define FUNCS_P_PER_BATCH 2
#define FUNCS_D_PER_BATCH 1

// max(6 * FUNCS_D_PER_BATCH, 3 * FUNCS_P_PER_BATCH, FUNCS_S_PER_BATCH)
#define FUNCS_PER_BATCH 6

template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density2(float* energy, float* factor, float4* point_weight_positions, uint points, float* rdm,
																		 uint* nuc, uint4 functions, uint spd, uint* contractions, float2* factor_ac)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;

	uint m = functions.w;
	
	float partial_density = 0.0f;
	
	/**** Load Point Information ****/
	bool valid_thread = (point < points);
	// TODO: ver el limite en la otra dimension

	__shared__ float w_sh[NCO_PER_BATCH * DENSITY_BLOCK_SIZE_X];
	__shared__ float f_sh[FUNCS_PER_BATCH * DENSITY_BLOCK_SIZE_X];
	__shared__ float rdm_sh[FUNCS_PER_BATCH * DENSITY_BLOCK_SIZE_X];

	float t, tg;
	float3 v;

	float4 point_weight_position;
	if (valid_thread) point_weight_position = point_weight_positions[point];	// TODO: mover a shared
	float3 point_position = make_float3(point_weight_position.x, point_weight_position.y, point_weight_position.z);
	float point_weight = point_weight_position.w;

	for (uint i = 0; i < gpu_nco; i += NCO_PER_BATCH) {
		for (uint ii = 0; ii < NCO_PER_BATCH; ii++) {
			w_sh[threadIdx.x * NCO_PER_BATCH + ii] = 0.0f;
		}	
	
		/********** compute functions **********/
		/* s functions */	
		for (uint j = 0; j < functions.x; j += FUNCS_S_PER_BATCH) {
			for (uint jj = 0; (jj < FUNCS_S_PER_BATCH) && (j+jj < functions.x); jj++) {
				compute_function<compute_derivs>(j+jj, point_position, contractions, factor_ac, nuc, spd, t, tg, v);
				f_sh[threadIdx.x * FUNCS_S_PER_BATCH + jj] = t;
			}
			
			/*** use functions ***/
			for (uint ii = 0; ii < NCO_PER_BATCH && (i+ii < gpu_nco); ii++) {
				for (uint jj = 0; (jj < FUNCS_S_PER_BATCH) && (j+jj < functions.x); jj++) {
					w_sh[threadIdx.x * NCO_PER_BATCH + ii] += f_sh[threadIdx.x * FUNCS_S_PER_BATCH + jj] * rdm[(i+ii) * m + (j + jj)];
				}
			}
		}

		/* p functions */	
		for (uint j = 0; j < functions.y; j += FUNCS_P_PER_BATCH) {
			for (uint jj = 0; (jj < FUNCS_P_PER_BATCH) && (j+jj < functions.y); jj++) {
				compute_function<compute_derivs>(functions.x + (j+jj), point_position, contractions, factor_ac, nuc, spd, t, tg, v);
				f_sh[(threadIdx.x * FUNCS_P_PER_BATCH + jj) * 3 + 0] = v.x * t;
				f_sh[(threadIdx.x * FUNCS_P_PER_BATCH + jj) * 3 + 1] = v.y * t;
				f_sh[(threadIdx.x * FUNCS_P_PER_BATCH + jj) * 3 + 2] = v.z * t;
			}
			
			/*** use functions ***/
			for (uint ii = 0; ii < NCO_PER_BATCH && (i+ii < gpu_nco); ii++) {
				for (uint jj = 0; (jj < FUNCS_P_PER_BATCH) && (j+jj < functions.y); jj++) {
					for (uint k = 0; k < 3; k++) {
						w_sh[threadIdx.x * NCO_PER_BATCH + ii] += f_sh[(threadIdx.x * FUNCS_P_PER_BATCH + jj) * 3 + k] * rdm[(i+ii) * m + (functions.x + (j+jj) * 3 + k)];
					}
				}
			}
		}

		/* d functions */	
		for (uint j = 0; j < functions.z; j += FUNCS_D_PER_BATCH) {
			for (uint jj = 0; (jj < FUNCS_D_PER_BATCH) && (j+jj < functions.z); jj++) {
				compute_function<compute_derivs>(functions.x + functions.y + (j+jj), point_position, contractions, factor_ac, nuc, spd, t, tg, v);

				float tx = t * v.x;
				float ty = t * v.y;
				float tz = t * v.z;

				f_sh[(threadIdx.x * FUNCS_D_PER_BATCH + jj) * 6 + 0] = tx * v.x * gpu_normalization_factor;
				f_sh[(threadIdx.x * FUNCS_D_PER_BATCH + jj) * 6 + 1] = ty * v.x;
				f_sh[(threadIdx.x * FUNCS_D_PER_BATCH + jj) * 6 + 2] = ty * v.y * gpu_normalization_factor;
				f_sh[(threadIdx.x * FUNCS_D_PER_BATCH + jj) * 6 + 3] = tz * v.x;
				f_sh[(threadIdx.x * FUNCS_D_PER_BATCH + jj) * 6 + 4] = tz * v.y;
				f_sh[(threadIdx.x * FUNCS_D_PER_BATCH + jj) * 6 + 5] = tz * v.z * gpu_normalization_factor;
			}
			
			/*** use functions ***/
			for (uint ii = 0; ii < NCO_PER_BATCH && (i+ii < gpu_nco); ii++) {
				for (uint jj = 0; (jj < FUNCS_D_PER_BATCH) && (j+jj < functions.z); jj++) {
					for (uint k = 0; k < 6; k++) {
						w_sh[threadIdx.x * NCO_PER_BATCH + ii] += f_sh[(threadIdx.x * FUNCS_D_PER_BATCH + jj) * 6 + k] *
																												rdm[(i+ii) * m + (functions.x + functions.y * 3 + (j+jj) * 6 + k)];
					}
				}
			}
		}

		/* density contribution */
		for (uint ii = 0; ii < NCO_PER_BATCH && (i+ii < gpu_nco); ii++) {
			// see if this can be loaded just once
			partial_density += w_sh[threadIdx.x * NCO_PER_BATCH + ii] * w_sh[threadIdx.x * NCO_PER_BATCH + ii];
		}
	}
	partial_density *= 2.0f;

	float exc, corr, y2a;
	if (compute_energy) {
		if (compute_derivs) {
			gpu_pot<true, true>(partial_density, exc, corr, y2a);
			if (valid_thread) factor[point] = point_weight * y2a;
		}
		else {
			gpu_pot<true, false>(partial_density, exc, corr, y2a);
		}
		if (valid_thread) energy[point] = (partial_density * point_weight) * (exc + corr);
	}
	else {
		gpu_pot<false, true>(partial_density, exc, corr, y2a);
		if (valid_thread) factor[point] = point_weight * y2a;
	}
}	
