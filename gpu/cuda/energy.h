#include "pot.h"

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

#if 0
// para BLOCK_SIZE de 32, va mejor con 21 (eso no coalescea)
#define NCO_BATCH_SIZE 16

/* TODO: coalescear contractions y demas */
template<bool do_forces>
__device__ __host__ void compute_function_(float3 point_position, uint contractions,
  float* factor_a, float* factor_c, uint nuc, float& fj, uint4 functions, uint idx)
{
	float3 atom_nuc_position = gpu_atom_positions[nuc]; // TODO: ver si al usar memoria compartida para esto, pago menos precio por todos los misses
	float3 v = point_position - atom_nuc_position;
	float dist = length2(v);

	float t = 0.0f;

	for (uint contraction = 0; contraction < contractions; contraction++) {
		float t0 = expf(-(factor_a[contraction] * dist)) * factor_c[contraction];
		t += t0;
		//if (do_forces) tg += t0 * curr_factor_ac.x;
	}

  if (idx < functions.x) { fj = t; }
  else if (idx < (functions.x + functions.y * 3)) {
    uint p_idx = (idx - functions.x) % 3;
    switch(p_idx) {
      case 0: fj = v.x * t; break;
      case 1: fj = v.y * t; break;
      case 2: fj = v.z * t; break;
    }
  }
  else {
    uint d_idx = (idx - functions.x - functions.y * 3) % 6;
    switch(d_idx) {
      case 0: fj = t * v.x * v.x * gpu_normalization_factor; break;
      case 1: fj = t * v.y * v.x;                            break;
      case 2: fj = t * v.y * v.y * gpu_normalization_factor; break;
      case 3: fj = t * v.z * v.x;                            break;
      case 4: fj = t * v.z * v.y;                            break;
      case 5: fj = t * v.z * v.z * gpu_normalization_factor; break;
    }
  }
}

template<bool compute_energy, bool compute_derivs>
__global__ void gpu_compute_density(float* energy, float* factor, float4* point_weight_positions, uint points, float* rdm, uint rdm_width,
  uint4 functions, float* w, uint2* nuc_contractions, float2* factor_ac)
{
  dim3 pos = index(blockDim, blockIdx, threadIdx);
  uint point = pos.x;
  uint m = functions.w;

  float partial_density = 0.0f;

  __shared__ float rdm_sh[DENSITY_BLOCK_SIZE];
  __shared__ uint nuc_sh[DENSITY_BLOCK_SIZE];
  __shared__ uint contractions_sh[DENSITY_BLOCK_SIZE];
  __shared__ float factor_a_sh[DENSITY_BLOCK_SIZE][MAX_CONTRACTIONS];
  __shared__ float factor_c_sh[DENSITY_BLOCK_SIZE][MAX_CONTRACTIONS];
  __shared__ float w_sh[NCO_BATCH_SIZE][DENSITY_BLOCK_SIZE + 1];   // TODO: poner en 0

  #pragma unroll 16
  for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) w_sh[ii][threadIdx.x] = 0.0f;

  bool valid_thread = (point < points);

  float4 point_weight_position;
  if (valid_thread) point_weight_position = point_weight_positions[point];
  float point_weight = point_weight_position.w;
  float3 point_position = to_float3(point_weight_position);

  for (uint i = 0; i < gpu_nco; i += NCO_BATCH_SIZE) {

    for (uint j = 0; j < m; j += DENSITY_BLOCK_SIZE) {
      __syncthreads();

      if (j + threadIdx.x < m) {
        uint2 nuc_contraction = nuc_contractions[j + threadIdx.x];
        nuc_sh[threadIdx.x] = nuc_contraction.x;
        contractions_sh[threadIdx.x] = nuc_contraction.y;
        for (uint contraction = 0; contraction < contractions_sh[threadIdx.x]; contraction++) {
          float2 factor_ac_local = factor_ac[COALESCED_DIMENSION(m) * contraction + (j + threadIdx.x)];
          factor_a_sh[threadIdx.x][contraction] = factor_ac_local.x;
          factor_c_sh[threadIdx.x][contraction] = factor_ac_local.y;
        }
      }

      #pragma unroll 16
      for (uint jj = 0; jj < DENSITY_BLOCK_SIZE; jj++) {
        if (j + jj < m) {
          float fj;
          compute_function_<compute_derivs>(point_position, contractions_sh[jj], factor_a_sh[jj], factor_c_sh[jj], nuc_sh[jj], fj, functions, j + jj);

          /* si threadIdx.x + i >= gpu_nco, esto carga ceros */
          if (i + threadIdx.x < gpu_nco) rdm_sh[threadIdx.x] = rdm[rdm_width * (j + jj) + (i + threadIdx.x)];
          //if (i + DENSITY_BLOCK_SIZE + threadIdx.x < gpu_nco) rdm_sh[DENSITY_BLOCK_SIZE + threadIdx.x] = rdm[rdm_width * (j + jj) + (i + DENSITY_BLOCK_SIZE + threadIdx.x)];
          //else rdm_sh[threadIdx.x] = 0.0f;

          __syncthreads();

          #pragma unroll 16
          for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) {
            if (i + ii < gpu_nco) w_sh[ii][threadIdx.x] += fj * rdm_sh[ii];
          }

          __syncthreads();
        }
      }
    }

    #pragma unroll 16
    for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) {
      partial_density += w_sh[ii][threadIdx.x] * w_sh[ii][threadIdx.x];
    }
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
#endif
