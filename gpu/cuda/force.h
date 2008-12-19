
/**
 * Called for each atom
 */

#define USE_SHARED_FORCE 0

// TODO: coalescear density_deriv y factor

__global__ void gpu_compute_forces(uint points, float* force_factors, float3* density_deriv, float3* forces, uint nucleii_count)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);
	uint atom = pos.x;
	
	bool valid_thread = (atom < nucleii_count);
	
	float3 atom_force = make_float3(0.0f,0.0f,0.0f);

  #if USE_SHARED_FORCE
  for (uint point = 0; point < points; point += RMM_BLOCK_SIZE_X * RMM_BLOCK_SIZE_Y) {
		__syncthreads();
		if (threadIdx.x == 0) factor_sh = 4.0f * force_factors[point];
		__syncthreads();

		if (valid_thread) atom_force = atom_force + density_deriv[point * nucleii_count + atom] * factor_sh;
	}
  #else
  
  __shared__ float factor_sh;

	for (uint point = 0; point < points; point++) {
		__syncthreads();
		if (threadIdx.x == 0) factor_sh = 4.0f * force_factors[point];
		__syncthreads();

		if (valid_thread) atom_force = atom_force + density_deriv[point * nucleii_count + atom] * factor_sh;
	}
  #endif

	if (valid_thread) forces[atom] = atom_force;
}
