
/**
 * Called for each atom
 */

#define USE_SHARED_FORCE 0

// TODO: coalescear density_deriv -> nucleii_count multiplo de 16

__global__ void gpu_compute_forces(uint points, float* force_factors, float3* density_deriv, float3* forces, uint nucleii_count)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);
	uint atom = pos.x;
	
	bool valid_thread = (atom < nucleii_count);
	
	float3 atom_force = make_float3(0.0f,0.0f,0.0f);

/*  #if USE_SHARED_FORCE
	__shared__ float factor_sh[FORCE_BLOCK_SIZE];
	
  for (uint point_base = 0; point_base < points; point_base += FORCE_BLOCK_SIZE) {
		__syncthreads();
		if (threadIdx.x == 0) {
			if (point_base + threadIdx.x < points) factor_sh[threadIdx.x] = force_factors[point_base + threadIdx.x];
			factor_sh *= 4.0f;
		}
		__syncthreads();

		for (uint point_sub = 0; point_sub < FORCE_BLOCK_SIZE && (point_base + point_sub < points); point_sub++) {
			if (valid_thread)
				atom_force = atom_force + density_deriv[point * nucleii_count + atom] * factor_sh[point_sub];
		}
	}
  #else*/
  
  __shared__ float factor_sh;

	for (uint point = 0; point < points; point++) {
		__syncthreads();
		if (threadIdx.x == 0) factor_sh = 4.0f * force_factors[point];
		__syncthreads();

		if (valid_thread) {
      _EMU(printf("atom: %i point: %i %.12e %.12e %.12e\n", atom, point, density_deriv[point * nucleii_count + atom].x, density_deriv[point * nucleii_count + atom].y, density_deriv[point * nucleii_count + atom].z));
      atom_force = atom_force + density_deriv[point * nucleii_count + atom] * factor_sh;
    }
	}
  //#endif

	if (valid_thread) forces[atom] = atom_force;
}
