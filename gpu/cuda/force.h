
/**
 * Called for each atom
 */

__global__ void gpu_compute_forces(uint points, float* force_factors, float3* density_deriv, float3* forces, uint nucleii_count)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);
	uint atom = pos.x;
	
	bool valid_thread = (atom < nucleii_count);

	float3 atom_force = make_float3(0.0f,0.0f,0.0f);

  __shared__ float factor_sh[FORCE_BLOCK_SIZE];
	
  for (uint point_base = 0; point_base < points; point_base += FORCE_BLOCK_SIZE) {
		__syncthreads();
		if (point_base + threadIdx.x < points) factor_sh[threadIdx.x] = force_factors[point_base + threadIdx.x];
		factor_sh[threadIdx.x] *= 4.0f;
		__syncthreads();

		for (uint point_sub = 0; point_sub < FORCE_BLOCK_SIZE && (point_base + point_sub < points); point_sub++) {
      uint point = point_base + point_sub;
      
			if (valid_thread) {
        float3 density_deriv_local = density_deriv[point * COALESCED_DIMENSION(nucleii_count) + atom];
        atom_force = atom_force + density_deriv_local * factor_sh[point_sub];
      }
		}
	}
  
	if (valid_thread) forces[atom] = atom_force;
}
