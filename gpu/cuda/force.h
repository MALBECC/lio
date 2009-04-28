
/**
 * Called for each atom
 */

__global__ void gpu_compute_forces(uint points, float* force_factors, float4* density_deriv, float4* forces, uint nucleii_count)
{
	uint atom = index_x(blockDim, blockIdx, threadIdx);
	
	bool valid_thread = (atom < nucleii_count);

	float4 atom_force = make_float4(0.0f,0.0f,0.0f,0.0f);

  __shared__ float factor_sh[FORCE_BLOCK_SIZE];
	
  for (uint point_base = 0; point_base < points; point_base += FORCE_BLOCK_SIZE) {
		__syncthreads();
		if (point_base + threadIdx.x < points) factor_sh[threadIdx.x] = force_factors[point_base + threadIdx.x];
		factor_sh[threadIdx.x] *= 4.0f;
		__syncthreads();

		if (valid_thread) {
      for (uint point_sub = 0; point_sub < FORCE_BLOCK_SIZE && (point_base + point_sub < points); point_sub++) {
        float4 density_deriv_local = density_deriv[COALESCED_DIMENSION(points) * atom + (point_base + point_sub)];
        atom_force += density_deriv_local * factor_sh[point_sub];
      }
		}
	}
  
	if (valid_thread) forces[atom] = atom_force;
}
