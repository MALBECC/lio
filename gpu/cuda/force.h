
/**
 * Called for each atom
 */
__global__ void gpu_compute_forces(uint points, float* force_factors, float3* density_deriv, float3* forces, uint nucleii_count)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);
	uint atom = pos.x;
	
	bool valid_thread = true;
	if (atom >= nucleii_count) valid_thread = false;
	
	float3 atom_force = make_float3(0.0f,0.0f,0.0f);
	
	__shared__ float factor_sh;

	for (uint point = 0; point < points; point++) {
		__syncthreads();
		if (threadIdx.x == 0) factor_sh = 4.0f * force_factors[point];
		__syncthreads();

		if (valid_thread) atom_force = atom_force + density_deriv[point * nucleii_count + atom] * factor_sh;
	}

	if (valid_thread) forces[atom] = atom_force;
}
