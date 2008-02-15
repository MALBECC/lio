
/**
 * Called for each atom
 */
template <uint grid_n, const uint* const curr_layers>
__global__ void calc_forces(const float3* atom_positions, const uint* types, const uint atoms_n, uint3 num_funcs,
													 const uint* nuc, const uint* contractions, bool normalize, const float2* factor_ac,
													 const float* factors, const float3* dds, float3* forces)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);
	uint atom_i = pos.x;
	
	bool valid_thread = true;
	if (atom_i >= atoms_n) valid_thread = false;
	
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);

	//float3 lost_numbers = make_float3(0.0f,0.0f,0.0f);	
	float3 atom_force = make_float3(0.0f,0.0f,0.0f);
	
	__shared__ float factor_local;

	for (uint atom_j = 0; atom_j < atoms_n; atom_j++) {
		uint atom_j_type = types[atom_j];
		uint atom_j_layers = curr_layers[atom_j_type];

		for (uint layer_atom_j = 0; layer_atom_j < atom_j_layers; layer_atom_j++) {
			for (uint point_atom_j = 0; point_atom_j < grid_n; point_atom_j++) {
				uint factor_idx = index_from3d(energySize, dim3(atom_j, layer_atom_j, point_atom_j));
				
				__syncthreads();				 // por si el escritor se adelanta a los lectores				

				/* cache into local memory */
				if (threadIdx.x == 0 && threadIdx.y == 0) {
					factor_local = factors[factor_idx];
				}
				
			  /* fill local variables from local cache */
				__syncthreads();	// por si los lectores se adelantan al escritor
				
				if (valid_thread) {
					/*float3 temp_force = ((dds[factor_idx * atoms_n + atom_i] * 4.0f * factor_local) - lost_numbers);
					float3 temp_force2 = atom_force + temp_force;
					lost_numbers = (temp_force2 - atom_force) - temp_force;
					atom_force = temp_force2;*/
					
					atom_force = atom_force + (dds[factor_idx * atoms_n + atom_i] * 4.0f * factor_local);
				}
			}			
		}
	}

	if (valid_thread) { forces[atom_i] = atom_force; }
}

template <uint grid_n>
	void calc_forces_cpu(const float3* atom_positions, const uint* types, const uint atoms_n, uint3 num_funcs,
										 const uint* nuc, const uint* contractions, bool normalize, const float2* factor_ac,
										 const float* factors, const float3* dds, double3* forces, const uint* curr_layers)
{
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);
	
	for (uint atom_i = 0; atom_i < atoms_n; atom_i++) {
		double3 atom_force;
		atom_force.x = atom_force.y = atom_force.z = 0.0f;
		
		for (uint atom_j = 0; atom_j < atoms_n; atom_j++) {
		uint atom_j_type = types[atom_j];
		uint atom_j_layers = curr_layers[atom_j_type];

			for (uint layer_atom_j = 0; layer_atom_j < atom_j_layers; layer_atom_j++) {
				for (uint point_atom_j = 0; point_atom_j < grid_n; point_atom_j++) {
					uint factor_idx = index_from3d(energySize, dim3(atom_j, layer_atom_j, point_atom_j));

					atom_force = atom_force + (dds[factor_idx * atoms_n + atom_i] * factors[factor_idx] * 4.0f);
				}
			}
		}
		
		forces[atom_i] = atom_force;
	}
}
