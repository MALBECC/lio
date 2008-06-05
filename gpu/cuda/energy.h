/**
 * Main Energy Kernel
 */

template <bool compute_forces, unsigned int grid_n, uint curr_layers>
__global__ void energy_kernel(float* energy, const uint atoms_n, uint nco, uint3 num_funcs,
		const uint2* nuc_contractions, bool normalize, const float2* factor_ac,
		const float* rmm, float* all_functions,  uint Ndens, float* output_factor, float3* dd, float3* Fg, float3* w3,
    bool compute_energy, bool update_rmm)
{
	const uint& m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;				 
	
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);
	dim3 pos2d = index(blockDim, blockIdx, threadIdx);
	
	const uint atom_i = pos2d.x;
	const uint point_atom_i = pos2d.y;
	
	// determine if thread is valid
	bool valid_thread = true;
	if (atom_i >= atoms_n) valid_thread = false;
	if (point_atom_i >= grid_n) valid_thread = false;
	if (!valid_thread) return;
	
	// not loaded from shared
	uint atom_i_type = gpu_types[atom_i]; // constant memory
	float3 rel_point_position = gpu_point_positions[point_atom_i];	// constant memory
	float wang_point_i = gpu_wang[point_atom_i]; // constant memory

	float rm = rm_factor[atom_i_type]; // constant memory
	uint atom_i_layers;
	if (curr_layers == GPU_LAYERS_1) atom_i_layers = gpu_layers_1[atom_i_type];	// constant memory
	else atom_i_layers = gpu_layers_2[atom_i_type];	// constant memory
	
	float3 atom_i_position = gpu_atom_positions[atom_i]; // constant memory
	float tmp0 = (PI / (atom_i_layers + 1.0f));	
	
	//_EMU(printf("atom %i type %i layers %i\n", atom_i, atom_i_type, atom_i_layers));

	// compute
	for (uint layer_atom_i = 0; layer_atom_i < atom_i_layers; layer_atom_i++) {
		dim3 pos = dim3(pos2d.x, layer_atom_i, pos2d.y);
		uint big_index = index_from3d(energySize, pos);

		float tmp1 = tmp0 * (layer_atom_i + 1.0f);
		float x = cosf(tmp1);
		float r1 = rm * (1.0f + x) / (1.0f - x);
		float w = tmp0 * fabsf(sinf(tmp1));
		float wrad = w * (r1 * r1) * rm * 2.0f / ((1.0f - x) * (1.0f - x));
		//if (point_atom_i == 0) printf("atom: %i layer: %i rm: %.12e tmp0: %.12e tmp1: %.12e x: %.12e\n", atom_i, layer_atom_i, rm, tmp0, tmp1, x);

		float3 abs_point_position = atom_i_position + rel_point_position * r1;
		float integration_weight = wang_point_i * wrad;

		float exc_curr = 1.0f;
		float corr_curr = 1.0f;
		float dens = 0.0f;
		float y2a = 0.0f;
		// float y2b = 0.0f; sin usar por ahora

		float* F = all_functions + big_index * m;

 		//if (Iexch < 3) {
		if (compute_forces) {
			float3* this_dd = dd + big_index * atoms_n;
			float3* this_Fg = &Fg[index_from3d(dim3(atoms_n, grid_n, m), dim3(atom_i, point_atom_i, 0))];
			float3* this_w3 = &w3[index_from3d(dim3(atoms_n, grid_n, atoms_n), dim3(atom_i, point_atom_i, 0))];
			density_deriv_kernel(dens, num_funcs, nuc_contractions, abs_point_position, 
													 normalize, factor_ac, rmm, nco, big_index, F, Ndens, this_dd, this_Fg, this_w3, atoms_n);
		}
		else {
			local_density_kernel(dens, num_funcs, nuc_contractions, abs_point_position,
													 normalize, factor_ac, rmm, nco, big_index, F, Ndens);
		}
		
		local_pot_kernel(dens, exc_curr, corr_curr, y2a, big_index);

		//printf("atomo: %i layer: %i punto: %i dens: %.12e\n", atom_i, layer_atom_i, point_atom_i, dens);

		/* Numerical Integration */
		float P_total = 0.0f;
		float P_atom_i = 1.0f;

		for (uint atomo_j = 0; atomo_j < atoms_n; atomo_j++) {
			float P_curr = 1.0f;

			float3 pos_atomo_j = gpu_atom_positions[atomo_j];
			float r_atomo_j = distance(abs_point_position,pos_atomo_j);
			float rm_atomo_j = rm_factor[gpu_types[atomo_j]];

			for (uint atomo_k = 0; atomo_k < atoms_n; atomo_k++) {
				if (atomo_k == atomo_j) continue;
				float3 pos_atomo_k = gpu_atom_positions[atomo_k];
				float rr = distance(pos_atomo_j, pos_atomo_k);
				float u = r_atomo_j - distance(abs_point_position, pos_atomo_k);
				u /= rr;

				float x;
				x = rm_atomo_j / rm_factor[gpu_types[atomo_k]];
				x = (x - 1.0f) / (x + 1.0f);
				u += (x / (x * x - 1.0f)) * (1.0f - u * u);

				u = 1.5f * u - 0.5f * (u * u * u);
				u = 1.5f * u - 0.5f * (u * u * u);
				u = 1.5f * u - 0.5f * (u * u * u);
				u = 0.5f * (1.0f - u);
				
				P_curr *= u;
			}

			if (atomo_j == atom_i) P_atom_i = P_curr;
			P_total += P_curr;
		}

		float atom_weight = (P_atom_i / P_total);
		float combined_weight = atom_weight * integration_weight;
		
		// store either the resulting energy or the factor needed to update RMM later
		if (update_rmm || compute_forces) {
			output_factor[big_index] = combined_weight * y2a;
			//_EMU(printf("factor %i %.12e\n", index_from3d(energySize, pos), output_factor[index_from3d(energySize, pos)]));
		}
		
		if (compute_energy)
		{
			float energy_curr = exc_curr + corr_curr;			
			float result = (dens * combined_weight) * energy_curr;
			energy[big_index] = result;
			//_EMU(printf("aca: %i %i %i %.12e %.12e %i\n", atom_i, layer_atom_i, point_atom_i, energy[big_index], result, big_index));
		}
	}
}
