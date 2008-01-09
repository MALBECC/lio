/**
 * Main Energy Kernel
 */

__device__ float3 get_atom_position(uint atom_number, const float3* atom_positions_shared, const float3* atom_positions_global)
{
	return (atom_number < ENERGY_SHARED_ATOM_POSITIONS ? atom_positions_shared[atom_number] : atom_positions_global[atom_number]);
}


template <uint grid_type/*, unsigned int grid_n, const uint* const curr_layers*/>
__global__ void energy_kernel(uint gridSizeZ, const float3* atom_positions, const uint* types,
		float* energy, const uint atoms_n, uint nco, uint3 num_funcs,
		const uint* nuc, const uint* contractions, bool normalize, const float2* factor_ac, /*const float* factor_c,*/
		const float* rmm, float* all_functions,  uint Ndens, float* output_factor, bool update_rmm)
{
	/** TODO: estos ifs se hacen porque no puedo pasar estos punteros por parametro ni por template
	 * (esto ultimo porque el compilador de NVIDIA muere) **/
	uint grid_n = 0;
	const uint* curr_layers = NULL;
	const float3* point_positions = NULL;
	const float* wang = NULL;
	
	if (grid_type == 0) {
		grid_n = EXCHNUM_SMALL_GRID_SIZE;
		curr_layers = layers2;
		point_positions = small_grid_positions;
		wang = small_wang;
	}
	else if (grid_type == 1) {
		grid_n = EXCHNUM_MEDIUM_GRID_SIZE;
		curr_layers = layers;
		point_positions = medium_grid_positions;
		wang = medium_wang;		
	}
	else if (grid_type == 2) {
		grid_n = EXCHNUM_BIG_GRID_SIZE;
		curr_layers = layers;
		point_positions = big_grid_positions;
		wang = big_wang;		
	}
	
	
	const uint& m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;				 
	
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);
	dim3 pos2d = index(blockDim, blockIdx, threadIdx);
	
	__shared__ float3 atom_positions_shared[ENERGY_SHARED_ATOM_POSITIONS];
	
	const uint atom_i = pos2d.x;
	const uint point_atom_i = pos2d.y;
	
	__syncthreads();	

	// load shared data
	if (threadIdx.y == 0) {
		for (uint i = 0; i < min(atoms_n, ENERGY_SHARED_ATOM_POSITIONS); i++) atom_positions_shared[i] = atom_positions[i];
	}
	
	// determine if thread is valid
	bool valid_thread = true;
	if (atom_i >= atoms_n) valid_thread = false;
	if (point_atom_i >= grid_n) valid_thread = false;
	
	// not loaded from shared
	float3 rel_point_position = point_positions[point_atom_i];	// constant memory
	float wang_point_i = wang[point_atom_i]; // constant memory

	// fill local variables from shared data
	float3 atom_i_position = make_float3(0.0f,0.0f,0.0f);
	
	float rm = 0.0f;
	uint atom_i_layers = 0;
	float tmp0 = 0.0f;

	if (valid_thread) {
		uint atom_i_type = types[atom_i];
		rm = rm_factor[atom_i_type]; // constant memory
		atom_i_layers = curr_layers[atom_i_type];	// constant memory
		tmp0 = (PI / (atom_i_layers + 1.0f));
	}
	
	__syncthreads();
	
	if (!valid_thread) return;
	else {
		atom_i_position = get_atom_position(atom_i, atom_positions_shared, atom_positions);
	}	
	
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
		local_density_kernel(dens, num_funcs, nuc, contractions, abs_point_position, atom_positions, atom_positions_shared,
												 normalize, factor_ac/*, factor_c*/, rmm, nco, big_index, F, Ndens);
		local_pot_kernel(dens, exc_curr, corr_curr, y2a, big_index);
		/*}
		 else {
		 float Dx,Dy,Dz,Dxx,Dxy,Dxz,Dyy,Dzz,Dyz;
		 nonlocal_density_kernel(dens, num_funcs, nuc, contractions, abs_point_position, atom_positions, normalize, factor_a, factor_c, rmm, nco, big_index, F, Ndens,
		 Dx,Dy,Dz,Dxx,Dxy,Dxz,Dyy,Dzz,Dyz);
		 nonlocal_pot_kernel(dens, exc_curr, corr_curr, y2a, big_index,
		 Dx,Dy,Dz,Dxx,Dxy,Dxz,Dyy,Dzz,Dyz);
		 }*/

		//printf("atomo: %i layer: %i punto: %i dens: %.12e\n", atom_i, layer_atom_i, point_atom_i, dens);

		/* Numerical Integration */
		float P_total = 0.0f;
		float P_atom_i = 1.0f;

		for (uint atomo_j = 0; atomo_j < atoms_n; atomo_j++) {
			float P_curr = 1.0f;

			float3 pos_atomo_j = get_atom_position(atomo_j, atom_positions_shared, atom_positions);
			float r_atomo_j = distance(abs_point_position,pos_atomo_j);
			float rm_atomo_j = rm_factor[types[atomo_j]];

			for (uint atomo_k = 0; atomo_k < atoms_n; atomo_k++) {
				if (atomo_k == atomo_j) continue;
				float3 pos_atomo_k = get_atom_position(atomo_k, atom_positions_shared, atom_positions);
				float rr = distance(pos_atomo_j, pos_atomo_k);
				float u = r_atomo_j - distance(abs_point_position, pos_atomo_k);
				u /= rr;

				float x = rm_atomo_j / rm_factor[types[atomo_k]];
				float x1 = (x - 1.0f) / (x + 1.0f);
				float Aij = x1 / (x1 * x1 - 1.0f);
				u += Aij * (1.0f - u * u);

				float p1 = 1.5f * u - 0.5f * (u * u * u);
				float p2 = 1.5f * p1 - 0.5f * (p1 * p1 * p1);
				float p3 = 1.5f * p2 - 0.5f * (p2 * p2 * p2);
				float s = 0.5f * (1.0f - p3);

				P_curr *= s;
			}

			if (atomo_j == atom_i) P_atom_i = P_curr;
			P_total += P_curr;
		}

		float atom_weight = (P_atom_i / P_total);
		float combined_weight = atom_weight * integration_weight;

		// store either the resulting energy or the factor needed to update RMM later
		if (update_rmm) {
			output_factor[big_index] = combined_weight * y2a;
			//_EMU(printf("factor %i %.12e\n", index_from3d(energySize, pos), output_factor[index_from3d(energySize, pos)]));
		}
#ifndef _DEBUG
		else
#endif
		{
			float energy_curr = exc_curr + corr_curr;			
			float result = (dens * combined_weight) * energy_curr;
			energy[big_index] = result;
			_EMU(printf("aca: %i %i %i %.12e %.12e %i\n", atom_i, layer_atom_i, point_atom_i, energy[big_index], result, big_index));
		}

	}
}
