/**
 * Main Energy Kernel
 */

template <unsigned int grid_n, const uint* const curr_layers>
__global__ void energy_kernel(uint gridSizeZ, const float3* atom_positions, const uint* types, const float3* point_positions,
		float* energy, const float* wang, const uint atoms_n, uint nco, uint3 num_funcs,
		const uint* nuc, const uint* contractions, bool normalize, const float* factor_a, const float* factor_c,
		const float* rmm, float* all_functions,  uint Ndens, float* output_factor, bool update_rmm)
{
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);
	dim3 pos = index(blockDim, dim3(blockIdx.x, blockIdx.y / gridSizeZ, blockIdx.y % gridSizeZ), threadIdx);
	uint big_index = index_from3d(energySize, pos);
	const uint& m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;	

	const uint& atom_i = pos.x;	
 	const uint& layer_atom_i = pos.y;
	const uint& point_atom_i = pos.z;
	
	// Hay varios lugares donde se podrian compartir entre threads, solo calculando una vez
		
	/* skip things that shouldn't be computed */
	// CUIDADO al hacer __syncthreads despues de esto
	if (atom_i >= atoms_n) return;
	if (point_atom_i >= grid_n) return;
	
	// Datos por atomo
	// printf("atomo: %i, punto: %i, layer: %i\n", atom_i, point_atom_i, layer_atom_i);
	uint atom_i_type = types[atom_i];
	uint atom_i_layers = curr_layers[atom_i_type];
	
	if (layer_atom_i >= atom_i_layers) return;	
		
	// Datos por capa
	float rm = rm_factor[atom_i_type];
	float tmp0 = (PI / (atom_i_layers + 1.0f));
	float tmp1 = tmp0 * (layer_atom_i + 1.0f);
	float x = cosf(tmp1);
	float r1 = rm * (1.0f + x) / (1.0f - x);
 	float w = tmp0 * fabsf(sinf(tmp1));
	float wrad = w * (r1 * r1) * rm * 2.0f / ((1.0f - x) * (1.0f - x));
	//if (point_atom_i == 0) printf("atom: %i layer: %i rm: %.12e tmp0: %.12e tmp1: %.12e x: %.12e\n", atom_i, layer_atom_i, rm, tmp0, tmp1, x);
	
	// Datos por punto
	float3 point_position = atom_positions[atom_i] + point_positions[point_atom_i] * r1;
	float integration_weight = wang[point_atom_i] * wrad;

	float exc_curr = 1.0f;
	float corr_curr = 1.0f;
	float dens = 0.0f;
	float y2a = 0.0f;
	// float y2b = 0.0f; sin usar por ahora
	
	float* F = all_functions + big_index * m;
	
	// float exc_curr, corr_current;
	density_kernel(dens, num_funcs, nuc, contractions, point_position, atom_positions, normalize, factor_a, factor_c, rmm, nco, big_index, F, Ndens);
	pot_kernel(dens, exc_curr, corr_curr, y2a,  big_index);
	
	//printf("atomo: %i layer: %i punto: %i dens: %.12e\n", atom_i, layer_atom_i, point_atom_i, dens);
	
	/* Numerical Integration */
	float P_total = 0.0f;
	float P_atom_i = 1.0f;
	
	for (uint atomo_j = 0; atomo_j < atoms_n; atomo_j++) {
		float P_curr = 1.0f;
		
		float3 pos_atomo_j = atom_positions[atomo_j];
		float r_atomo_j = distance(point_position,atom_positions[atomo_j]);

		for (uint atomo_k = 0; atomo_k < atoms_n; atomo_k++) {
			if (atomo_k == atomo_j) continue;
			float rr = distance(pos_atomo_j, atom_positions[atomo_k]);
			float u = r_atomo_j - distance(point_position, atom_positions[atomo_k]);
			u /= rr;

			float x = rm_factor[types[atomo_j]] / rm_factor[types[atomo_k]];
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
	
	// store result
	float energy_curr = exc_curr + corr_curr;
	tmp0 = atom_weight * integration_weight;
	float result = (dens * tmp0) * energy_curr;

	// store either the resulting energy or the factor needed to update RMM later
	if (update_rmm) {
		output_factor[index_from3d(energySize, pos)] = tmp0 * y2a;
		//_EMU(printf("factor %i %.12e\n", index_from3d(energySize, pos), output_factor[index_from3d(energySize, pos)]));
 	}	
	#ifndef _DEBUG
	else
	#endif
		energy[big_index] = result;

  //_EMU(printf("idx: %i dens: %.12e e: %.12e r: %.12e\n", big_index,  dens, energy_curr, result));
}
