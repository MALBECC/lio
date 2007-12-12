/* -*- mode: c -*- */
#include <cstdio>
#include "cuda_extra.h"
#include "../matrix.h"
#include "accum.h"
#include "exchnum.h"
#include "exchnum_constants.h"
using namespace G2G;
using namespace std;

#ifdef __DEVICE_EMULATION__
#include <fpu_control.h>
#endif


/**
 * TODO: revisar distance / distance2 cuando sea necesario
 */

template <unsigned int grid_n, const uint* const curr_layers>
	__global__ void energy_kernel(uint gridSizeZ, const float3* atom_positions, const uint* types, const float3* point_positions,
																float* energy, const float* wang, const uint atoms_n, uint Iexch, uint nco, uint3 num_funcs,
																const uint* nuc, const uint* contractions, bool normalize, const float* factor_a, const float* factor_c,
																const float* rmm, float* output_rmm, uint Ndens, bool is_int3lu);

/**
 * Fortran interface
 */

/** TODO: ----- ARREGLAR ESTO ---------- */
#define FORTRAN_MAX_ATOMS 1845
#define FORTRAN_NG 900
#define FORTRAN_NL 10

/**
 * Parametros innecesarios: m (es sum(num_funcs))
 */
extern "C" void exchnum_gpu_(const unsigned int& norm, const unsigned int& natom, const double* r, const unsigned int* Iz, const unsigned int* Nuc,
														 const unsigned int& m, const unsigned int* ncont, const unsigned int* nshell, const double* c, const double* a,
														 double* RMM, const unsigned int& m18, const unsigned int& m5, const unsigned int& nco, double& Exc, const unsigned int& nopt,
														 const unsigned int& Iexch, const unsigned int& igrid,
														 const double* e, const double* e2, const double* e3,
														 const double* fort_wang, const double* fort_wang2, const double* fort_wang3,
														 const unsigned int& Ndens, const unsigned int& is_int3lu)
{
	
/*	#ifdef __DEVICE_EMULATION__
	printf("Enabling 32bit floats\n");
	unsigned int original_control_word;
	_FPU_GETCW(original_control_word);
	unsigned int new_control_word = (original_control_word & ~0x300) | 0x000;
	_FPU_SETCW(new_control_word);
	#endif*/
	
	printf("<======= exchnum_gpu (from %s) ============>\n", is_int3lu ? "int3lu" : "SCF");
	printf("Ndens: %i\n", Ndens);
	uint3 num_funcs = make_uint3(nshell[0], nshell[1], nshell[2]);
	uint3 num_funcs_div = num_funcs / make_uint3(1, 3, 6);
	
	uint total_funcs = sum(num_funcs);
	uint total_funcs_div = sum(num_funcs_div);
	
	uint points = EXCHNUM_SMALL_GRID_SIZE;
	switch (igrid) {
		case 0: points = EXCHNUM_SMALL_GRID_SIZE; 	break;
		case 1: points = EXCHNUM_MEDIUM_GRID_SIZE;	break;
		case 2: points = EXCHNUM_BIG_GRID_SIZE;			break;
	}
	
	dim3 threads(natom, MAX_LAYERS, points);
	dim3 blockSize(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 gridSize3d = divUp(threads, blockSize);	
	
	HostMatrixFloat3 atom_positions(natom), point_positions(points);
	/* output_rmm size: TODO: divUp(m * (m - 1),2) */
	HostMatrixFloat factor_a(total_funcs, MAX_CONTRACTIONS), factor_c(total_funcs, MAX_CONTRACTIONS), wang(points);
	HostMatrixUInt types(natom), nuc(total_funcs_div), contractions(total_funcs_div);
	
	// REVISAR: nuc: imagen y dominio (especialmente por la parte de * 3 y * 6)

	printf("%i atoms\n", natom);
	for (unsigned int i = 0; i < natom; i++) {
		atom_positions.data[i] = make_float3(r[FORTRAN_MAX_ATOMS * 0 + i], r[i + FORTRAN_MAX_ATOMS * 1], r[i + FORTRAN_MAX_ATOMS * 2]);
		//printf("Pos(%i): %f %f %f, Types(%i): %i\n", i, atom_positions.data[i].x, atom_positions.data[i].y, atom_positions.data[i].z, i, Iz[i]);		
		types.data[i] = Iz[i] - 1;
	}
	
	printf("ns: %i, np: %i, nd: %i, Total_Funcs: %i\n", num_funcs.x, num_funcs.y, num_funcs.z, total_funcs);
	{
		uint inc = 1;
		uint i, j;
		for (i = 0, j = 0; i < total_funcs; i += inc, j++) {
			if (i == num_funcs.x) inc = 3;
			else if (i == num_funcs.x + num_funcs.y) inc = 6;

			//printf("i: %i, j: %i\n", i, j);
			//printf("Nuc(%i) = %i\n", i, Nuc[i] - 1);
			//printf("ncont(%i) = %i\n", i, ncont[i]);
			nuc.data[j] = Nuc[i] - 1;
			contractions.data[j] = ncont[i];
			
			for (unsigned int k = 0; k < ncont[i]; k++) {
				factor_a.data[j * MAX_CONTRACTIONS + k] = a[FORTRAN_NG * k + i];
				factor_c.data[j * MAX_CONTRACTIONS + k] = c[FORTRAN_NG * k + i];
				//printf("cont: %i, a: %f, c: %f\n", k, factor_a.data[j * MAX_CONTRACTIONS + k], factor_c.data[j * MAX_CONTRACTIONS + k]);
			}			
		}
	}
	
	HostMatrixFloat rmm;
	printf("NCO: %i, M: %i, Iexch: %i\n", nco, total_funcs, Iexch);
	{
		if (Ndens == 1) {
			rmm.resize(m * m);
			uint k = 0;
			for (unsigned int i = 0; i < m; i++) {
				for (unsigned int j = i; j < m; j++) {
					rmm.data[k] = RMM[k];
					//printf("rmm(%i): %.30e\n", k, RMM[m5 + k - 1]);
					k++;
				}
			}
		}
		else {
			rmm.resize(m, nco);
			uint k = m18 - 1;
			for (unsigned int i = 0; i < m; i++) {
				for (unsigned int j = 0; j < nco; j++) {
					rmm.data[i * nco + j] = RMM[k];
					//printf("rmm(%i,%i): %.30e (%i)\n", i, j, RMM[k], k);
					k++;
				}
			}
		}
	}

	const double* real_e = NULL;
	const double* real_wang = NULL;
	switch (igrid) {
		case 0: real_e = e;  real_wang = fort_wang;  	break;
		case 1: real_e = e2; real_wang = fort_wang2; 	break;
		case 2: real_e = e3; real_wang = fort_wang3; 	break;		
	}

	printf("Puntos (grilla %i):\n", igrid);	
	for (unsigned int i = 0; i < points; i++) {
		wang.data[i] = real_wang[i];
		point_positions.data[i] = make_float3(real_e[0 * points + i], real_e[1 * points + i], real_e[2 * points + i]);
		//printf("wang: %f, e: (%f,%f,%f)\n", wang.data[i], point_positions.data[i].x, point_positions.data[i].y, point_positions.data[i].z);
	}
	
	HostMatrixDouble rmm_partial_out(m * m);
	rmm_partial_out.fill(0.0f);
		
	HostMatrixFloat energy(1);
	calc_energy(atom_positions, types, igrid, point_positions, energy, wang,
							Iexch, Ndens, nco, num_funcs_div, nuc, contractions, norm, factor_a, factor_c, rmm, &RMM[m5-1],
							is_int3lu, threads, blockSize, gridSize3d);

	/* update fortran variables */
	Exc = energy.data[0];
	printf("Exc: %f\n", energy.data[0]);
	
/*	#ifdef __DEVICE_EMULATION__
	printf("Restoring normal floats\n");
	_FPU_SETCW(original_control_word);
  #endif*/
}

/**
 * Host <-> CUDA Communication function
 */

void calc_energy(const HostMatrixFloat3& atom_positions, const HostMatrixUInt& types, uint grid_type,
								 const HostMatrixFloat3& point_positions, HostMatrixFloat& energy, const HostMatrixFloat& wang,
								 uint Iexch, uint Ndens, uint nco, uint3 num_funcs, const HostMatrixUInt& nuc,
								 const HostMatrixUInt& contractions, bool normalize, const HostMatrixFloat& factor_a, const HostMatrixFloat& factor_c,
								 const HostMatrixFloat& rmm, double* cpu_rmm_output, bool is_int3lu, const dim3& threads, const dim3& blockSize, const dim3& gridSize3d)
{	
	const CudaMatrixFloat3 gpu_atom_positions(atom_positions);
	const CudaMatrixUInt gpu_types(types), gpu_nuc(nuc), gpu_contractions(contractions);
	
	uint gridSizeZ = gridSize3d.z;
	dim3 gridSize = dim3(gridSize3d.x, gridSize3d.y * gridSize3d.z, 1);

	CudaMatrixFloat gpu_energy(threads.x * threads.y * threads.z), gpu_total_energy, gpu_wang(wang), gpu_factor_a(factor_a), gpu_factor_c(factor_c),
									gpu_rmm(rmm);
	CudaMatrixFloat3 gpu_point_positions(point_positions);
	
	// optional update of RMM(M5)
	CudaMatrixFloat gpu_rmm_output;
	if (is_int3lu) gpu_rmm_output.resize(threads.x * threads.y * threads.z * MAX_FUNCTIONS * MAX_FUNCTIONS);
	
	printf("threads: %i %i %i, blockSize: %i %i %i, gridSize: %i %i %i\n", threads.x, threads.y, threads.z, blockSize.x, blockSize.y, blockSize.z, gridSize.x, gridSize.y / gridSizeZ, gridSizeZ);
	if (is_int3lu) printf("GPU RMM SIZE: %i (%i bytes)\n", gpu_rmm_output.elements(), gpu_rmm_output.bytes());
	// TODO: is_int3lu should be a template parameter
	const uint* curr_cpu_layers = NULL;

	switch(grid_type) {
		case 0:
		{
			energy_kernel<EXCHNUM_SMALL_GRID_SIZE, layers2><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																								 gpu_wang.data, gpu_atom_positions.width, Iexch, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																								 normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_rmm_output.data, Ndens, is_int3lu);
			curr_cpu_layers = cpu_layers2;
		}
		break;
		case 1:
		{
			energy_kernel<EXCHNUM_MEDIUM_GRID_SIZE, layers><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																								 gpu_wang.data, gpu_atom_positions.width, Iexch, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																								 normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_rmm_output.data, Ndens, is_int3lu);
			curr_cpu_layers = cpu_layers;
		}
		break;
		case 2:
		{
			energy_kernel<EXCHNUM_BIG_GRID_SIZE, layers><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																							gpu_wang.data, gpu_atom_positions.width, Iexch, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																							normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_rmm_output.data, Ndens, is_int3lu);
			curr_cpu_layers = cpu_layers;
		}
		break;
	}
	

	/** CPU Accumulation */
	energy = gpu_energy;
	
	HostMatrixFloat gpu_rmm_output_copy(gpu_rmm_output);
	uint m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;
	
	double energy_double = 0.0;
	for (unsigned int i = 0; i < threads.x; i++) {
		//printf("atomo: %i, capas: %i\n", i, curr_cpu_layers[i]);
		for (unsigned int j = 0; j < curr_cpu_layers[types.data[i]]; j++) {
			for (unsigned int k = 0; k < threads.z; k++) {
				uint idx = index_from3d(threads, dim3(i, j, k));
				double energy_curr = energy.data[idx];
				printf("atomo: %i, capa: %i, punto: %i, valor: %.12e idx: %i\n", i, j, k, energy_curr, idx);
				energy_double += energy_curr;
				energy.data[0] += energy_curr;
				
				uint big_rmm_idx = index_from4d(threads + make_uint4(0,0,0,MAX_FUNCTIONS * MAX_FUNCTIONS), make_uint4(i, j, k, 0));

				if (is_int3lu) {
					uint rmm_idx = 0;
					for (uint func_i = 0; func_i < m; func_i++) {
						for (uint func_j = func_i; func_j < m; func_j++) {
							if (idx == 996) printf("idx: %i rmm_out(%i): %.12e %.12e %.12e\n", idx, rmm_idx, gpu_rmm_output_copy.data[big_rmm_idx + rmm_idx], cpu_rmm_output[rmm_idx],
																		 cpu_rmm_output[rmm_idx] + gpu_rmm_output_copy.data[big_rmm_idx + rmm_idx]);
							cpu_rmm_output[rmm_idx] += gpu_rmm_output_copy.data[big_rmm_idx + rmm_idx];
							
							rmm_idx++;
						}
					}
				}				
			}
		}
	}
	
	//for (unsigned int i = 1; i < energy.elements(); i++) { energy_double += energy.data[i]; energy.data[0] += energy.data[i]; }
	printf("Energy (double): %.12e\n", energy_double);
	
	// calc_accum_cuda(gpu_energy, gpu_total_energy);

	// TODO: esta copia es redundante con la que hay en calc_acuum (esa es GPU<->GPU)
	// energy.copy_submatrix(gpu_total_energy, 1);
	// 
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) printf("CUDA ERROR: %s\n", cudaGetErrorString(error));
}

/**
 * Main Energy Kernel
 */

/*, bool Ndens, unsigned int Iexch*/
template <unsigned int grid_n, const uint* const curr_layers>
	__global__ void energy_kernel(uint gridSizeZ, const float3* atom_positions, const uint* types, const float3* point_positions,
																float* energy, const float* wang, const uint atoms_n, uint Iexch, uint nco, uint3 num_funcs,
																const uint* nuc, const uint* contractions, bool normalize, const float* factor_a, const float* factor_c,
																const float* rmm, float* rmm_output, uint Ndens, bool update_rmm)
{
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);	
	dim3 pos = index(blockDim, dim3(blockIdx.x, blockIdx.y / gridSizeZ, blockIdx.y % gridSizeZ), threadIdx);
	uint big_index = index_from3d(energySize, pos);
	const uint& m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;	

	const uint& atom_i = pos.x;	
 	const uint& layer_atom_i = pos.y;
	const uint& point_atom_i = pos.z;
	
	// Hay varios lugares donde se podrian compartir entre threads, solo calculando una vez
		
	/* load atom and point positions */
	/*__shared__ float3 atom_positions_local[MAX_ATOMS];
	__shared__ float3 point_positions_local[grid_n];
	
	if (layer_atom_i == 0 && point_atom_i == 0) atom_positions_local[atom_i] = atom_positions[atom_i];
	if (atom_i == 0 && layer_atom_i == 0) point_positions_local[point_atom_i] = point_positions[point_atom_i];*/

	// __syncthreads();
			
	// atom_positions = atom_positions_local;
	// point_positions = point_positions_local;

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
	
	float F[MAX_FUNCTIONS];
	
	// float exc_curr, corr_current;
	density_kernel(dens, num_funcs, nuc, contractions, point_position, atom_positions, normalize, factor_a, factor_c, rmm, nco, big_index, F, Ndens);
	pot_kernel(dens, exc_curr, corr_curr, y2a, Iexch, big_index);
	
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
			float u = r_atomo_j - distance(point_position, atom_positions[atomo_k]); // revisar que se haga igual que abajo
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
	
	// store result
	float energy_curr = exc_curr + corr_curr;
	tmp0 = (P_atom_i / P_total) * integration_weight;
	float result = (dens * tmp0) * energy_curr;

	energy[index_from3d(energySize, pos)] = result;

	if (update_rmm) {
		if (dens == 0.0f) {
			tmp0 = 0.0f;
			energy_curr = 0.0f;
			y2a = 0.0f;
			// y2b = 0.0f;
		}
		
		float tmp1 = tmp0 * y2a;
		
		uint big_k = index_from4d(energySize + make_uint4(0, 0, 0, MAX_FUNCTIONS * MAX_FUNCTIONS), pos + make_uint4(0,0,0,0));
		uint k = 0;
		//if (tmp1 != 0.0f) {
			for (uint i = 0; i < m; i++) {
				float tmp2 = tmp1 * F[i];
				
				/*if (tmp2 == 0.0f) {
					for (uint j = i; j < m; j++) {
						rmm_output[k] = 0.0f;
						k++;
					}
				}
				else {*/
					for (uint j = i; j < m; j++) {
						// rmm_output[k] = F[j] * tmp2;
						if (big_index == 996)
						  _EMU(printf("rmm_out_calculo: %i %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n", k, F[j], tmp2, tmp1, y2a, tmp0, P_atom_i / P_total, integration_weight, F[j] * tmp2));
						rmm_output[big_k + k] = F[j] * tmp2;
						k++;
					}
				//}
			}
		//}
 	}	
  //_EMU(printf("idx: %i dens: %.12e e: %.12e r: %.12e\n", big_index,  dens, energy_curr, result));
}

__device__ void calc_function_s(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float* factor_a, const float* factor_c, uint big_index, uint func_index, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	*func_value = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;
		*func_value += expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];
	}	
}

__device__ void calc_function_p(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float* factor_a, const float* factor_c, uint big_index, uint func_index, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	for (uint i = 0; i < 3; i++) func_value[i] = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];
		float3 v = (point_position - atom_positions[atom_nuc]) * t;
		
		for (uint i = 0; i < 3; i++) func_value[i] += elem(v,i);
	}
}

__device__ void calc_function_d(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float* factor_a, const float* factor_c, uint big_index, uint func_index,
																float normalization_factor, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	for (uint i = 0; i < 6; i++) func_value[i] = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];

		float3 v = point_position - atom_positions[atom_nuc];

		uint ij = 0;
		for (uint i = 0; i < 3; i++) {
			float t1 = elem(v, i);
			for (uint j = 0; j <= i; j++) {
				float t2 = (i == j ? elem(v, j) * normalization_factor : elem(v, j));
				func_value[ij] += t * t1 * t2;
				ij++;
			}
		}
	}
}


/* Local Density Functionals */
__device__ void density_kernel(float& density, uint3 num_funcs, const uint* nuc, const uint* contractions, float3 point_position,
															 const float3* atom_positions, bool normalize, const float* factor_a,
															 const float* factor_c, const float* rmm, uint nco, uint big_index, float* F, uint Ndens)
{
	/*
	 * now we should evaluate all same loops as the ones used for
	 * 1 electron matrix elements, but doing only products
	 * then, the particular density functional wanted is calculated
	 */

	density = 0.0f;

	const uint& funcs_s = num_funcs.x;
	const uint& funcs_p = num_funcs.y;
	uint funcs_total = sum(num_funcs);
	const uint& m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;
	uint func_real = 0;
	
	/* s functions */
	for (uint func = 0; func < funcs_s; func++, func_real++)
		calc_function_s(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, big_index, func, &F[func_real]);
	
	/* p functions */
	for (uint func = funcs_s; func <  funcs_s + funcs_p; func++, func_real+=3)
		calc_function_p(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, big_index, func, &F[func_real]);
	
	/* d functions */
	float normalization_factor = (normalize ? rsqrtf(3.0f) : 1.0f);
	
	for (uint func = (funcs_s + funcs_p); func < funcs_total; func++, func_real+=6)
		calc_function_d(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, big_index, func, normalization_factor, &F[func_real]);

	/* density */
	if (Ndens == 1) {
		uint k = 0;
		
		for (uint i = 0; i < m; i++) {
			if (F[i] == 0.0f) {
				k += m - i;
				continue;
			}
			
			float w = 0.0f;			
			for (uint j = i; j < m; j++) {
				w += rmm[k] * F[j];
				k++;				
			}
			
			density += F[i] * w;
		}
	}
	else {
		uint k = 0;
		for (uint i = 0; i < nco; i++) {
			float w = 0.0f;
			for (uint j = 0; j < m; j++) {
				w += rmm[k] * F[j];
				k++;
			}
			density += w * w;
		}
		density *= 2.0f;		
	}
}

/*
 * function for evaluating exchange correlation density and
 * potential, depends of Iexch, index that says which potential
 * will be used
 * - 1 X alpha
 * - 2 Gunnarson-Lundqvist
 * - 3 Vosko, Wilk and Nusair
 * Exchange part is in all cases the same, for the time being LD
 * Self interaction corrections are used in correlation part
 */
	
/*template<unsigned int Iexch>*/ __device__ void pot_kernel(float dens, float& ex, float& ec, float& y2a, uint Iexch, uint big_index)
{
	// data X alpha

	if (dens == 0) { ex = 0.0f; ec = 0.0f; y2a = 0.0f; return; }
	
	float y = powf(dens, 0.333333333333333333f);
	float e0 = pot_alpha * y;	
	float v0 = -0.984745021842697f * y; // 4/3 * pot_alpha * y
	
	switch(Iexch) {
		case 1:
		{
			ex = e0;
			ec = 0;
			y2a = v0;
		}
		break;
		case 2:
		{
			ex = e0;
			float rs = pot_gl / y;
			float x1 = rs / 11.4f;
			float vc;
			
			if (x1 > 1.0f) {
				ec = -0.0333f * (0.5f * x1 - 0.33333333333333f);
				vc = 0.0111f * x1 * 0.5f;
			}
			else {
				float t1 = (1.0f + x1 * x1 * x1);
				float t2 = logf(1.0f + 1.0f / x1);
				float t3 = x1 * x1;
        ec = -0.0333f * (t1 * t2 - t3 + 0.5f * x1 - 0.33333333333333f);
        vc = 0.0111f * x1 * (3.0f * t3 * t2 - t1 / (x1 * (x1 + 1.0f)) - 2.0f * x1 + 0.5f);
			}
			y2a = v0 + ec + vc;
		}
		break;
		case 3:
		{
			ex = e0;
			float rs = pot_gl / y;
			float x1 = sqrtf(rs);
			float Xx = rs + pot_vosko_b1 * x1 + pot_vosko_c1;
			float Xxo = pot_vosko_x0 * pot_vosko_x0 + pot_vosko_b1 * pot_vosko_x0 + pot_vosko_c1;
			float t1 = 2 * x1 + pot_vosko_b1;
			float t2 = logf(Xx);
			float t3 = atanf(pot_vosko_q/t1);
			float t4 = pot_vosko_b1 * pot_vosko_x0 / Xxo;
      float t5 = (pot_vosko_b1 * x1 + 2.0f * pot_vosko_c1) / x1;
			float t6 = pot_vosko_x0 / Xxo;			
			
      ec = pot_vosko_a1 * (2 * logf(x1) - t2 + 2 * pot_vosko_b1 / pot_vosko_q * t3 - t4 *
								 (2 * logf(x1 - pot_vosko_x0) - t2 + 2 * (pot_vosko_b1 + 2 * pot_vosko_x0) / pot_vosko_q * t3));
			
      float vc = ec - pot_vosko_a16 * x1 * (t5 / Xx - 4.0f * pot_vosko_b1 / (t1 * t1 + pot_vosko_q * pot_vosko_q) * (1.0f - t6 * (pot_vosko_b1 - 2.0f * pot_vosko_x0)) -
																						t4 * (2.0f / (x1 - pot_vosko_x0) - t1 / Xx));
			
			y2a = v0 + vc;
		}
		break;		
	}
}
