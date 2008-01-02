/* -*- mode: c -*- */
#include <cstdio>
#include "cuda_extra.h"
#include "../matrix.h"
#include "accum.h"
#include "exchnum.h"
#include "exchnum_constants.h"
#include "../timer.h"

#include <cassert>

using namespace G2G;
using namespace std;

/**
 * TODO: revisar distance / distance2 cuando sea necesario
 */

template <const uint* const curr_layers, uint grid_n>
__global__ void calc_new_rmm(const float3* atom_positions, const uint* types, const float3* point_positions,
														 const float* wang, const uint atoms_n, uint nco, uint3 num_funcs,
														 const uint* nuc, const uint* contractions, bool normalize, const float* factor_a, const float* factor_c,
														 const float* rmm, float* rmm_output, float* factors);


__device__ void calc_function(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
															const float3* atom_positions, const float* factor_a, const float* factor_c, uint big_func_index,
															bool normalize, float& func_value);

/***************************************** ENERGY KERNEL ******************************************/
#include "energy.h"

/*************************************** ACTUALIZACION DE RMM *************************************/
#include "rmm.h"


/**
 * Fortran interface
 */

/** TODO: ----- ARREGLAR ESTO ---------- */
#define FORTRAN_MAX_ATOMS 1860
#define FORTRAN_NG 1200
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
	
	#ifdef ENERGY3D
	dim3 threads(natom, MAX_LAYERS, points);
	dim3 blockSize(1, 8, 16);
	dim3 gridSize3d = divUp(threads, blockSize);
	#else
	dim3 threads(natom, points);
	dim3 blockSize(1, 128);
	dim3 gridSize3d = divUp(threads, blockSize);
	#endif
	
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
	assert(Iexch == 1);	
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
	
	HostMatrixFloat energy(1);
	calc_energy(atom_positions, types, igrid, point_positions, energy, wang,
							Ndens, nco, num_funcs_div, nuc, contractions, norm, factor_a, factor_c, rmm, &RMM[m5-1],
							is_int3lu, threads, blockSize, gridSize3d);

	if (!is_int3lu) {
		/* update fortran variables */
		Exc = energy.data[0];
		printf("Exc: %f\n", energy.data[0]);
	}
}

/**
 * Host <-> CUDA Communication function
 */

void calc_energy(const HostMatrixFloat3& atom_positions, const HostMatrixUInt& types, uint grid_type,
								 const HostMatrixFloat3& point_positions, HostMatrixFloat& energy, const HostMatrixFloat& wang,
								 uint Ndens, uint nco, uint3 num_funcs, const HostMatrixUInt& nuc,
								 const HostMatrixUInt& contractions, bool normalize, const HostMatrixFloat& factor_a, const HostMatrixFloat& factor_c,
								 const HostMatrixFloat& rmm, double* cpu_rmm_output, bool update_rmm, const dim3& threads, const dim3& blockSize, const dim3& gridSize3d)
{
	const CudaMatrixFloat3 gpu_atom_positions(atom_positions);
	const CudaMatrixUInt gpu_types(types), gpu_nuc(nuc), gpu_contractions(contractions);
	
	#ifdef ENERGY3D
	uint gridSizeZ = gridSize3d.z;
	dim3 gridSize = dim3(gridSize3d.x, gridSize3d.y * gridSize3d.z, 1);
	#else
	uint gridSizeZ = 1;
	dim3 gridSize = gridSize3d;
	#endif
	
	uint m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;	
	uint small_m = sum(num_funcs);
	
	uint natoms = atom_positions.width;
	uint npoints = point_positions.width;
	
	CudaMatrixFloat gpu_energy, gpu_total_energy, gpu_wang(wang), gpu_factor_a(factor_a), gpu_factor_c(factor_c),
									gpu_rmm(rmm), gpu_functions(m *  natoms * MAX_LAYERS * npoints);
	CudaMatrixFloat3 gpu_point_positions(point_positions);
	
	printf("creando espacio para funcs output: size: %i (%i bytes) data: %i\n", gpu_functions.elements(), gpu_functions.bytes(), (bool)gpu_functions.data);	

	// optional update of RMM(M5)
	CudaMatrixFloat gpu_rmm_output;
	if (update_rmm) {
		// gpu_rmm_output.resize((m * (m + 1)) / 2);
		gpu_rmm_output.resize(m * m); // TODO: ajustar
		printf("creando espacio para rmm output: size: %i (%i bytes) data: %i\n", gpu_rmm_output.elements(), gpu_rmm_output.bytes(), (bool)gpu_rmm_output.data);
	}
	#ifndef _DEBUG
	else
	#endif
		gpu_energy.resize(natoms * MAX_LAYERS * npoints);
	
	printf("threads: %i %i %i, blockSize: %i %i %i, gridSize: %i %i %i\n", threads.x, threads.y, threads.z, blockSize.x, blockSize.y, blockSize.z, gridSize.x, gridSize.y / gridSizeZ, gridSizeZ);
	if (update_rmm) printf("GPU RMM SIZE: %i (%i bytes)\n", gpu_rmm_output.elements(), gpu_rmm_output.bytes());
	printf("energy data elements: %i data: %i\n", gpu_energy.elements(), (bool)gpu_energy.data);
	// TODO: update_rmm should be a template parameter
	const uint* curr_cpu_layers = NULL;
	
	float* factor_output = NULL;
	
	CudaMatrixFloat gpu_factor_output(natoms * MAX_LAYERS * npoints);
	factor_output = gpu_factor_output.data;

	#if 0
	dim3 rmmThreads(divUp(m, 2), m + 1);
	dim3 rmmBlockSize(8,16);
	dim3 rmmGridSize = divUp(rmmThreads, rmmBlockSize);
	#else
	/*dim3 rmmThreads(m,m);
	dim3 rmmBlockSize(1);
	dim3 rmmGridSize(m,m);*/
	dim3 rmmThreads(m, m);
	dim3 rmmBlockSize(8,16);
	dim3 rmmGridSize = divUp(rmmThreads, rmmBlockSize);
	#endif
	
	printf("rmm threads: %i %i, blockSize: %i %i, gridSize: %i %i\n", rmmThreads.x, rmmThreads.y, rmmBlockSize.x, rmmBlockSize.y, rmmGridSize.x, rmmGridSize.y);		
	
	switch(grid_type) {
		case 0:
		{
			energy_kernel<EXCHNUM_SMALL_GRID_SIZE, layers2><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																								 gpu_wang.data, gpu_atom_positions.width, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																								 normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_functions.data,
																																								 Ndens, factor_output, update_rmm);
			
			if (update_rmm) {				
				cudaError_t error = cudaGetLastError();
				if (error != cudaSuccess) fprintf(stderr, "=!=!=!=!=====> CUDA ERROR <=====!=!=!=!=: %s\n", cudaGetErrorString(error));
				
				calc_new_rmm<layers2, EXCHNUM_SMALL_GRID_SIZE><<<rmmGridSize, rmmBlockSize>>>(gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_wang.data,
																																											gpu_atom_positions.width, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																											normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_rmm_output.data,
																																											factor_output, gpu_functions.data);
			}

			curr_cpu_layers = cpu_layers2;
		}
		break;
		case 1:
		{
			energy_kernel<EXCHNUM_MEDIUM_GRID_SIZE, layers><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																								 gpu_wang.data, gpu_atom_positions.width, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																								 normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_functions.data, 
																																								 Ndens, factor_output, update_rmm);
			if (update_rmm) {
				cudaError_t error = cudaGetLastError();
				if (error != cudaSuccess) fprintf(stderr, "=!=!=!=!=====> CUDA ERROR <=====!=!=!=!=: %s\n", cudaGetErrorString(error));

				calc_new_rmm<layers, EXCHNUM_MEDIUM_GRID_SIZE><<<rmmGridSize, rmmBlockSize>>>(gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_wang.data,
																																											 gpu_atom_positions.width, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																											 normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_rmm_output.data,
																																											 factor_output, gpu_functions.data);
			}

			curr_cpu_layers = cpu_layers;
		}
		break;
		case 2:
		{
			energy_kernel<EXCHNUM_BIG_GRID_SIZE, layers><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																							gpu_wang.data, gpu_atom_positions.width, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																							normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_functions.data,
																																							Ndens, factor_output, update_rmm);
			if (update_rmm) {
				cudaError_t error = cudaGetLastError();
				if (error != cudaSuccess) fprintf(stderr, "=!=!=!=!=====> CUDA ERROR <=====!=!=!=!=: %s\n", cudaGetErrorString(error));
				
				calc_new_rmm<layers, EXCHNUM_BIG_GRID_SIZE><<<rmmGridSize, rmmBlockSize>>>(gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_wang.data,
																																										gpu_atom_positions.width, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																										normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, gpu_rmm_output.data,
																																										factor_output, gpu_functions.data);
			}
			curr_cpu_layers = cpu_layers;
		}
		break;
	}
	

	/** CPU Accumulation */
	#ifndef _DEBUG
	if (!update_rmm)
	#endif
	{
		double energy_double = 0.0;		
		energy = gpu_energy;
		
		for (unsigned int i = 0; i < natoms; i++) {
			for (unsigned int j = 0; j < curr_cpu_layers[types.data[i]]; j++) {
				for (unsigned int k = 0; k < npoints; k++) {
					#ifdef ENERGY3D
					uint idx = index_from3d(threads, dim3(i, j, k));
					#else
					uint idx = index_from3d(dim3(threads.x, MAX_LAYERS, threads.y), dim3(i, j, k));
					#endif
					printf("idx: %i size: %i\n", idx, energy.elements());

					double energy_curr = energy.data[idx];
					printf("atomo: %i, capa: %i, punto: %i, valor: %.12e idx: %i\n", i, j, k, energy_curr, idx);
					energy_double += energy_curr;
					energy.data[0] += energy_curr;
				}
			}
		}
		printf("Energy (double): %.12e\n", energy_double);		
	}

	/** RMM update **/
	if (update_rmm) {
		HostMatrixFloat gpu_rmm_output_copy(gpu_rmm_output);
		
		uint rmm_idx = 0;
		for (uint func_i = 0; func_i < m; func_i++) {
			for (uint func_j = func_i; func_j < m; func_j++) {
				printf("rmm_output(%i): %.12e\n", rmm_idx, gpu_rmm_output_copy.data[rmm_idx]);
				cpu_rmm_output[rmm_idx] += gpu_rmm_output_copy.data[rmm_idx];
				rmm_idx++;
			}
		}
	}
	
	// calc_accum_cuda(gpu_energy, gpu_total_energy);

	// TODO: esta copia es redundante con la que hay en calc_acuum (esa es GPU<->GPU)
	// energy.copy_submatrix(gpu_total_energy, 1);
	// 
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) fprintf(stderr, "=!=!=!=!=====> CUDA ERROR <=====!=!=!=!=: %s\n", cudaGetErrorString(error));
}


/************************************************** FUNCTIONS ****************************************/
#include "functions.h"

/************************************* DENSITY KERNEL ******************************/

#include "density.h"

/******************************** POT KERNEL ***********************************/

#include "pot.h"
