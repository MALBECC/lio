/* -*- mode: c -*- */
#include <cstdio>
#include "cuda_extra.h"
#include "../matrix.h"
#include "accum.h"
#include "exchnum.h"
#include "exchnum_constants.h"
using namespace G2G;
using namespace std;


/**
 * TODO: revisar distance / distance2 cuando sea necesario
 */

template <unsigned int grid_n>
	__global__ void energy_kernel(uint gridSizeZ, const float3* atom_positions, const uint* types, const float3* point_positions,
																float* energy, const float* wang, const uint atoms_n, uint Iexch, uint nco, uint3 num_funcs,
																const uint* nuc, const uint* contractions, bool normalize, const float* factor_a, const float* factor_c,
																const float* rmm, bool Ndens);

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
														 const double* RMM, const unsigned int& m18, const unsigned int& nco, double& Exc, const unsigned int& nopt,
														 const unsigned int& Iexch, const unsigned int& igrid, const double* e, const double* e3, const double* fort_wang, const double* fort_wang3,
														 const unsigned int& Ndens)
{
	printf("<======= exchnum_gpu ============>\n");
	printf("Ndens: %i\n", Ndens);
	uint3 num_funcs = make_uint3(nshell[0], nshell[1], nshell[2]);
	uint3 num_funcs_div = num_funcs / make_uint3(1, 3, 6);
	
	uint total_funcs = sum(num_funcs);
	uint total_funcs_div = sum(num_funcs_div);
	
	uint points = (igrid == 1 ? EXCHNUM_GRIDA_SIZE : EXCHNUM_GRIDB_SIZE);
	
	HostMatrixFloat3 atom_positions(natom), point_positions(points);
	HostMatrixFloat factor_a(total_funcs, MAX_CONTRACTIONS), factor_c(total_funcs, MAX_CONTRACTIONS), rmm(m, nco), wang(points);
	HostMatrixUInt types(natom), nuc(total_funcs_div), contractions(total_funcs_div);
	
	// REVISAR: nuc: imagen y dominio (especialmente por la parte de * 3 y * 6)

	printf("%i atoms\n", natom);
	for (unsigned int i = 0; i < natom; i++) {
		atom_positions.data[i] = make_float3(r[FORTRAN_MAX_ATOMS * 0 + i], r[i + FORTRAN_MAX_ATOMS * 1], r[i + FORTRAN_MAX_ATOMS * 2]);
		printf("Pos(%i): %f %f %f, Types(%i): %i\n", i, atom_positions.data[i].x, atom_positions.data[i].y, atom_positions.data[i].z, i, Iz[i]);		
		types.data[i] = Iz[i] - 1;
	}
	
	printf("ns: %i, np: %i, nd: %i, Total_Funcs: %i\n", num_funcs.x, num_funcs.y, num_funcs.z, total_funcs);
	{
		uint inc = 1;
		uint i, j;
		for (i = 0, j = 0; i < total_funcs; i += inc, j++) {
			if (i == num_funcs.x) inc = 3;
			else if (i == num_funcs.x + num_funcs.y) inc = 6;

			printf("i: %i, j: %i\n", i, j);
			printf("Nuc(%i) = %i\n", i, Nuc[i] - 1);
			printf("ncont(%i) = %i\n", i, ncont[i]);
			nuc.data[j] = Nuc[i] - 1;
			contractions.data[j] = ncont[i];
			
			for (unsigned int k = 0; k < ncont[i]; k++) {
				factor_a.data[j * MAX_CONTRACTIONS + k] = a[FORTRAN_NG * k + i];
				factor_c.data[j * MAX_CONTRACTIONS + k] = c[FORTRAN_NG * k + i];
				printf("cont: %i, a: %f, c: %f\n", k, factor_a.data[j * MAX_CONTRACTIONS + k], factor_c.data[j * MAX_CONTRACTIONS + k]);
			}			
		}
	}
	
	printf("NCO: %i, M: %i, Iexch: %i\n", nco, total_funcs, Iexch);
	{
		uint k = m18 - 1;
		for (unsigned int i = 0; i < total_funcs; i++) {
			for (unsigned int j = 0; j < nco; j++) {
				rmm.data[i * nco + j] = RMM[k];
				printf("rmm(%i,%i): %.30e (%i)\n", i, j, RMM[k], k);
				k++;	
			}
		}
	}

	printf("Puntos (grilla %i):\n", igrid);
	const double* real_e = (igrid == 1 ? e : e3);
	const double* real_wang = (igrid == 1 ? fort_wang : fort_wang3);
	for (unsigned int i = 0; i < points; i++) {
		wang.data[i] = real_wang[i];
		point_positions.data[i] = make_float3(real_e[0 * points + i], real_e[1 * points + i], real_e[2 * points + i]);
		printf("wang: %f, e: (%f,%f,%f)\n", wang.data[i], point_positions.data[i].x, point_positions.data[i].y, point_positions.data[i].z);
	}
	
	HostMatrixFloat energy(1);
	calc_energy(atom_positions, types, igrid, point_positions, energy, wang,
							Iexch, Ndens, nco, num_funcs_div, nuc, contractions, norm, factor_a, factor_c, rmm);
	Exc = energy.data[0];
	printf("Exc: %f\n", energy.data[0]);
}

/**
 * Host <-> CUDA Communication function
 */

void calc_energy(const HostMatrixFloat3& atom_positions, const HostMatrixUInt& types, uint grid_type,
								 const HostMatrixFloat3& point_positions, HostMatrixFloat& energy, const HostMatrixFloat& wang,
								 uint Iexch, bool Ndens, uint nco, uint3 num_funcs, const HostMatrixUInt& nuc,
								 const HostMatrixUInt& contractions, bool normalize, const HostMatrixFloat& factor_a, const HostMatrixFloat& factor_c,
								 const HostMatrixFloat& rmm)
{	
	const CudaMatrixFloat3 gpu_atom_positions(atom_positions);
	const CudaMatrixUInt gpu_types(types), gpu_nuc(nuc), gpu_contractions(contractions);
	
	dim3 threads(atom_positions.width, MAX_LAYERS, grid_type == 1 ? EXCHNUM_GRIDA_SIZE : EXCHNUM_GRIDB_SIZE);
	dim3 blockSize(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
	dim3 gridSize = divUp(threads, blockSize);
	uint gridSizeZ = gridSize.z;
	gridSize = dim3(gridSize.x, gridSize.y * gridSize.z, 1);

	CudaMatrixFloat gpu_energy(threads.x * threads.y * threads.z), gpu_total_energy, gpu_wang(wang), gpu_factor_a(factor_a), gpu_factor_c(factor_c),
									gpu_rmm(rmm);
	CudaMatrixFloat3 gpu_point_positions(point_positions);

	printf("threads: %i %i %i, blockSize: %i %i %i, gridSize: %i %i %i\n", threads.x, threads.y, threads.z, blockSize.x, blockSize.y, blockSize.z, gridSize.x, gridSize.y / gridSizeZ, gridSizeZ);
		
	if (grid_type == 1) {
		energy_kernel<EXCHNUM_GRIDA_SIZE><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																 gpu_wang.data, gpu_atom_positions.width, Iexch, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																 normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, Ndens);
	}
	else {
		energy_kernel<EXCHNUM_GRIDB_SIZE><<< gridSize, blockSize >>>(gridSizeZ, gpu_atom_positions.data, gpu_types.data, gpu_point_positions.data, gpu_energy.data,
																																 gpu_wang.data, gpu_atom_positions.width, Iexch, nco, num_funcs, gpu_nuc.data, gpu_contractions.data,
																																 normalize, gpu_factor_a.data, gpu_factor_c.data, gpu_rmm.data, Ndens);
	}
	
	/** CPU Accumulation */
	printf("energy %i gpu_energy %i\n", energy.data, gpu_energy.data);
	energy = gpu_energy;
	printf("energy %i gpu_energy %i\n", energy.data, gpu_energy.data);	
	
	double energy_double = 0.0;
	for (unsigned int i = 0; i < threads.x; i++) {
		printf("atomo: %i, capas: %i\n", i, cpu_layers[i]);
		for (unsigned int j = 0; j < cpu_layers[types.data[i]]; j++) {
			for (unsigned int k = 0; k < threads.z; k++) {
				uint idx = index_from3d(threads, dim3(i, j, k));
				double energy_curr = energy.data[idx];
				printf("atomo: %i, capa: %i, punto: %i, valor: %.12e idx: %i\n", i, j, k, energy_curr, idx);
				energy_double += energy_curr;
				energy.data[0] += energy_curr;
			}
		}
	}
	
	//for (unsigned int i = 1; i < energy.elements(); i++) { energy_double += energy.data[i]; energy.data[0] += energy.data[i]; }
	printf("Energy (double): %f\n", energy_double);
	
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
template <unsigned int grid_n>
	__global__ void energy_kernel(uint gridSizeZ, const float3* atom_positions, const uint* types, const float3* point_positions,
																float* energy, const float* wang, const uint atoms_n, uint Iexch, uint nco, uint3 num_funcs,
																const uint* nuc, const uint* contractions, bool normalize, const float* factor_a, const float* factor_c,
																const float* rmm, bool Ndens)
{
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);	
	dim3 pos = index(blockDim, dim3(blockIdx.x, blockIdx.y / gridSizeZ, blockIdx.y % gridSizeZ), threadIdx);
	uint big_index = index_from3d(energySize, pos);

	const uint& atom_i = pos.x;
 	const uint& layer_atom_i = pos.y;
	const uint& point_atom_i = pos.z;
	
	// Hay varios lugares donde se podrian compartir entre threads, solo calculando una vez
	
	// Datos por atomo
	uint atom_i_type = types[atom_i];
	uint atom_i_layers = layers[atom_i_type];
	
	/* load atom and point positions */
	/*__shared__ float3 atom_positions_local[MAX_ATOMS];
	__shared__ float3 point_positions_local[grid_n];
	
	if (layer_atom_i == 0 && point_atom_i == 0) atom_positions_local[atom_i] = atom_positions[atom_i];
	if (atom_i == 0 && layer_atom_i == 0) point_positions_local[point_atom_i] = point_positions[point_atom_i];*/

	__syncthreads();
			
	// atom_positions = atom_positions_local;
	// point_positions = point_positions_local;

	/* skip things that shouldn't be computed */
	// CUIDADO al hacer __syncthreads despues de esto
	if (atom_i >= atoms_n) return;
	if (layer_atom_i >= atom_i_layers) return;	
	if (point_atom_i >= grid_n) return;
		
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
	//printf("atomo: %i, punto: %i, layers: %i, w: %f\n", atom_i, point_atom_i, atom_i_layers, integration_weight);

	float exc_curr = 1.0f;
	float corr_curr = 1.0f;
	float dens = 0.0f;
	
	// float exc_curr, corr_current;
	density_kernel/*<Ndens>*/(dens, num_funcs, nuc, contractions, point_position, atom_positions, normalize, factor_a, factor_c, rmm, nco, big_index);
	pot_kernel(dens, exc_curr, corr_curr, Iexch);
	
	
	/* Numerical Integration */
	float P_total = 0.0f;
	float P_atom_i = 1.0f;
	
	for (uint atomo_j = 0; atomo_j < atoms_n; atomo_j++) {
		float P_curr = 1.0f;

		for (uint atomo_k = 0; atomo_k < atoms_n; atomo_k++) {
			if (atomo_k == atomo_j) continue;
			float rr = distance(atom_positions[atomo_j],atom_positions[atomo_k]);
			float u = distance(point_position,atom_positions[atomo_j]) - distance(point_position, atom_positions[atomo_k]); // revisar que se haga igual que abajo			
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
	
	// output index

	tmp0 = (P_atom_i / P_total) * dens * integration_weight;
		
	energy[index_from3d(energySize, pos)] = tmp0 * exc_curr + tmp0 * corr_curr;
  _EMU(printf("idx: %i dens: %.12e t: %.12e e: %.12e\n", big_index,  dens, tmp0, tmp0 * exc_curr + tmp0 * corr_curr));
}


/* Local Density Functionals */
/*template<bool Ndens>*/ __device__ void density_kernel(float& density, uint3 num_funcs, /*array<double3> distancias_punto_atomos,*/
																												const uint* nuc, const uint* contractions, float3 point_position,
																												const float3* atom_positions, bool normalize, const float* factor_a,
																												const float* factor_c, const float* rmm, uint nco, uint big_index)
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
	const uint& m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6; // sacar luego
	uint func_real = 0;
	
	// DUDA: ng?
	// array<double> W(ng);
	// array<double> F(m);
	
	float W[MAX_NCO];

	for (unsigned int i = 0; i < nco; i++) W[i] = 0.0f;

	// funciones s
	for (uint func = 0; func < funcs_s; func++, func_real++) {
		uint atom_nuc = nuc[func];
		//printf("func: %i, nuc: %i\n", func, atom_nuc);
		float dist = distance2(point_position, atom_positions[atom_nuc]);

		uint func_contractions = contractions[func];
		float func_value = 0.0f;
		
		for (uint contraction = 0; contraction < func_contractions; contraction++) {
		  float rexp = factor_a[func * MAX_CONTRACTIONS + contraction] * dist;
			//printf("func: %i, cont: %i, rexp: %.12e, c: %.12e\n", func, contraction, rexp, factor_c[func * MAX_CONTRACTIONS + contraction]);
			if (rexp > 30.0f) continue;
			func_value += expf(-rexp) * factor_c[func * MAX_CONTRACTIONS + contraction];
		}
		
		//if (big_index == 13624) printf("F(%i) s: %.12e\n", func_real, func_value);
		
		// printf("s: func: %i, value: %.12e\n", func, func_value);
		for (uint i = 0; i < nco; i++) {
			//if (big_index == 13624) printf("W[%i], rmm[%i * %i + %i] = rmm[%i] = %.12e\n", i, i, m, func_real, i * m + func_real, rmm[i * m + func_real]);
			W[i] += rmm[i * m + func_real] * func_value;
		}
	}
	
	// funciones p
	for (uint func = funcs_s; func <  funcs_s + funcs_p; func++, func_real+=3) {
		uint atom_nuc = nuc[func];
		//printf("func: %i, nuc: %i\n", func, atom_nuc);		
		float dist = distance2(point_position, atom_positions[atom_nuc]);
		// float dist = distancias_punto_atomos[atomo_nucleo];
		
		uint func_contractions = contractions[func];
		float func_value[3];
		func_value[0] = func_value[1] = func_value[2] = 0.0f;

		//if (big_index == 3624) printf("cont: %i\n", func_contractions);
		for (uint contraction = 0; contraction < func_contractions; contraction++) {
      //float rexp = factor_a[func][contraction] * dist;
			float rexp = factor_a[func * MAX_CONTRACTIONS + contraction] * dist;
			//printf("func: %i, cont: %i, rexp: %.12e, c: %.12e\n", func, contraction, rexp, factor_c[func * MAX_CONTRACTIONS + contraction]);
      if (rexp > 30.0f) continue;
			
			float t = expf(-rexp) * factor_c[func * MAX_CONTRACTIONS + contraction];
				
			float3 v = (point_position - atom_positions[atom_nuc]) * t;
			for (uint i = 0; i < 3; i++) {
			//if (big_index == 13624) printf("F(%i) p: %.12e\n", func_real + i, elem(v,i));
				func_value[i] += elem(v,i);
			}

      /*for (uint dim = 0; dim < 3; dim++) {
				uint ii = func + dim;
				F[ii] += t * (posicion_punto(dim) - (posiciones_atomos[atomo_nucleo])(dim));
			}*/
			
		}

		//if (big_index == 13624) printf("p: func: %i, value: %.12e %.12e %.12e\n", func, func_value[0], func_value[1], func_value[2]);
		for (uint i = 0; i < nco; i++) {
			for (uint j = 0; j < 3; j++) {
				//if (big_index == 13624) printf("W[%i], rmm[%i * %i + %i + %i] = rmm[%i] = %.12e, f %.12e prod %.12e\n", i, i, m, func_real, j, i * m + func_real + j, rmm[i * m + func_real + j], func_value[j], rmm[i * m + func_real + j] * func_value[j]);
				W[i] += rmm[i * m + func_real + j] * func_value[j];
			}
		}
	}
	
	//float fc = (normalize ? rsqrtf(3.0f) : 1.0f);	// TODO: ver si efectivamente es 1 / sqrtf(3);
	float fc = (normalize ? 1 / sqrtf(3.0f) : 1.0f);
	
	// funciones d
	for (uint func = (funcs_s + funcs_p); func < funcs_total; func++, func_real+=6) {
		uint atom_nuc = nuc[func];
		//printf("func: %i, nuc: %i\n", func, atom_nuc);		
		float dist = distance2(point_position, atom_positions[atom_nuc]);
		//float dist = distancias_punto_atomos(atomo_nucleo);
		
		uint func_contractions = contractions[func];
		float func_value[6];
		func_value[0] = func_value[1] = func_value[2] = func_value[3] = func_value[4] = func_value[5] = 0.0f;
		
		for (uint contraction = 0; contraction < func_contractions; contraction++) {
			float rexp = factor_a[func * MAX_CONTRACTIONS + contraction] * dist;
			//printf("func: %i, cont: %i, rexp: %.12e, c: %.12e\n", func, contraction, rexp, factor_c[func * MAX_CONTRACTIONS + contraction]);
			if (rexp > 30.0f) continue;
			
			float t = expf(-rexp) * factor_c[func * MAX_CONTRACTIONS + contraction];
			
			// TODO: esto quizas se puede reescribir sin tanto ciclo
			float3 v = point_position - atom_positions[atom_nuc];
			
			uint ij = 0;
			for (uint i = 0; i < 3; i++) {
				float t1 = elem(v, i);
				for (uint j = 0; j <= i; j++) {
					float t2 = (i == j ? elem(v, j) * fc : elem(v, j));
					//if (big_index == 13624) printf("F(%i) d: %.12e\n", func_real + i + j, t * t1 * t2);
					func_value[ij] += t * t1 * t2;
					ij++;
					//printf("d(%i)=%.12e\n", i + j, t * t1 * t2);
				}
			}
		}
		
		//if (big_index == 13624) printf("d: func: %i, value: %.12e %.12e %.12e %.12e %.12e %.12e\n", func, func_value[0], func_value[1], func_value[2], func_value[3], func_value[4], func_value[5]);
		for (uint i = 0; i < nco; i++) {
			for (uint j = 0; j < 6; j++) {
				//printf("nco: %i, func: %i, d: rmm(%i)=%.12e\n", i, func_real + j, i * m + func_real + j, rmm[i * m + func_real + j]);
				W[i] += rmm[i * m + func_real + j] * func_value[j];
			}
		}
	}

	for (uint i = 0; i < nco; i++) {
		//if (big_index == 13624) printf("W(%i)=%f\n", i, W[i]);
		density += W[i] * W[i];
	}

	density *= 2.0f;
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
	
/*template<unsigned int Iexch>*/ __device__ void pot_kernel(float dens, float& ex, float& ec, uint Iexch)
{
	// data X alpha

	if (dens == 0) { ex = 0; ec = 0; return; }
	
	float y = powf(dens, 0.333333333333333333);
	float e0 = pot_alpha * y;
	// float v0 = 1.33333333333333 * e0;
	
	switch(Iexch) {
		case 1:
		{
			ex = e0;
			ec = 0;
		}
		break;
		case 2:
		{
			ex = e0;
			float rs = pot_gl / y;
			float x1 = rs / 11.4f;
			
			if (x1 > 1.0f) ec = -0.0333f * (0.5f * x1 - 0.33333333333333f);
			else {
				float t1 = (1.0f + x1 * x1 * x1);
				float t2 = logf(1.0f + 1.0f / x1);
				float t3 = x1 * x1;
        ec = -0.0333f * (t1 * t2 - t3 + 0.5f * x1 - 0.33333333333333f);
			}
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
			
      ec = pot_vosko_a1 * (2 * logf(x1) - t2 + 2 * pot_vosko_b1 / pot_vosko_q * t3 - t4 *
								 (2 * logf(x1 - pot_vosko_x0) - t2 + 2 * (pot_vosko_b1 + 2 * pot_vosko_x0) / pot_vosko_q * t3));
		}
		break;		
	}
}
