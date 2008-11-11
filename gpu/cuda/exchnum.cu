/* -*- mode: c -*- */
#include <cassert>
#include <iostream>
#include <fstream>
#include "common.h"
#include "init.h"
#include "cuda_extra.h"
#include "../matrix.h"
#include "exchnum.h"
#include "gpu_variables.h"
#include "../timer.h"
#include "double.h"
#include "cubes.h"

/** KERNELS **/
#include "functions.h"
#include "energy.h"
#include "rmm.h"
//#include "force.h"

using namespace G2G;
using namespace std;


/**
 * Local Variables
 */

enum ComputationType { COMPUTE_NONE = 0, COMPUTE_ENERGY = 1, UPDATE_RMM = 2, COMPUTE_FORCES = 4 };

/**
 * Methods
 */

void gpu_compute_cube_functions(void)
{
	cout << "<===== computing functions ========>" << endl;
	CudaMatrixFloat3 points_position_gpu;
	CudaMatrixFloat2 factor_ac_gpu;
	CudaMatrixUInt nuc_gpu;
	CudaMatrixUInt contractions_gpu;
	CudaMatrixUInt atom_of_point_gpu;	
	
	Timer t1;
	t1.sync();
	t1.start();
	
	for (list<LittleCube>::iterator it = final_cube.begin(); it != final_cube.end(); ++it) {
		LittleCube& cube = *it;

		/** Load cube points **/
		{
			HostMatrixFloat3 points_position_cpu(cube.number_of_points, 1);
						
			uint i = 0;		
			for (list<Point>::const_iterator p = cube.points.begin(); p != cube.points.end(); ++p, ++i) {
				points_position_cpu.get(i) = make_float3(p->position.x, p->position.y, p->position.z);
			}
			points_position_gpu = points_position_cpu;
		}
		
		/* Load cube functions */
		uint cube_m = cube.s_functions + cube.p_functions * 3 + cube.d_functions * 6;
		uint cube_spd = cube.s_functions + cube.p_functions + cube.d_functions;
		uint4 cube_functions = make_uint4(cube.s_functions, cube.p_functions, cube.d_functions, cube_m);		
		{
			HostMatrixFloat2 factor_ac_cpu(cube_spd, MAX_CONTRACTIONS);
			HostMatrixUInt nuc_cpu(cube_spd, 1), contractions_cpu(cube_spd, 1);
			
			uint i = 0;
			for (set<uint>::const_iterator func = cube.functions.begin(); func != cube.functions.end(); ++func, ++i) {
				nuc_cpu.get(i) = fortran_vars.nucleii.get(*func) - 1;
				contractions_cpu.get(i) = fortran_vars.contractions.get(*func);
				assert(contractions_cpu.get(i) <= MAX_CONTRACTIONS);
				
				for (unsigned int k = 0; k < contractions_cpu.get(i); k++)
					factor_ac_cpu.get(i, k) = make_float2(fortran_vars.a_values.get(*func, k), fortran_vars.c_values.get(*func, k));
			}

			factor_ac_gpu = factor_ac_cpu;
			nuc_gpu = nuc_cpu;
			contractions_gpu = contractions_cpu;
		}
		
		/** Compute Functions **/		
		cube.function_values.resize(cube.number_of_points, cube_functions.w);
		
		dim3 threads(cube.number_of_points);
		dim3 threadBlock(FUNCTIONS_BLOCK_SIZE);
		dim3 threadGrid = divUp(threads, threadBlock);		
		
		gpu_compute_functions<<<threadGrid, threadBlock>>>(points_position_gpu.data, cube.number_of_points,
				contractions_gpu.data, factor_ac_gpu.data, nuc_gpu.data, cube.function_values.data, cube_functions, cube_spd);
		cudaAssertNoError("compute_functions");

		/*HostMatrixFloat functions_cpu(cube.function_values);
		uint i = 0;		
		for (list<Point>::const_iterator p = cube.points.begin(); p != cube.points.end(); ++p, ++i) {
			for (uint j = 0; j < cube_m; j++) {
				cout << "* point (" << p->position.x << "," << p->position.y << "," << p->position.z << ") - F(" << j << ")=" << functions_cpu.get(i, j) << endl;
			}
		}*/
	}	
	
	t1.sync();
	t1.stop();
	cout << "TIMER: funcs: " << t1 << endl;
}

extern "C" void gpu_solve_cubes_(const unsigned int& computation_type_n, double& Exc, double* fort_forces)
{
	/*** determine what to compute ***/
	uint computation_type;
	cout << "<======= GPU Solve Cubes ========>" << endl;
	switch(computation_type_n) {
		case 0: computation_type = COMPUTE_ENERGY; break;
		case 1: computation_type = UPDATE_RMM; break;
		case 2: computation_type = COMPUTE_ENERGY | COMPUTE_FORCES; break;		
	}	

//	#ifdef _DEBUG
//	computation_type |= COMPUTE_ENERGY;
//	#endif

	cout << "calcular energias: " << ((computation_type & COMPUTE_ENERGY) ? "si" : "no") << " ";
	cout << "actualizar rmm: " << ((computation_type & UPDATE_RMM) ? "si" : "no") << " ";
 	cout << "calcular fuerzas: " << ((computation_type & COMPUTE_FORCES) ? "si" : "no") << endl;
		
	/*** Computo sobre cada cubo ****/
	CudaMatrixFloat4 points_position_weight_gpu;
	CudaMatrixFloat rdm_gpu;
	
	double total_energy = 0.0;
	
	/*ofstream points_salida("points", ios_base::app);
	points_salida << "paso" << endl;*/
	Timer t1;
	t1.sync();
	t1.start();								
	
	for (list<LittleCube>::const_iterator it = final_cube.begin(); it != final_cube.end(); ++it) {
		const LittleCube& cube = *it;
				
		/** Load cube points **/
		{
			HostMatrixFloat4 points_position_weight_cpu(cube.number_of_points, 1);
						
			uint i = 0;		
			for (list<Point>::const_iterator p = cube.points.begin(); p != cube.points.end(); ++p, ++i) {
				points_position_weight_cpu.get(i) = make_float4(p->position.x, p->position.y, p->position.z, p->weight);
			}
			points_position_weight_gpu = points_position_weight_cpu;
		}
		
		/** Load cube functions **/
		uint cube_m = cube.s_functions + cube.p_functions * 3 + cube.d_functions * 6;
		
		//cout << "func " << cube_m << " ptos " << cube.number_of_points << endl;
		//points_salida << cube.number_of_points * cube_m * cube_m << endl;
		
		/* load RDM */
		{
			HostMatrixFloat rdm_cpu;
			rdm_cpu.resize(cube_m, fortran_vars.nco);
			for (unsigned int i = 0; i < fortran_vars.nco; i++) {
				uint j = 0;
				for (set<uint>::const_iterator func = cube.functions.begin(); func != cube.functions.end(); ++func) {
				//i * fortran_vars.m + *func + k
					if (*func < fortran_vars.s_funcs) {
						rdm_cpu.get(j, i) = fortran_vars.rmm_input.get(*func, i);
						j++;
					}
					else if (*func < (fortran_vars.s_funcs + fortran_vars.p_funcs * 3)) {
						for (uint k = 0; k < 3; k++, j++) { rdm_cpu.get(j, i) = fortran_vars.rmm_input.get(*func + k, i); }
					}
					else {
						for (uint k = 0; k < 6; k++, j++) { rdm_cpu.get(j, i) = fortran_vars.rmm_input.get(*func + k, i); }
					}
				}
			}
			rdm_gpu = rdm_cpu;
		}
							
		dim3 threads(cube.number_of_points);
		dim3 threadBlock, threadGrid;
		threadBlock = dim3(DENSITY_BLOCK_SIZE);
		threadGrid = divUp(threads, threadBlock);
		
		/*** Compute Energy ***/
		if (computation_type & COMPUTE_ENERGY) {			
			CudaMatrixFloat energy_gpu(cube.number_of_points);
			gpu_compute_density<true><<<threadGrid, threadBlock>>>(energy_gpu.data, points_position_weight_gpu.data, cube.number_of_points, rdm_gpu.data, cube.function_values.data, cube_m);
			cudaAssertNoError("compute_density");
					
					
			HostMatrixFloat energy_cpu(energy_gpu);
			for (uint i = 0; i < cube.number_of_points; i++) total_energy += energy_cpu.get(i);
		}
		else if (computation_type & UPDATE_RMM) {
			/*** Compute Factor necessary for RMM update ***/
			CudaMatrixFloat rmm_factor_gpu(cube.number_of_points);
			gpu_compute_density<false><<<threadGrid, threadBlock>>>(rmm_factor_gpu.data, points_position_weight_gpu.data, cube.number_of_points, rdm_gpu.data, cube.function_values.data, cube_m);
			cudaAssertNoError("compute_density");
			
			HostMatrixFloat rmm_factor_cpu(rmm_factor_gpu);
					
			/*** Compute RMM update ***/
			threads = dim3(cube_m, cube_m);
			threadBlock = dim3(RMM_BLOCK_SIZE_X, RMM_BLOCK_SIZE_Y);
			threadGrid = divUp(threads, threadBlock);
			
			CudaMatrixFloat rmm_output_gpu((cube_m * (cube_m + 1))/2);
			gpu_update_rmm<<<threadGrid, threadBlock>>>(rmm_factor_gpu.data, cube.number_of_points, rmm_output_gpu.data, cube.function_values.data, cube_m);
			cudaAssertNoError("update_rmm");
			
			HostMatrixFloat rmm_output_cpu(rmm_output_gpu);

			/*** Contribute this RMM to the total RMM ***/
			uint small_index = 0;
			for (set<uint>::iterator it_fi = cube.functions.begin(); it_fi != cube.functions.end(); ++it_fi) {
				uint fi_advance;
				if (*it_fi < fortran_vars.s_funcs) fi_advance = 1;
				else if (*it_fi < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) fi_advance = 3;
				else fi_advance = 6;
				
				for (uint i = 0; i < fi_advance; i++) {
					for (set<uint>::iterator it_fj = cube.functions.begin(); it_fj != cube.functions.end(); ++it_fj) {					
						uint fj_advance;
						if (*it_fj < fortran_vars.s_funcs) fj_advance = 1;
						else if (*it_fj < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) fj_advance = 3;
						else fj_advance = 6;
					
						for (uint j = 0; j < fj_advance; j++) {
							uint fi = *it_fi + i; uint fj = *it_fj + j;
							if (fi > fj) continue;
							uint big_index = (fi * fortran_vars.m - (fi * (fi - 1)) / 2) + (fj - fi);
							fortran_vars.rmm_output.get(big_index) += rmm_output_cpu.get(small_index);
							small_index++;
						}					
					}
				}
			}
		}		
	}
	
	t1.sync();
	t1.stop();
	cout << "density/rmm " << t1 << endl;
	
	if (computation_type & COMPUTE_ENERGY) {
		cout << "total energy: " << total_energy << endl;
		Exc = total_energy;
	}		
}
