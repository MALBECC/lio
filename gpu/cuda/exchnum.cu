/* -*- mode: c -*- */
#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "common.h"
#include "init.h"
#include "cuda_extra.h"
#include "../matrix.h"
#include "exchnum.h"
#include "gpu_variables.h"
#include "../timer.h"
#include "double.h"
#include "cubes.h"

#define OLD_DENSITY_KERNEL 1

/** KERNELS **/
#include "functions.h"
#include "energy.h"
#include "rmm.h"
#include "force.h"

using namespace G2G;
using namespace std;


#define COMPUTE_RMM 					0
#define COMPUTE_ENERGY_ONLY		1
#define COMPUTE_ENERGY_FORCE	2
#define COMPUTE_FORCE_ONLY		3


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
		cube.function_values.resize(COALESCED_DIMENSION(cube_functions.w), cube.number_of_points);
		if (fortran_vars.do_forces) cube.gradient_values.resize(cube_functions.w, cube.number_of_points);
		
		dim3 threads(cube.number_of_points);
		dim3 threadBlock(FUNCTIONS_BLOCK_SIZE);
		dim3 threadGrid = divUp(threads, threadBlock);		

		//cout << "points: " << threads.x << " " << threadGrid.x << " " << threadBlock.x << endl;
		
		if (fortran_vars.do_forces)
			gpu_compute_functions<true><<<threadGrid, threadBlock>>>(points_position_gpu.data, cube.number_of_points, contractions_gpu.data, factor_ac_gpu.data, nuc_gpu.data, cube.function_values.data, cube.gradient_values.data, cube_functions, cube_spd);
		else
			gpu_compute_functions<false><<<threadGrid, threadBlock>>>(points_position_gpu.data, cube.number_of_points, contractions_gpu.data, factor_ac_gpu.data, nuc_gpu.data, cube.function_values.data, cube.gradient_values.data, cube_functions, cube_spd);

		cudaAssertNoError("compute_functions");

#if 0
		if (fortran_vars.grid_type == BIG_GRID) {
			cout << "s_funcs: " << cube.s_functions << " p_funcs " << cube.p_functions << " d_funcs " << cube.d_functions << endl;
			HostMatrixFloat functions_cpu(cube.function_values);
			HostMatrixFloat3 gradients_cpu(cube.gradient_values);
			uint i = 0;		
			for (list<Point>::const_iterator p = cube.points.begin(); p != cube.points.end(); ++p, ++i) {
				uint func_idx = 0;
				for (set<uint>::const_iterator func = cube.functions.begin(); func != cube.functions.end(); ++func, ++func_idx) {
					if (fortran_vars.nucleii.get(*func) - 1 != 0) continue;
					if (func_idx < cube.s_functions)
						cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(func_idx, i).x << " "  << gradients_cpu.get(func_idx, i).y << " " << gradients_cpu.get(func_idx, i).z << " F " << functions_cpu.get(func_idx, i) << " " << func_idx << endl;
					else if (func_idx < cube.p_functions + cube.s_functions) {
						uint p_idx = 3 * (func_idx - cube.s_functions) + cube.s_functions;
						for (uint j = 0; j < 3; j++)
							cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(p_idx + j, i).x << " "  << gradients_cpu.get(p_idx + j, i).y << " " << gradients_cpu.get(p_idx + j, i).z << " F " << functions_cpu.get(p_idx + j, i) << " " << p_idx + j << endl;
					}
					else {
						uint s_idx = cube.s_functions + cube.p_functions * 3 + 6 * (func_idx - cube.s_functions - cube.p_functions);
						for (uint j = 0; j < 6; j++)
							cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(s_idx + j, i).x << " "  << gradients_cpu.get(s_idx + j, i).y << " " << gradients_cpu.get(s_idx + j, i).z << " F " << functions_cpu.get(s_idx + j, i) << " " << s_idx + j << endl;

					}
//				cout << "* point " << p->position.x << " " << p->position.y << " " << p->position.z << " " << functions_cpu.get(p_idx, i) << endl;
				}
			}
		}
#endif
	}	
	
	t1.sync();
	t1.stop();
	cout << "TIMER: funcs: " << t1 << endl;
}

extern "C" void gpu_solve_cubes_(uint& computation_type, double* fort_energy_ptr, double* fort_forces_ptr)
{
	cout << "<================ calculo de: [";
	switch(computation_type) {
		case COMPUTE_ENERGY_ONLY: cout << "energia"; break;
		case COMPUTE_RMM: cout << "rmm"; break;
		case COMPUTE_FORCE_ONLY: cout << "fuerzas"; break;
		case COMPUTE_ENERGY_FORCE: cout << "energia+fuerzas"; break;
	}
	cout << "] ==========>" << endl;
	
	Timer t1;
	t1.sync();
	t1.start();									
		
	/*** Computo sobre cada cubo ****/
	CudaMatrixFloat4 points_position_weight_gpu;
	CudaMatrixFloat rdm_gpu;
	CudaMatrixUInt nuc_gpu;
	CudaMatrixFloat2 factor_ac_gpu;
	CudaMatrixUInt contractions_gpu;
	

	FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, FORTRAN_MAX_ATOMS);
	
	double total_energy = 0.0;
		
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
		uint cube_spd = cube.s_functions + cube.p_functions + cube.d_functions;
		uint4 cube_functions = make_uint4(cube.s_functions, cube.p_functions, cube.d_functions, cube_m);
		//cout << "s: " << cube_functions.x << " p: " << cube_functions.y << " d: " << cube_functions.z << endl;
		
		/* load RDM */
		{
			HostMatrixFloat rdm_cpu;
			//cout << "cube_m " << (cube_m + (16 - cube_m % 16)) << endl;
			rdm_cpu.resize(COALESCED_DIMENSION(cube_m), fortran_vars.nco);
			for (unsigned int i = 0; i < fortran_vars.nco; i++) {
				uint j = 0;
				for (set<uint>::const_iterator func = cube.functions.begin(); func != cube.functions.end(); ++func) {
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

#if !OLD_DENSITY_KERNEL
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
#endif
							
		dim3 threads(cube.number_of_points);
		dim3 threadBlock, threadGrid;
#if OLD_DENSITY_KERNEL
		threadBlock = dim3(DENSITY_BLOCK_SIZE);
#else
		threadBlock = dim3(DENSITY_BLOCK_SIZE_X);
#endif
		threadGrid = divUp(threads, threadBlock);

		/* compute energy */
		if (computation_type == COMPUTE_ENERGY_ONLY) {
			CudaMatrixFloat energy_gpu(cube.number_of_points);
#if OLD_DENSITY_KERNEL
			gpu_compute_density<true, false><<<threadGrid, threadBlock>>>(energy_gpu.data, NULL, points_position_weight_gpu.data, cube.number_of_points, rdm_gpu.data, cube.function_values.data, NULL, NULL, NULL, 0, cube_functions);
#else
			gpu_compute_density2<true, false><<<threadGrid, threadBlock>>>(energy_gpu.data, NULL, points_position_weight_gpu.data, cube.number_of_points, rdm_gpu.data, nuc_gpu.data, cube_functions, cube_spd, contractions_gpu.data, factor_ac_gpu.data);
#endif
			cudaAssertNoError("compute_density");

			HostMatrixFloat energy_cpu(energy_gpu);
			for (uint i = 0; i < cube.number_of_points; i++) { total_energy += energy_cpu.get(i); }
		}
		/* compute necessary factor **/
		else if (computation_type == COMPUTE_RMM) {
			CudaMatrixFloat rmm_factor_gpu(cube.number_of_points);
			gpu_compute_density<false, false><<<threadGrid, threadBlock>>>(NULL, rmm_factor_gpu.data, points_position_weight_gpu.data, cube.number_of_points, rdm_gpu.data, cube.function_values.data, NULL, NULL, NULL, 0, cube_functions);
			cudaAssertNoError("compute_density");

			/*** Compute RMM update ***/
			threads = dim3(cube_m, cube_m);
			threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
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
		/* compute forces */
		else {
			map<uint, uint> nuc_map;
			{
				HostMatrixUInt nuc_cpu(cube_spd, 1);			
				uint i = 0;
				uint curr_idx = 0;

				for (set<uint>::const_iterator func = cube.functions.begin(); func != cube.functions.end(); ++func, ++i) {
					uint func_nuc = fortran_vars.nucleii.get(*func) - 1;
					map<uint, uint>::iterator nuc_idx = nuc_map.find(func_nuc);
					if (nuc_idx == nuc_map.end()) {
						nuc_map[func_nuc] = curr_idx;
						nuc_cpu.get(i) = curr_idx;
						curr_idx++;
					}
					else nuc_cpu.get(i) = nuc_idx->second;
				}
				nuc_gpu = nuc_cpu;
			}

			CudaMatrixFloat force_factor_gpu(cube.number_of_points);
			CudaMatrixFloat energy_gpu;

			/* energy may be needed at this step */
			CudaMatrixFloat3 density_deriv(COALESCED_DIMENSION(cube.nucleii.size()), cube.number_of_points);
			if (computation_type == COMPUTE_ENERGY_FORCE) {
				energy_gpu.resize(cube.number_of_points);
				gpu_compute_density<true, true><<<threadGrid, threadBlock>>>(energy_gpu.data, force_factor_gpu.data, points_position_weight_gpu.data, cube.number_of_points, rdm_gpu.data, cube.function_values.data, cube.gradient_values.data, density_deriv.data, nuc_gpu.data, cube.nucleii.size(), cube_functions);

				HostMatrixFloat energy_cpu(energy_gpu);
				for (uint i = 0; i < cube.number_of_points; i++) { total_energy += energy_cpu.get(i); }
			}
			else
				gpu_compute_density<false, true><<<threadGrid, threadBlock>>>(energy_gpu.data, force_factor_gpu.data, points_position_weight_gpu.data, cube.number_of_points, rdm_gpu.data, cube.function_values.data, cube.gradient_values.data, density_deriv.data, nuc_gpu.data, cube.nucleii.size(), cube_functions);

			cudaAssertNoError("compute_density");

			threads = dim3(cube.nucleii.size());
			threadBlock = dim3(FORCE_BLOCK_SIZE);
			threadGrid = divUp(threads, threadBlock);

#if 0
			/* TODO: remover DEBUG */
			HostMatrixFloat3 density_deriv_cpu(density_deriv);
			uint point_idx = 0;
			for (list<Point>::const_iterator it = cube.points.begin(); it != cube.points.end(); ++it, point_idx++) {
				for (map<uint, uint>::iterator nuc_it = nuc_map.begin(); nuc_it != nuc_map.end(); ++nuc_it) {
					cout << "factor: " << point_idx << " " << nuc_it->second << " " <<
            density_deriv_cpu.get(nuc_it->second, point_idx).x << " " <<
            density_deriv_cpu.get(nuc_it->second, point_idx).y << " " <<
            density_deriv_cpu.get(nuc_it->second, point_idx).z << " " <<
          endl;
				}
			}
			/* DEBUG */
#endif
			
			CudaMatrixFloat3 gpu_forces(cube.nucleii.size());
			gpu_compute_forces<<<threadGrid, threadBlock>>>(cube.number_of_points, force_factor_gpu.data, density_deriv.data, gpu_forces.data, cube.nucleii.size());
			cudaAssertNoError("gpu_compute_forces");

			HostMatrixFloat3 cpu_forces(gpu_forces);
			for (map<uint, uint>::iterator nuc_it = nuc_map.begin(); nuc_it != nuc_map.end(); ++nuc_it) {
				float3 atom_force = cpu_forces.get(nuc_it->second);
        //cout << "atom force: " << atom_force.x << " " << atom_force.y << " " << atom_force.z << endl;
				fort_forces.get(nuc_it->first, 0) += atom_force.x;
				fort_forces.get(nuc_it->first, 1) += atom_force.y;
				fort_forces.get(nuc_it->first, 2) += atom_force.z;
			}
		}
	}
		
	/** pass results to fortran */
	if (computation_type == COMPUTE_ENERGY_ONLY || computation_type == COMPUTE_ENERGY_FORCE) {
		cout << "total energy: " << total_energy << endl;
		*fort_energy_ptr = total_energy;
	}
	
	t1.sync();
	t1.stop();
	cout << "TIMER: gpu_solve_cubes " << t1 << endl;
}
