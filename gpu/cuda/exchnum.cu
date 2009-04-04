/* -*- mode: c -*- */
#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "../common.h"
#include "../init.h"
#include "cuda_extra.h"
#include "../matrix.h"
#include "gpu_variables.h"
#include "../timer.h"
#include "double.h"
#include "../partition.h"

#define OLD_DENSITY_KERNEL 1

/** KERNELS **/
#include "functions.h"
#include "energy.h"
#include "rmm.h"
#include "force.h"
#include "weight.h"

using namespace G2G;
using namespace std;


#define COMPUTE_RMM 					0
#define COMPUTE_ENERGY_ONLY		1
#define COMPUTE_ENERGY_FORCE	2
#define COMPUTE_FORCE_ONLY		3


/*******************************
 * Cube Functions
 *******************************/
void gpu_compute_group_functions(void)
{
	cout << "<===== computing functions ========>" << endl;
	CudaMatrixFloat3 points_position_gpu;
	CudaMatrixFloat2 factor_ac_gpu;
	CudaMatrixUInt nuc_gpu;
	CudaMatrixUInt contractions_gpu;
	
	Timer t1;
	t1.sync();
	t1.start();
	
	for (list<PointGroup>::iterator it = final_partition.begin(); it != final_partition.end(); ++it) {
		PointGroup& group = *it;
		/** Load points from group **/
		{
			HostMatrixFloat3 points_position_cpu(group.number_of_points, 1);
						
			uint i = 0;		
			for (list<Point>::const_iterator p = group.points.begin(); p != group.points.end(); ++p, ++i) {
				points_position_cpu.get(i) = make_float3(p->position.x, p->position.y, p->position.z);
			}
			points_position_gpu = points_position_cpu;
		}
		
		/* Load group functions */
		uint group_m = group.s_functions + group.p_functions * 3 + group.d_functions * 6;
		uint group_spd = group.s_functions + group.p_functions + group.d_functions;
		uint4 group_functions = make_uint4(group.s_functions, group.p_functions, group.d_functions, group_m);
    //cout << "m: " << group_m << " spd " << group_spd << endl;
		{
			HostMatrixFloat2 factor_ac_cpu(group_spd, MAX_CONTRACTIONS);
			HostMatrixUInt nuc_cpu(group_spd, 1), contractions_cpu(group_spd, 1);
			
			uint i = 0;
			for (set<uint>::const_iterator func = group.functions.begin(); func != group.functions.end(); ++func, ++i) {
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
    group.function_values.resize(COALESCED_DIMENSION(group.number_of_points), group_functions.w);
		if (fortran_vars.do_forces) group.gradient_values.resize(group_functions.w, group.number_of_points);
		
		dim3 threads(group.number_of_points);
		dim3 threadBlock(FUNCTIONS_BLOCK_SIZE);
		dim3 threadGrid = divUp(threads, threadBlock);		

		//cout << "points: " << threads.x << " " << threadGrid.x << " " << threadBlock.x << endl;
		
		if (fortran_vars.do_forces)
			gpu_compute_functions<true><<<threadGrid, threadBlock>>>(points_position_gpu.data, group.number_of_points, contractions_gpu.data, factor_ac_gpu.data, nuc_gpu.data, group.function_values.data, group.gradient_values.data, group_functions, group_spd);
		else
			gpu_compute_functions<false><<<threadGrid, threadBlock>>>(points_position_gpu.data, group.number_of_points, contractions_gpu.data, factor_ac_gpu.data, nuc_gpu.data, group.function_values.data, group.gradient_values.data, group_functions, group_spd);

		cudaAssertNoError("compute_functions");

#if 0
		if (fortran_vars.grid_type == BIG_GRID) {
			cout << "s_funcs: " << group.s_functions << " p_funcs " << group.p_functions << " d_funcs " << group.d_functions << endl;
			HostMatrixFloat functions_cpu(group.function_values);
			HostMatrixFloat3 gradients_cpu(group.gradient_values);
			uint i = 0;		
			for (list<Point>::const_iterator p = group.points.begin(); p != group.points.end(); ++p, ++i) {
				uint func_idx = 0;
				for (set<uint>::const_iterator func = group.functions.begin(); func != group.functions.end(); ++func, ++func_idx) {
					if (fortran_vars.nucleii.get(*func) - 1 != 0) continue;
					if (func_idx < group.s_functions)
						cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(func_idx, i).x << " "  << gradients_cpu.get(func_idx, i).y << " " << gradients_cpu.get(func_idx, i).z << " F " << functions_cpu.get(func_idx, i) << " " << func_idx << endl;
					else if (func_idx < group.p_functions + group.s_functions) {
						uint p_idx = 3 * (func_idx - group.s_functions) + group.s_functions;
						for (uint j = 0; j < 3; j++)
							cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(p_idx + j, i).x << " "  << gradients_cpu.get(p_idx + j, i).y << " " << gradients_cpu.get(p_idx + j, i).z << " F " << functions_cpu.get(p_idx + j, i) << " " << p_idx + j << endl;
					}
					else {
						uint s_idx = group.s_functions + group.p_functions * 3 + 6 * (func_idx - group.s_functions - group.p_functions);
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

/*******************************
 * Cube Weights
 *******************************/

void gpu_compute_group_weights(PointGroup& group)
{
  CudaMatrixFloat4 point_positions_gpu;
  {
    HostMatrixFloat4 points_positions_cpu(group.number_of_points, 1);

		uint i = 0;
		for (list<Point>::const_iterator p = group.points.begin(); p != group.points.end(); ++p, ++i) {
			points_positions_cpu.get(i) = make_float4(p->position.x, p->position.y, p->position.z, p->atom);
		}
		point_positions_gpu = points_positions_cpu;
	}

  CudaMatrixFloat weights_gpu(group.number_of_points);
  dim3 threads(group.number_of_points);
  dim3 blockSize(WEIGHT_BLOCK_SIZE);
  dim3 gridSize = divUp(threads, blockSize);
  gpu_compute_weights<<<gridSize,blockSize>>>(group.number_of_points, point_positions_gpu.data, weights_gpu.data);
  cudaAssertNoError("compute_weights");

  HostMatrixFloat weights_cpu(weights_gpu);
  uint i = 0;
  for (list<Point>::iterator p = group.points.begin(); p != group.points.end(); ++p, ++i) {
    p->weight *= weights_cpu.get(i);
  }
}

/********************************
 * Solve Cubes
 ********************************/
extern "C" void gpu_solve_groups_(uint& computation_type, double* fort_energy_ptr, double* fort_forces_ptr)
{
	cout << "<================ calculo de: [";
	switch(computation_type) {
		case COMPUTE_ENERGY_ONLY: cout << "energia"; break;
		case COMPUTE_RMM: cout << "rmm"; break;
		case COMPUTE_FORCE_ONLY: cout << "fuerzas"; break;
		case COMPUTE_ENERGY_FORCE: cout << "energia+fuerzas"; break;
	}
	cout << "] ==========>" << endl;
	
	Timer t_total;
	t_total.sync();
	t_total.start();
		
	/*** Computo sobre cada cubo ****/
	CudaMatrixFloat point_weights_gpu;
	CudaMatrixFloat rdm_gpu;
	CudaMatrixUInt nuc_gpu;

  Timer t_density, t_rmm, t_energy;

	FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, FORTRAN_MAX_ATOMS);
	
	double total_energy = 0.0;
		
	for (list<PointGroup>::const_iterator it = final_partition.begin(); it != final_partition.end(); ++it) {
		const PointGroup& group = *it;
				
		/** Load points from group **/
		{
			HostMatrixFloat point_weights_cpu(group.number_of_points, 1);

			uint i = 0;		
			for (list<Point>::const_iterator p = group.points.begin(); p != group.points.end(); ++p, ++i) {
				point_weights_cpu.get(i) = p->weight;
			}
			point_weights_gpu = point_weights_cpu;
		}
		
		/** Load functions from group **/
		uint group_m = group.s_functions + group.p_functions * 3 + group.d_functions * 6;
		uint group_spd = group.s_functions + group.p_functions + group.d_functions;
		uint4 group_functions = make_uint4(group.s_functions, group.p_functions, group.d_functions, group_m);
		
		/* load RDM */
		{
			HostMatrixFloat rdm_cpu(COALESCED_DIMENSION(group_m), fortran_vars.nco);
			for (unsigned int i = 0; i < fortran_vars.nco; i++) {
				uint j = 0;
				for (set<uint>::const_iterator func = group.functions.begin(); func != group.functions.end(); ++func) {
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
							
		dim3 threads(group.number_of_points);
		dim3 threadBlock, threadGrid;
		threadBlock = dim3(DENSITY_BLOCK_SIZE);
		threadGrid = divUp(threads, threadBlock);

		/* compute energy */
		if (computation_type == COMPUTE_ENERGY_ONLY) {
      t_energy.start_and_sync();
			CudaMatrixFloat energy_gpu(group.number_of_points);
			gpu_compute_density<true, false><<<threadGrid, threadBlock>>>(energy_gpu.data, NULL, point_weights_gpu.data, group.number_of_points, rdm_gpu.data, group.function_values.data, NULL, NULL, NULL, 0, group_functions);
			cudaAssertNoError("compute_density");
      t_energy.pause_and_sync();

			HostMatrixFloat energy_cpu(energy_gpu);
			for (uint i = 0; i < group.number_of_points; i++) { total_energy += energy_cpu.get(i); }
		}
		/* compute necessary factor **/
		else if (computation_type == COMPUTE_RMM) {
			CudaMatrixFloat rmm_factor_gpu(group.number_of_points);
      t_density.start_and_sync();
			gpu_compute_density<false, false><<<threadGrid, threadBlock>>>(NULL, rmm_factor_gpu.data, point_weights_gpu.data, group.number_of_points, rdm_gpu.data, group.function_values.data, NULL, NULL, NULL, 0, group_functions);
			cudaAssertNoError("compute_density");
      t_density.pause_and_sync();

			/*** Compute RMM update ***/
			threads = dim3(group_m, group_m);
			threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
			threadGrid = divUp(threads, threadBlock);

      CudaMatrixFloat rmm_output_gpu(COALESCED_DIMENSION(group_m), group_m);

      t_rmm.start_and_sync();
			gpu_update_rmm<<<threadGrid, threadBlock>>>(rmm_factor_gpu.data, group.number_of_points, rmm_output_gpu.data, group.function_values.data, group_m);
			cudaAssertNoError("update_rmm");
      t_rmm.pause_and_sync();
			HostMatrixFloat rmm_output_cpu(rmm_output_gpu);

      /*** Contribute this RMM to the total RMM ***/
      uint small_fi = 0;

			for (set<uint>::iterator it_fi = group.functions.begin(); it_fi != group.functions.end(); ++it_fi) {
				uint fi_advance;
				if (*it_fi < fortran_vars.s_funcs) fi_advance = 1;
				else if (*it_fi < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) fi_advance = 3;
				else fi_advance = 6;
				
				for (uint i = 0; i < fi_advance; i++) {

          uint small_fj = 0;
					for (set<uint>::iterator it_fj = group.functions.begin(); it_fj != group.functions.end(); ++it_fj) {
						uint fj_advance;
						if (*it_fj < fortran_vars.s_funcs) fj_advance = 1;
						else if (*it_fj < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) fj_advance = 3;
						else fj_advance = 6;
					
						for (uint j = 0; j < fj_advance; j++) {
							uint fi = *it_fi + i; uint fj = *it_fj + j;
							if (fi > fj) continue;
							uint big_index = (fi * fortran_vars.m - (fi * (fi - 1)) / 2) + (fj - fi);
              //fortran_vars.rmm_output.get(big_index) += rmm_output_cpu.get(small_index);
              //cout << "small (" << cube_m << "): " << small_fj + small_fi << " " << small_fi<< ": " << rmm_output_cpu.get(small_fi, small_fj + small_fi) << endl;
              fortran_vars.rmm_output.get(big_index) += rmm_output_cpu.get(small_fi, small_fj + small_fi);
              small_fj++;
						}					
					}
          small_fi++;
				}
			}
		}
		/* compute forces */
		else {
			map<uint, uint> nuc_map;
			{
				HostMatrixUInt nuc_cpu(group_spd, 1);
				uint i = 0;
				uint curr_idx = 0;

				for (set<uint>::const_iterator func = group.functions.begin(); func != group.functions.end(); ++func, ++i) {
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

			CudaMatrixFloat force_factor_gpu(group.number_of_points);
			CudaMatrixFloat energy_gpu;

			/* energy may be needed at this step */
			CudaMatrixFloat3 density_deriv(COALESCED_DIMENSION(group.nucleii.size()), group.number_of_points);
			if (computation_type == COMPUTE_ENERGY_FORCE) {
				energy_gpu.resize(group.number_of_points);
				gpu_compute_density<true, true><<<threadGrid, threadBlock>>>(energy_gpu.data, force_factor_gpu.data, point_weights_gpu.data, group.number_of_points, rdm_gpu.data, group.function_values.data, group.gradient_values.data, density_deriv.data, nuc_gpu.data, group.nucleii.size(), group_functions);

				HostMatrixFloat energy_cpu(energy_gpu);
				for (uint i = 0; i < group.number_of_points; i++) { total_energy += energy_cpu.get(i); }
			}
			else
				gpu_compute_density<false, true><<<threadGrid, threadBlock>>>(energy_gpu.data, force_factor_gpu.data, point_weights_gpu.data, group.number_of_points, rdm_gpu.data, group.function_values.data, group.gradient_values.data, density_deriv.data, nuc_gpu.data, group.nucleii.size(), group_functions);

			cudaAssertNoError("compute_density");

			threads = dim3(group.nucleii.size());
			threadBlock = dim3(FORCE_BLOCK_SIZE);
			threadGrid = divUp(threads, threadBlock);

#if 0
			/* TODO: remover DEBUG */
			HostMatrixFloat3 density_deriv_cpu(density_deriv);
			uint point_idx = 0;
			for (list<Point>::const_iterator it = group.points.begin(); it != group.points.end(); ++it, point_idx++) {
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
			
			CudaMatrixFloat3 gpu_forces(group.nucleii.size());
			gpu_compute_forces<<<threadGrid, threadBlock>>>(group.number_of_points, force_factor_gpu.data, density_deriv.data, gpu_forces.data, group.nucleii.size());
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
	t_total.stop_and_sync();

	cout << "TIMER: gpu_solve_cubes " << t_total << endl;
  cout << "TIMER: density " << t_density << endl;
  cout << "TIMER: energy " << t_energy << endl;
  cout << "TIMER: rmm: " << t_rmm << endl;
}
