/* -*- mode: c -*- */
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <math_constants.h>
#include "../common.h"
#include "../init.h"
#include "cuda_extra.h"
#include "../matrix.h"
#include "../timer.h"
#include "../partition.h"

/** KERNELS **/
#include "gpu_variables.h"
#include "kernels/pot.h"
#include "kernels/energy.h"
#include "kernels/energy_derivs.h"
#include "kernels/rmm.h"
#include "kernels/weight.h"
#include "kernels/functions.h"
#include "kernels/force.h"

using namespace G2G;
using namespace std;

template<bool lda, bool compute_forces>
void g2g_iteration(bool compute_energy, bool compute_rmm, double* fort_energy_ptr, double* fort_forces_ptr)
{
  double total_energy = 0.0;

  Timer t_total, t_rmm, t_density, t_density_derivs, t_forces, t_resto, t_functions;
	t_total.sync();
	t_total.start();

  //uint max_used_memory = 0;

	/*** Computo sobre cada cubo ****/
	CudaMatrixFloat point_weights_gpu;
	CudaMatrixFloat rmm_input_gpu;
	FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, FORTRAN_MAX_ATOMS);

	for (list<PointGroup>::iterator it = final_partition.begin(); it != final_partition.end(); ++it) {
		PointGroup& group = *it;

    t_functions.start_and_sync();
    /** Compute this group's functions **/
    group.compute_functions<compute_forces, !lda>();
    t_functions.pause_and_sync();

    uint group_m = group.total_functions();
    //cout << "group: " << group.number_of_points << endl;

		/** Load points from group **/
    HostMatrixFloat point_weights_cpu(group.number_of_points, 1);

		uint i = 0;
		for (list<Point>::const_iterator p = group.points.begin(); p != group.points.end(); ++p, ++i) {
			point_weights_cpu(i) = p->weight;
		}
		point_weights_gpu = point_weights_cpu;

    /* compute density/factors */
		dim3 threads(group.number_of_points);
		dim3 threadBlock, threadGrid;
		threadBlock = dim3(DENSITY_BLOCK_SIZE);
		threadGrid = divUp(threads, threadBlock);

    CudaMatrixFloat factors_gpu;
    if (compute_rmm || compute_forces) factors_gpu.resize(group.number_of_points);

    t_density.start_and_sync();
    HostMatrixFloat rmm_input_cpu(COALESCED_DIMENSION(group_m), group_m);
    group.get_rmm_input(rmm_input_cpu);
    rmm_input_gpu = rmm_input_cpu;

		if (compute_energy) {
      CudaMatrixFloat energy_gpu(group.number_of_points);
      gpu_compute_density<true, compute_forces, lda><<<threadGrid, threadBlock>>>(energy_gpu.data, factors_gpu.data, point_weights_gpu.data, group.number_of_points,
        rmm_input_gpu.data, group.function_values.data, group.gradient_values.data, group.hessian_values.data, group_m);
      cudaAssertNoError("compute_density");

      HostMatrixFloat energy_cpu(energy_gpu);
			for (uint i = 0; i < group.number_of_points; i++) { total_energy += energy_cpu(i); } // TODO: hacer con un kernel?
    }
    else {
      gpu_compute_density<false, compute_forces, lda><<<threadGrid, threadBlock>>>(NULL, factors_gpu.data, point_weights_gpu.data, group.number_of_points,
        rmm_input_gpu.data, group.function_values.data, group.gradient_values.data, group.hessian_values.data, group_m);
      cudaAssertNoError("compute_density");
      {
        //HostMatrixFloat factors_cpu(factors_gpu);
        //for (uint i = 0; i < group.number_of_points; i++) { cout << "factor " << factors_cpu(i) << endl;; }
        //factors_cpu.check_values();
      }
    }

    t_density.pause_and_sync();

    /* compute forces */
    if (compute_forces) {
      t_density_derivs.start_and_sync();
      threads = dim3(group.number_of_points);
			threadBlock = dim3(DENSITY_DERIV_BLOCK_SIZE);
			threadGrid = divUp(threads, threadBlock);

      CudaMatrixFloat4 dd_gpu(COALESCED_DIMENSION(group.number_of_points), group.total_nucleii()); dd_gpu.zero();
      CudaMatrixUInt nuc_gpu(group.func2local_nuc);  // TODO: esto en realidad se podria guardar una sola vez durante su construccion

      gpu_compute_density_derivs<<<threadGrid, threadBlock>>>(group.function_values.data, group.gradient_values.data, rmm_input_gpu.data, nuc_gpu.data, dd_gpu.data, group.number_of_points, group_m, group.total_nucleii());
      cudaAssertNoError("density_derivs");
      t_density_derivs.pause_and_sync();

      t_forces.start_and_sync();
      CudaMatrixFloat4 forces_gpu(group.total_nucleii());

      threads = dim3(group.total_nucleii());
			threadBlock = dim3(FORCE_BLOCK_SIZE);
			threadGrid = divUp(threads, threadBlock);
      gpu_compute_forces<<<threadGrid, threadBlock>>>(group.number_of_points, factors_gpu.data, dd_gpu.data, forces_gpu.data, group.total_nucleii());
      cudaAssertNoError("forces");

      HostMatrixFloat4 forces_cpu(forces_gpu);

      for (uint i = 0; i < group.total_nucleii(); ++i) {
				float4 atom_force = forces_cpu(i);
        uint global_nuc = group.local2global_nuc[i];
				fort_forces(global_nuc, 0) += atom_force.x;
				fort_forces(global_nuc, 1) += atom_force.y;
				fort_forces(global_nuc, 2) += atom_force.z;
      }
      t_forces.pause_and_sync();
    }

    /* compute RMM */
    t_rmm.start_and_sync();
    if (compute_rmm) {
			threads = dim3(group_m, group_m);
			threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
			threadGrid = divUp(threads, threadBlock);

      CudaMatrixFloat rmm_output_gpu(COALESCED_DIMENSION(group_m), group_m);
			gpu_update_rmm<<<threadGrid, threadBlock>>>(factors_gpu.data, group.number_of_points, rmm_output_gpu.data, group.function_values.data, group_m);
			cudaAssertNoError("update_rmm");

      /*** Contribute this RMM to the total RMM ***/
      HostMatrixFloat rmm_output_cpu(rmm_output_gpu);
      group.add_rmm_output(rmm_output_cpu);
		}
    t_rmm.pause_and_sync();

    /* clear functions */
    group.function_values.deallocate();
    group.gradient_values.deallocate();
    group.hessian_values.deallocate();
	}

	/** pass results to fortran */
	if (compute_energy) {
		cout << "total energy: " << total_energy << endl;
		*fort_energy_ptr = total_energy;
	}
	t_total.stop_and_sync();
  cout << "iteration: " << t_total << endl;
  cout << "rmm: " << t_rmm << " density: " << t_density << " density_derivs: " << t_density_derivs << " force: " << t_forces << " resto: " << t_resto;
  cout << " funcs: " << t_functions << endl;

  //uint free_memory, total_memory;
  //cudaGetMemoryInfo(free_memory, total_memory);
  //cout << "Maximum used memory: " << (double)max_used_memory / (1024 * 1024) << "MB (" << ((double)max_used_memory / total_memory) * 100.0 << "%)" << endl;
  cudaPrintMemoryInfo();
}

template void g2g_iteration<true, true>(bool compute_energy, bool compute_rmm, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<true, false>(bool compute_energy, bool compute_rmm, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<false, true>(bool compute_energy, bool compute_rmm, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<false, false>(bool compute_energy, bool compute_rmm, double* fort_energy_ptr, double* fort_forces_ptr);

/*******************************
 * Cube Functions
 *******************************/
template <bool forces, bool gga>
void PointGroup::compute_functions(void)
{
  CudaMatrixFloat4 points_position_gpu;
  CudaMatrixFloat2 factor_ac_gpu;
  CudaMatrixUInt nuc_gpu;
  CudaMatrixUInt contractions_gpu;

  /** Load points from group **/
  HostMatrixFloat4 points_position_cpu(number_of_points, 1);
  uint i = 0;
  for (list<Point>::const_iterator p = points.begin(); p != points.end(); ++p, ++i) {
    points_position_cpu(i) = make_float4(p->position.x, p->position.y, p->position.z, 0);
  }
  points_position_gpu = points_position_cpu;

  /* Load group functions */
  uint group_m = s_functions + p_functions * 3 + d_functions * 6;
  uint4 group_functions = make_uint4(s_functions, p_functions, d_functions, group_m);
  HostMatrixFloat2 factor_ac_cpu(COALESCED_DIMENSION(group_m), MAX_CONTRACTIONS);
  HostMatrixUInt nuc_cpu(group_m, 1), contractions_cpu(group_m, 1);

  // TODO: hacer que functions.h itere por total_small_functions()... asi puedo hacer que func2global_nuc sea de tamaño total_functions() y directamente copio esa matriz aca y en otros lados

  uint ii = 0;
  for (i = 0; i < total_functions_simple(); ++i) {
    uint inc = small_function_type(i);

    uint func = local2global_func[i];
    uint this_nuc = func2global_nuc(i);
    uint this_cont = fortran_vars.contractions(func);

    for (uint j = 0; j < inc; j++) {
      nuc_cpu(ii) = this_nuc;
      contractions_cpu(ii) = this_cont;
      for (unsigned int k = 0; k < this_cont; k++)
        factor_ac_cpu(ii, k) = make_float2(fortran_vars.a_values(func, k), fortran_vars.c_values(func, k));
      ii++;
    }
  }
  factor_ac_gpu = factor_ac_cpu;
  nuc_gpu = nuc_cpu;
  contractions_gpu = contractions_cpu;

  /** Compute Functions **/

  function_values.resize(COALESCED_DIMENSION(number_of_points), group_functions.w);
  if (fortran_vars.do_forces || fortran_vars.gga) gradient_values.resize(COALESCED_DIMENSION(number_of_points), group_functions.w);
  if (fortran_vars.gga) hessian_values.resize(COALESCED_DIMENSION(number_of_points), group_functions.w * 2);

  dim3 threads(number_of_points);
  dim3 threadBlock(FUNCTIONS_BLOCK_SIZE);
  dim3 threadGrid = divUp(threads, threadBlock);

  //cout << "points: " << threads.x << " " << threadGrid.x << " " << threadBlock.x << endl;

  gpu_compute_functions<forces, gga><<<threadGrid, threadBlock>>>(points_position_gpu.data, number_of_points, contractions_gpu.data, factor_ac_gpu.data,
    nuc_gpu.data, function_values.data, gradient_values.data, hessian_values.data, group_functions);
  cudaAssertNoError("compute_functions");
}

template void PointGroup::compute_functions<true, false>(void);
template void PointGroup::compute_functions<true, true>(void);
template void PointGroup::compute_functions<false, false>(void);
template void PointGroup::compute_functions<false, true>(void);

/*******************************
 * Cube Weights
 *******************************/
void PointGroup::compute_weights(void)
{
  CudaMatrixFloat4 point_positions_gpu;
  CudaMatrixFloat4 atom_position_rm_gpu;
  {
    HostMatrixFloat4 points_positions_cpu(number_of_points, 1);
		uint i = 0;
		for (list<Point>::const_iterator p = points.begin(); p != points.end(); ++p, ++i) {
			points_positions_cpu(i) = make_float4(p->position.x, p->position.y, p->position.z, p->atom);
		}
    point_positions_gpu = points_positions_cpu;

    HostMatrixFloat4 atom_position_rm_cpu(fortran_vars.atoms, 1);
    for (uint i = 0; i < fortran_vars.atoms; i++) {
      double3 atom_pos = fortran_vars.atom_positions(i);
      atom_position_rm_cpu(i) = make_float4(atom_pos.x, atom_pos.y, atom_pos.z, fortran_vars.rm(i));
    }
    atom_position_rm_gpu = atom_position_rm_cpu;
	}

  CudaMatrixUInt nucleii_gpu(local2global_nuc);

  CudaMatrixFloat weights_gpu(number_of_points);
  dim3 threads(number_of_points);
  dim3 blockSize(WEIGHT_BLOCK_SIZE);
  dim3 gridSize = divUp(threads, blockSize);
  gpu_compute_weights<<<gridSize,blockSize>>>(number_of_points, point_positions_gpu.data, atom_position_rm_gpu.data,
                                              weights_gpu.data, nucleii_gpu.data, total_nucleii());
  cudaAssertNoError("compute_weights");

  #if REMOVE_ZEROS
  std::list<Point> nonzero_points;
  uint nonzero_number_of_points = 0;
  #endif

  uint ceros = 0;

  HostMatrixFloat weights_cpu(weights_gpu);
  uint i = 0;
  for (list<Point>::iterator p = points.begin(); p != points.end(); ++p, ++i) {
    p->weight *= weights_cpu(i);

    if (p->weight == 0.0) {
      ceros++;
    }
    #if REMOVE_ZEROS
    else {
      nonzero_points.push_back(*p);
      nonzero_number_of_points++;
    }
    #endif
  }

  //cout << "ceros: " << ceros << "/" << group.number_of_points << " (" << (ceros / (double)group.number_of_points) * 100 << "%)" << endl;

  #if REMOVE_ZEROS
  points = nonzero_points;
  number_of_points = nonzero_number_of_points;
  #endif
}
