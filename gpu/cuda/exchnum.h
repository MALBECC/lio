#ifndef __EXCHNUM_H__
#define __EXCHNUM_H__

void calc_energy(uint grid_type, uint npoints, uint Ndens, uint3 num_funcs, bool normalize, double& energy, double* cpu_rmm_output, double* cpu_forces_output,
								 bool compute_energy, bool update_rmm, bool compute_forces);

#endif
