#ifndef __EXCHNUM_H__
#define __EXCHNUM_H__

/*#define EXCHNUM_GRIDA_SIZE 116
#define EXCHNUM_GRIDB_SIZE 194*/

#define EXCHNUM_SMALL_GRID_SIZE		50
#define EXCHNUM_MEDIUM_GRID_SIZE	116
#define EXCHNUM_BIG_GRID_SIZE 		194

void calc_energy(uint grid_type, uint npoints, uint Ndens, uint3 num_funcs, bool normalize, double& energy, double* cpu_rmm_output, double* cpu_forces_output,
								 bool compute_energy, bool update_rmm, bool compute_forces);

#endif
