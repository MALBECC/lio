#ifndef __EXCHNUM_H__
#define __EXCHNUM_H__

/*#define EXCHNUM_GRIDA_SIZE 116
#define EXCHNUM_GRIDB_SIZE 194*/

#define EXCHNUM_SMALL_GRID_SIZE		50
#define EXCHNUM_MEDIUM_GRID_SIZE	116
#define EXCHNUM_BIG_GRID_SIZE 		194

void calc_energy(const G2G::HostMatrixFloat3& atom_positions, const G2G::HostMatrixUInt& types, uint grid_type,
								 uint npoits, G2G::HostMatrixFloat& energy, uint Ndens, unsigned int nco, uint3 num_funcs, const G2G::HostMatrixUInt& nuc,
								 const G2G::HostMatrixUInt& contractions, bool normalize, const G2G::HostMatrixFloat& factor_a, const G2G::HostMatrixFloat& factor_c,
								 const G2G::HostMatrixFloat& rmm, double* cpu_rmm_output, bool is_int3lu, const dim3& threads, const dim3& blockSize,
								 const dim3& gridSize3d);

#endif
