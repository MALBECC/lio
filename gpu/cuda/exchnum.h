#ifndef __EXCHNUM_H__
#define __EXCHNUM_H__

/*#define EXCHNUM_GRIDA_SIZE 116
#define EXCHNUM_GRIDB_SIZE 194*/

#define EXCHNUM_SMALL_GRID_SIZE		50
#define EXCHNUM_MEDIUM_GRID_SIZE	116
#define EXCHNUM_BIG_GRID_SIZE 		194

void calc_energy(const G2G::CudaMatrixFloat3& atom_positions, const G2G::HostMatrixUInt& types, uint grid_type,
								 uint npoits, double& energy_double, uint Ndens, unsigned int nco, uint3 num_funcs, const G2G::CudaMatrixUInt& nuc,
								 const G2G::CudaMatrixUInt& contractions, bool normalize, const G2G::CudaMatrixFloat2& factor_ac,
								 const G2G::CudaMatrixFloat& rmm, double* cpu_rmm_output, bool is_int3lu, const dim3& threads, const dim3& blockSize,
								 const dim3& gridSize3d);

#endif
