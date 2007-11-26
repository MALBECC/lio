#ifndef __EXCHNUM_H__
#define __EXCHNUM_H__

#define EXCHNUM_GRIDA_SIZE 116
#define EXCHNUM_GRIDB_SIZE 194

void calc_energy(const G2G::HostMatrixFloat3& atom_positions, const G2G::HostMatrixUInt& types, unsigned int grid_type,
								 const G2G::HostMatrixFloat3& point_positions, G2G::HostMatrixFloat& energy, const G2G::HostMatrixFloat& wang,
								 uint Iexch, bool Ndens, unsigned int nco, uint3 num_funcs, const G2G::HostMatrixUInt& nuc,
								 const G2G::HostMatrixUInt& contractions, bool normalize, const G2G::HostMatrixFloat& factor_a, const G2G::HostMatrixFloat& factor_c,
								 const G2G::HostMatrixFloat& rmm);

#endif
