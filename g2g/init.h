#ifndef __INIT_H__
#define __INIT_H__

#include "matrix.h"

#if USE_LIBINT
#include <unordered_map>
#include "libint/libintproxy.h"
using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
using shellpair_data_t = std::vector<std::vector<std::shared_ptr<libint2::ShellPair>>>;
#endif

using std::vector;

namespace G2G {

enum GridType { SMALL_GRID, MEDIUM_GRID, BIG_GRID };

enum GridSize {
  SMALL_GRID_SIZE = 50,
  MEDIUM_GRID_SIZE = 116,
  BIG_GRID_SIZE = 194
};

struct FortranVars {
  uint atoms, max_atoms, gaussians;
  bool normalize;
#if FULL_DOUBLE
  double normalization_factor;
#else
  float normalization_factor;
#endif
  uint s_funcs, p_funcs, d_funcs, spd_funcs, m;

  uint nco;
  bool OPEN;
  uint nunp;
  uint iexch;
  bool do_forces;
  bool gga, lda;
  GridType grid_type;
  GridSize grid_size;
  FortranMatrix<double> atom_positions_pointer;
  HostMatrix<double3> atom_positions;
  HostMatrix<uint> atom_types;
  HostMatrix<uint> atom_Z;
  HostMatrix<uint> shells, shells1, shells2;
  HostMatrix<double> rm, rm_base;
  HostMatrix<double> atom_atom_dists, nearest_neighbor_dists;
  FortranMatrix<uint> nucleii, contractions;
  FortranMatrix<double> a_values, c_values;
  FortranMatrix<double> rmm_input_ndens1, rmm_output;
  FortranMatrix<double> rmm_dens_a, rmm_dens_b, rmm_output_a, rmm_output_b;
  FortranMatrix<double> e, e1, e2, e3, wang, wang1, wang2, wang3;
  uint dim;
  uint nvirt;
  uint ncolr;

  // If using Becke partitioning.
  bool becke;
  HostMatrix<double> becke_atom_dens;
  HostMatrix<double> becke_atom_spin;

  /////////////////////////////////////
  // Agregado para integrar con Libxc
  bool use_libxc; // Si usa o no libxc
  uint ex_functional_id; // Identificador del funcional de intercambio (exchange)
  uint ec_functional_id; // Identificador del funcional de correlacion (correlation)
  /////////////////////////////////////

  // PBE0 factor
  double fexc;

  // Extern Functional Variables
  int nx_func, nc_func, nsr_id;
  int* func_id;
  double screen;
  double* func_coef;

  // LIBINT VARIABLES //
#if USE_LIBINT
  vector<libint2::Shell> obs; // Basis (in libint format)
  vector<int> shell2bf;       // first basis function
  vector<int> shell2atom;     // atom centre of shell
  uint center4Recalc;  // Method
  double* integrals;   // Integrals in memory
  shellpair_list_t obs_shellpair_list;  // shellpair list for precalculated Integral
  shellpair_data_t obs_shellpair_data;  // shellpair data for precalculated Integral
#endif
  /////////////////////
  
};

struct CDFTVars {
  bool do_chrg;
  bool do_spin;
  uint regions;
  uint max_nat;
  HostMatrix<double> Vc;
  HostMatrix<double> Vs;
  HostMatrix<uint>   natom;
  HostMatrix<uint>   atoms;
};

extern FortranVars fortran_vars;
extern CDFTVars    cdft_vars;

extern uint max_function_exponent;
extern double
    little_cube_size;  // largo de una arista de los cubos chiquitos [Angstrom]
extern uint min_points_per_cube;
extern double becke_cutoff;
extern bool assign_all_functions;
extern double sphere_radius;  // between 0 and 1!
extern bool remove_zero_weights;
extern bool energy_all_iterations;
extern double big_function_cutoff;
extern double free_global_memory;

extern bool timer_sum;
extern bool timer_single;
extern uint verbose;
}

#endif
