/* headers */
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <fenv.h>
#include <signal.h>
#include "common.h"
#include "cuda/cuda_extra.h"
#include "init.h"
#include "timer.h"
#include "partition.h"
#include "matrix.h"
using std::cout;
using std::endl;
using std::boolalpha;
using std::runtime_error;
using std::ifstream;
using std::string;
using namespace G2G;

/* external function prototypes */
template<bool, bool> void compute_functions(void);

/* internal function prototypes */
void read_options(void);

/* global variables */
namespace G2G {
	FortranVars fortran_vars;
}

/* methods */
extern "C" void g2g_init_(void)
{
  cout << "<====== Initializing G2G ======>";

  #if !CPU_KERNELS
  cuInit(0);
  int devnum = -1;
  cudaDeviceProp devprop;
  if (cudaGetDevice(&devnum) != cudaSuccess) throw runtime_error("Could not get device number!");
  if (cudaGetDeviceProperties(&devprop, devnum) != cudaSuccess) throw runtime_error("Could not get device propierties!");
  cout << "GPU Device used: " << devprop.name << endl;
  cout << "Kernels: gpu" << endl;
  #else
  cout << "Kernels: cpu" << endl;
  #endif

  cout.precision(10);
}

namespace G2G {
void gpu_set_variables(void);
template<class T> void gpu_set_atom_positions(const HostMatrix<T>& m);
}

extern "C" void g2g_parameter_init_(const unsigned int& norm, const unsigned int& natom, const unsigned int& max_atoms, const unsigned int& ngaussians, double* r, double* Rm, const unsigned int* Iz, const unsigned int* Nr,
										 const unsigned int* Nr2, unsigned int* Nuc, const unsigned int& M, unsigned int* ncont, const unsigned int* nshell, double* c,
					           double* a, double* RMM, const unsigned int& M18, const unsigned int& M5, const unsigned int& nco, const unsigned int& nopt,
                     const unsigned int& Iexch, double* e, double* e2, double* e3, double* wang, double* wang2, double* wang3)
{
	printf("<======= GPU Code Initialization ========>\n");
	fortran_vars.atoms = natom;
        fortran_vars.max_atoms = max_atoms;
        fortran_vars.gaussians = ngaussians;

	cout << "atoms: " << fortran_vars.atoms << endl;
        cout << "max atoms: " << fortran_vars.max_atoms << endl;
        cout << "number of gaussians: " << fortran_vars.gaussians << endl;
	
	fortran_vars.do_forces = (nopt == 2);
	cout << "do_forces: " << boolalpha << fortran_vars.do_forces << endl;

	fortran_vars.normalize = norm;
	fortran_vars.normalization_factor = (fortran_vars.normalize ? (1.0/sqrt(3)) : 1.0);

  #ifdef _DEBUG
  // trap floating point exceptions on debug
  signal(SIGFPE, SIG_DFL);
  feenableexcept(FE_INVALID);
  #endif
	
	fortran_vars.s_funcs = nshell[0];
	fortran_vars.p_funcs = nshell[1] / 3;
	fortran_vars.d_funcs = nshell[2] / 6;
  cout << "s: " << fortran_vars.s_funcs  << " p: " << fortran_vars.p_funcs << " d: " << fortran_vars.d_funcs << endl;
	fortran_vars.spd_funcs = fortran_vars.s_funcs + fortran_vars.p_funcs + fortran_vars.d_funcs;
	fortran_vars.m = M;	
	fortran_vars.nco = nco;
  cout << "m: " << fortran_vars.m  << " nco: " << fortran_vars.nco << endl;
		
	fortran_vars.iexch = Iexch;
  if (Iexch == 4 || Iexch == 5) cout << "***** WARNING ***** : Iexch 4 y 5 no andan bien todavia" << endl;
  fortran_vars.lda = (Iexch <= 3);
  fortran_vars.gga = !fortran_vars.lda;
  assert(0 < Iexch && Iexch <= 9);
	
	fortran_vars.atom_positions_pointer = FortranMatrix<double>(r, fortran_vars.atoms, 3, max_atoms);
	fortran_vars.atom_types.resize(fortran_vars.atoms);
	for (uint i = 0; i < fortran_vars.atoms; i++) { fortran_vars.atom_types(i) = Iz[i] - 1; }
	
	fortran_vars.shells1.resize(fortran_vars.atoms);
	fortran_vars.shells2.resize(fortran_vars.atoms);
	fortran_vars.rm.resize(fortran_vars.atoms);

	/* ignore the 0th element on these */
	for (uint i = 0; i < fortran_vars.atoms; i++) { fortran_vars.shells1(i) = Nr[Iz[i]];  }
	for (uint i = 0; i < fortran_vars.atoms; i++) { fortran_vars.shells2(i) = Nr2[Iz[i]]; }		
	for (uint i = 0; i < fortran_vars.atoms; i++) { fortran_vars.rm(i) = Rm[Iz[i]]; }
  
	fortran_vars.nucleii = FortranMatrix<uint>(Nuc, fortran_vars.m, 1, 1);	
	fortran_vars.contractions = FortranMatrix<uint>(ncont, fortran_vars.m, 1, 1);
  for (uint i = 0; i < fortran_vars.m; i++) {
    if ((fortran_vars.contractions(i) - 1) > MAX_CONTRACTIONS)  throw runtime_error("Maximum functions per contraction reached!");
}
	fortran_vars.a_values = FortranMatrix<double>(a, fortran_vars.m, MAX_CONTRACTIONS, ngaussians);
	fortran_vars.c_values = FortranMatrix<double>(c, fortran_vars.m, MAX_CONTRACTIONS, ngaussians);
	fortran_vars.rmm_input = FortranMatrix<double>(RMM + (M18 - 1), fortran_vars.m, fortran_vars.nco, fortran_vars.m);
  fortran_vars.rmm_input_ndens1 = FortranMatrix<double>(RMM, fortran_vars.m, fortran_vars.m, fortran_vars.m);
	fortran_vars.rmm_output = FortranMatrix<double>(RMM + (M5 - 1), (fortran_vars.m * (fortran_vars.m + 1)) / 2);

	fortran_vars.e1 = FortranMatrix<double>(e, SMALL_GRID_SIZE, 3, SMALL_GRID_SIZE);
	fortran_vars.e2 = FortranMatrix<double>(e2, MEDIUM_GRID_SIZE, 3, MEDIUM_GRID_SIZE);
	fortran_vars.e3 = FortranMatrix<double>(e3, BIG_GRID_SIZE, 3, BIG_GRID_SIZE);
	fortran_vars.wang1 = FortranMatrix<double>(wang, SMALL_GRID_SIZE, 1, SMALL_GRID_SIZE);
	fortran_vars.wang2 = FortranMatrix<double>(wang2, MEDIUM_GRID_SIZE, 1, MEDIUM_GRID_SIZE);
	fortran_vars.wang3 = FortranMatrix<double>(wang3, BIG_GRID_SIZE, 1, BIG_GRID_SIZE);

	fortran_vars.atom_atom_dists = HostMatrix<double>(fortran_vars.atoms, fortran_vars.atoms);
	fortran_vars.nearest_neighbor_dists = HostMatrix<double>(fortran_vars.atoms);
  
#if !CPU_KERNELS
  G2G::gpu_set_variables();
#endif

	read_options();
}

extern "C" void g2g_deinit_(void) {
  cout << "<====== Deinitializing G2G ======>" << endl;
  partition.clear();
}

void compute_new_grid(const unsigned int grid_type) {
	switch(grid_type) {
		case 0:
			fortran_vars.grid_type = SMALL_GRID; fortran_vars.grid_size = SMALL_GRID_SIZE;
			fortran_vars.e = fortran_vars.e1; fortran_vars.wang = fortran_vars.wang1;
			fortran_vars.shells = fortran_vars.shells1;
		break;
		case 1:
			fortran_vars.grid_type = MEDIUM_GRID; fortran_vars.grid_size = MEDIUM_GRID_SIZE;
			fortran_vars.e = fortran_vars.e2; fortran_vars.wang = fortran_vars.wang2;
			fortran_vars.shells = fortran_vars.shells2;
		break;
		case 2:
			fortran_vars.grid_type = BIG_GRID; fortran_vars.grid_size = BIG_GRID_SIZE;
			fortran_vars.e = fortran_vars.e3; fortran_vars.wang = fortran_vars.wang3;
			fortran_vars.shells = fortran_vars.shells2;
		break;
	}	

  Timer t_grilla;
  t_grilla.start_and_sync();
  partition.regenerate();
  t_grilla.stop_and_sync();
  //cout << "timer grilla: " << t_grilla << endl;

#if CPU_KERNELS && !CPU_RECOMPUTE
  /** compute functions **/
  //if (fortran_vars.do_forces) cout << "<===== computing all functions [forces] =======>" << endl;
  //else cout << "<===== computing all functions =======>" << endl;

  partition.compute_functions(fortran_vars.do_forces, fortran_vars.gga);
#endif
}

extern "C" void g2g_reload_atom_positions_(const unsigned int& grid_type) {
//	cout  << "<======= GPU Reload Atom Positions (" << grid_type << ")========>" << endl;

	HostMatrixFloat3 atom_positions(fortran_vars.atoms);	// gpu version (float3)
	fortran_vars.atom_positions.resize(fortran_vars.atoms);	// cpu version (double3)
	for (uint i = 0; i < fortran_vars.atoms; i++) {
		double3 pos = make_double3(fortran_vars.atom_positions_pointer(i, 0), fortran_vars.atom_positions_pointer(i, 1), fortran_vars.atom_positions_pointer(i, 2));
		fortran_vars.atom_positions(i) = pos;
		atom_positions(i) = make_float3(pos.x, pos.y, pos.z);
    //cout << atom_positions(i) << endl;
	}
  
#if !CPU_KERNELS
#if FULL_DOUBLE
  G2G::gpu_set_atom_positions(fortran_vars.atom_positions);
#else
  G2G::gpu_set_atom_positions(atom_positions);
#endif
#endif

	compute_new_grid(grid_type);
}

extern "C" void g2g_new_grid_(const unsigned int& grid_type) {
//	cout << "<======= GPU New Grid (" << grid_type << ")========>" << endl;
	if (grid_type == (uint)fortran_vars.grid_type)
                                                             ;
//		cout << "not loading, same grid as loaded" << endl;
	else
		compute_new_grid(grid_type);
}

template<bool compute_rmm, bool lda, bool compute_forces> void g2g_iteration(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr) 
{
  Timers timers;
  timers.total.start();

  partition.solve(timers, compute_rmm, lda, compute_forces, compute_energy, fort_energy_ptr, fort_forces_ptr);

  timers.total.stop();
  cout << timers << endl;  
}

extern "C" void g2g_solve_groups_(const uint& computation_type, double* fort_energy_ptr, double* fort_forces_ptr)
{
 /*	cout << "<================ iteracion [";
	switch(computation_type) {
    case COMPUTE_RMM: cout << "rmm"; break;
		case COMPUTE_ENERGY_ONLY: cout << "energia"; break;
		case COMPUTE_ENERGY_FORCE: cout << "energia+fuerzas"; break;
		case COMPUTE_FORCE_ONLY: cout << "fuerzas"; break;
	}
	cout << "] ==========>" << endl;
*/
  bool compute_energy = (computation_type == COMPUTE_ENERGY_ONLY || computation_type == COMPUTE_ENERGY_FORCE);
  bool compute_forces = (computation_type == COMPUTE_FORCE_ONLY || computation_type == COMPUTE_ENERGY_FORCE);
  bool compute_rmm = (computation_type == COMPUTE_RMM);

  if (energy_all_iterations) compute_energy = true;

  if (compute_rmm) {
    if (fortran_vars.lda) {
      if (compute_forces) g2g_iteration<true, true, true>(compute_energy, fort_energy_ptr, fort_forces_ptr);
      else g2g_iteration<true, true, false>(compute_energy, fort_energy_ptr, fort_forces_ptr);
    }
    else {
      if (compute_forces) g2g_iteration<true, false, true>(compute_energy, fort_energy_ptr, fort_forces_ptr);
      else g2g_iteration<true, false, false>(compute_energy, fort_energy_ptr, fort_forces_ptr);
    }
  }
  else {
    if (fortran_vars.lda) {
      if (compute_forces) g2g_iteration<false, true, true>(compute_energy, fort_energy_ptr, fort_forces_ptr);
      else g2g_iteration<false, true, false>(compute_energy, fort_energy_ptr, fort_forces_ptr);
    }
    else {
      if (compute_forces) g2g_iteration<false, false, true>(compute_energy, fort_energy_ptr, fort_forces_ptr);
      else g2g_iteration<false, false, false>(compute_energy, fort_energy_ptr, fort_forces_ptr);
    }
  }
  if (compute_energy) cout << "XC energy: " << *fort_energy_ptr << endl;
}

/* general options */
namespace G2G {
	uint max_function_exponent = 10;
	double little_cube_size = 8.0;
	uint min_points_per_cube = 1;
	double becke_cutoff = 1e-7;
  bool assign_all_functions = false;
  double sphere_radius = 0.6;
  bool remove_zero_weights = true;
  bool energy_all_iterations = false;
  double big_function_cutoff = 1;
  double free_global_memory = 0.0;
}

void read_options(void) {
	cout << "<====== read_options ========>" << endl;
	ifstream f("gpu_options");
	if (!f) { cout << "No \"gpu_options\" file: using defaults" << endl; return; }
	
	string option;
	while (f >> option) {
		cout << option << " ";

		if (option == "max_function_exponent")
			{ f >> max_function_exponent; cout << max_function_exponent; }
		else if (option == "little_cube_size")
			{ f >> little_cube_size; cout << little_cube_size; }
		else if (option == "min_points_per_cube")
			{ f >> min_points_per_cube; cout << min_points_per_cube; }
		else if (option == "becke_cutoff")
			{ f >> becke_cutoff; cout << becke_cutoff; }
    else if (option == "assign_all_functions")
      { f >> assign_all_functions; cout << assign_all_functions; }
    else if (option == "sphere_radius")
      { f >> sphere_radius; cout << sphere_radius; }
    else if (option == "remove_zero_weights")
      { f >> remove_zero_weights; cout << remove_zero_weights; }
    else if (option == "energy_all_iterations")
      { f >> energy_all_iterations; cout << energy_all_iterations; }
    else if (option == "big_function_cutoff")
      { f >> big_function_cutoff; cout << big_function_cutoff; }
    else if (option == "free_global_memory")
      { f >> free_global_memory; cout << free_global_memory; }
		else throw runtime_error(string("Invalid option: ") + option);

		cout << endl;
		if (!f) throw runtime_error(string("Error reading gpu options file (last option: ") + option + string(")"));
	}
}
