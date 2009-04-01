/* includes */
#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "common.h"
#include "init.h"
#include "cuda/double.h"
#include "matrix.h"
#include "cuda/exchnum.h"
#include "cubes.h"
#include "functions.h"
#include "weight.h"
#include "timer.h"
using namespace std;
using namespace G2G;

/* local variables */


/* global variables */
list<LittleCube> final_cube;

void gpu_compute_cube_weights(LittleCube& cube);


/* methods */
void regenerate_cubes(void)
{
	cout << "<============ GPU Generate Cubes ============>" << endl;

	/* determina el exponente minimo para cada tipo de atomo */
	cout << "determining minimum exponents" << endl;	
	vector<double> min_exps(120, numeric_limits<double>::max());	// uno por elemento de la tabla periodica

	for (uint i = 0; i < fortran_vars.m; i++) {
		uint contractions = fortran_vars.contractions.get(i);
		uint nuc = fortran_vars.nucleii.get(i) - 1;
		uint nuc_type = fortran_vars.atom_types.get(nuc);
		for (uint j = 0; j < contractions; j++) {
			min_exps[nuc_type] = min(min_exps[nuc_type], fortran_vars.a_values.get(i, j));
		}
	}

	vector<double> min_exps_func(fortran_vars.m, numeric_limits<double>::max());	// uno por elemento de la tabla periodica
	for (uint i = 0; i < fortran_vars.m; i++) {
		uint contractions = fortran_vars.contractions.get(i);
		for (uint j = 0; j < contractions; j++) {
			min_exps_func[i] = min(min_exps_func[i], fortran_vars.a_values.get(i, j));
		}
	}
	
	/* permite encontrar el cubo conteniendo a todos los puntos */
	cout << "determining x0 and x1" << endl;
	double3 x0, x1;	
	for (uint atom = 0; atom < fortran_vars.atoms; atom++) {
		double3 atom_position(fortran_vars.atom_positions.get(atom));

		uint atom_type = fortran_vars.atom_types.get(atom);
		double max_radius = sqrt(max_function_exponent / min_exps[atom_type]);
		_DBG(cout << "tipo: " << atom_type << " " << min_exps[atom_type] << " radio: " << max_radius << endl);
		if (atom == 0) { x0 = x1 = atom_position; }
		else {
			if ((atom_position.x - max_radius) < x0.x) x0.x = (atom_position.x - max_radius);
			if ((atom_position.y - max_radius) < x0.y) x0.y = (atom_position.y - max_radius);
			if ((atom_position.z - max_radius) < x0.z) x0.z = (atom_position.z - max_radius);

			if ((atom_position.x + max_radius) > x1.x) x1.x = (atom_position.x + max_radius);
			if ((atom_position.y + max_radius) > x1.y) x1.y = (atom_position.y + max_radius);
			if ((atom_position.z + max_radius) > x1.z) x1.z = (atom_position.z + max_radius);
		}		
	}

	/* determina el tamanio del cubo grande en cubitos de tamanio "little_cube_size" */
	cout << "x0 " << x0.x << " " << x0.y << " " << x0.z << endl;
	cout << "x1 " << x1.x << " " << x1.y << " " << x1.z << endl;
	uint3 big_cube_size = ceil_uint3((x1 - x0) / little_cube_size);
	cout << "generating cube (" << big_cube_size.x << "," << big_cube_size.y << "," << big_cube_size.z << ")..." << endl;

	vector< vector < vector< LittleCube > > > big_cube(big_cube_size.x,
					vector< vector < LittleCube > >(big_cube_size.y,
									vector < LittleCube >(big_cube_size.z)));

	cout << "precomputing distances..." << endl;
	for (uint i = 0; i < fortran_vars.atoms; i++) {
		const double3& atom_i_position(fortran_vars.atom_positions.get(i));
		double nearest_neighbor_dist = numeric_limits<double>::max();
		for (uint j = 0; j < fortran_vars.atoms; j++) {
			const double3& atom_j_position(fortran_vars.atom_positions.get(j));
			double dist = (atom_i_position - atom_j_position).length();
			fortran_vars.atom_atom_dists.get(i, j) = dist;
			if (i != j) nearest_neighbor_dist = min(nearest_neighbor_dist, dist);
		}
		fortran_vars.nearest_neighbor_dists.get(i) = nearest_neighbor_dist;
	}
	
	cout << "computing points..." << endl;	
	uint puntos_totales = 0;
	
	Timer t_total;
	t_total.start();
	/* computa las posiciones de los puntos (y los guarda) */
	for (uint atom = 0; atom < fortran_vars.atoms; atom++) {		
		uint atom_shells = fortran_vars.shells.get(atom);
		
		const double3& atom_position(fortran_vars.atom_positions.get(atom));

		double t0 = M_PI / (atom_shells + 1);
		double rm = fortran_vars.rm.get(atom);

		for (uint shell = 0; shell < atom_shells; shell++) {			
			double t1 = t0 * (shell + 1);
			double x = cos(t1);
			double w = t0 * abs(sin(t1));
			double r1 = rm * (1.0 + x) / (1.0 - x);
      double wrad = w * (r1 * r1) * rm * 2.0 / ((1.0 - x) * (1.0 - x));
			
			for (uint point = 0; point < (uint)fortran_vars.grid_size; point++) {
				puntos_totales++;
				
				double3 rel_point_position(fortran_vars.e.get(point,0), fortran_vars.e.get(point,1), fortran_vars.e.get(point,2));
				double3 point_position = atom_position + rel_point_position * r1;
				bool inside_cube = ((x0.x <= point_position.x && point_position.x <= x1.x) &&
														(x0.y <= point_position.y && point_position.y <= x1.y) &&
														(x0.z <= point_position.z && point_position.z <= x1.z));

				if (inside_cube) {					
#if WEIGHT_GPU || WEIGHT_CUTOFFS
					double point_weight = wrad * fortran_vars.wang.get(point); // integration weight
#else
					double point_weight = compute_point_weight(point_position, wrad, atom, point);
					//if (point_weight == 0.0) continue;
#endif

          /* insert into corresponding cube */
					uint3 little_cube_coordinate = floor_uint3((point_position - x0) / little_cube_size);
					LittleCube& little_cube = big_cube[little_cube_coordinate.x][little_cube_coordinate.y][little_cube_coordinate.z];
				
					little_cube.points.push_back(Point(atom, shell, point, point_position, point_weight));
					little_cube.number_of_points++;
				}
			}
		}		
	}
	t_total.stop();
  
	cout << "Time: " << t_total << endl;
	cout << "Grilla original: " << puntos_totales << " funciones totales: " << puntos_totales * fortran_vars.m << endl;	
	
	puntos_totales = 0;
	uint funciones_totales = 0;

	final_cube.clear();

	cout << "filling cube..." << endl;
	t_total.start();
	for (uint i = 0; i < big_cube_size.x; i++) {
		for (uint j = 0; j < big_cube_size.y; j++) {
			for (uint k = 0; k < big_cube_size.z; k++) {
				LittleCube& little_cube = big_cube[i][j][k];

				double3 cube_coord_abs = x0 + make_uint3(i,j,k) * little_cube_size;
				if (little_cube.number_of_points < min_points_per_cube) { /*cout << "cubo vacio" << endl;*/ continue; }
				
				assign_significative_functions(little_cube, cube_coord_abs, min_exps_func);
				if (little_cube.functions.empty()) { /*cout << "cubo sin funciones" << endl;*/ continue; }

#if WEIGHT_GPU
        gpu_compute_cube_weights(little_cube);
#else
  #if WEIGHT_CUTOFFS
				assign_cube_weights(little_cube);
				if (little_cube.number_of_points < min_points_per_cube) { /*cout << "cubo vacio" << endl;*/ continue; }
  #endif
#endif
				final_cube.push_back(little_cube);

        // para hacer histogramas
#ifdef HISTOGRAM
        cout << "[" << fortran_vars.grid_type << "] cubo: " << little_cube.number_of_points << " puntos; " << little_cube.s_functions + little_cube.p_functions * 3 + little_cube.d_functions * 6 << " funciones" << endl;
#endif

				puntos_totales += little_cube.number_of_points;
				funciones_totales += little_cube.number_of_points * (little_cube.s_functions + little_cube.p_functions * 3 + little_cube.d_functions * 6);
			}
		}
	}
	t_total.stop();
	cout << "total: " << t_total << endl;
	cout << "Grilla final: " << puntos_totales << " funciones totales: " << funciones_totales << " puntos x cubo: " << (double)puntos_totales / final_cube.size() << " funciones x cubo: " << (double)funciones_totales / final_cube.size() << endl;
  cout << "Cubo final: " << final_cube.size() << " cubo original: " << big_cube_size.x * big_cube_size.y * big_cube_size.z << endl;
}
