#include <iostream>
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
#include "weight.h"
using namespace std;
using namespace G2G;
 
double compute_point_weight(const double3& point_position, double wrad, uint atom, uint point)
{
	double P_total = 0.0;
	double P_atom = 1.0;

#if !BECKE
	double a = 0.64; // sacado de paper de DFT lineal
	if ((point_position - fortran_vars.atom_positions.get(atom)).length() < (0.5 * (1 - a) * fortran_vars.nearest_neighbor_dists.get(atom)))
		return 1.0;
#endif
					
	for (uint atom_j = 0; atom_j < fortran_vars.atoms; atom_j++) {
		double P_curr = 1.0;
		const double3& pos_atom_j(fortran_vars.atom_positions.get(atom_j));
		double rm_atom_j = fortran_vars.rm.get(atom_j);

		for (uint atom_k = 0; atom_k < fortran_vars.atoms; atom_k++) {
			if (atom_k == atom_j) continue;
			const double3& pos_atom_k(fortran_vars.atom_positions.get(atom_k));
			double u = ((point_position - pos_atom_j).length() - (point_position - pos_atom_k).length()) / fortran_vars.atom_atom_dists.get(atom_j, atom_k);

			double x;
			x = rm_atom_j / fortran_vars.rm.get(atom_k);
			x = (x - 1.0) / (x + 1.0);
			u += (x / (x * x - 1.0)) * (1.0 - u * u);

#if BECKE
			u = 1.5 * u - 0.5 * (u * u * u);
			u = 1.5 * u - 0.5 * (u * u * u);
			u = 1.5 * u - 0.5 * (u * u * u);
			u = 0.5 * (1.0 - u);
#else
			if (u <= -a) u = 1;
			else if (u >= a) u = 0;
			else {
				double ua = u / a;
				u = 2.1875 * (ua * (1 - ua) * (1 + ua)) + 1.3125 * pow(ua, 5) - 0.3125 * pow(ua, 7);
				u = 0.5 * (1.0 - u);
			}
#endif

			P_curr *= u;
			if (P_curr == 0.0) break;
		}

		if (atom_j == atom) {
			P_atom = P_curr;
			if (P_atom == 0.0) break;
		}
		P_total += P_curr;
	}

	double atom_weight = (P_total == 0.0 ? 0.0 : (P_atom / P_total));
	double integration_weight = wrad * fortran_vars.wang.get(point);
	double point_weight = atom_weight * integration_weight;
	
	return point_weight;
}

void assign_cube_weights(LittleCube& cube)
{
	list<Point>::iterator it = cube.points.begin();
	while (it != cube.points.end()) {
		uint atom = it->atom;

#if !BECKE
		if (cube.nucleii.find(atom) == cube.nucleii.end()) {
			it = cube.points.erase(it);
			cube.number_of_points--;
			continue;
		}
#endif

		const double3& point_position = it->position;
		double a = 0.64; // sacado de paper de DFT lineal

#if !BECKE
		if ((point_position - fortran_vars.atom_positions.get(atom)).length() < (0.5 * (1 - a) * fortran_vars.nearest_neighbor_dists.get(atom))) {
			it->weight = 1.0;
			++it;
			continue;
		}
#endif

		double P_total = 0.0;
		double P_atom = 0.0;

		for (set<uint>::iterator atom_j_it = cube.nucleii.begin(); atom_j_it != cube.nucleii.end(); ++atom_j_it) {
			double P_curr = 1.0;
			uint atom_j = *atom_j_it;
			const double3& pos_atom_j(fortran_vars.atom_positions.get(atom_j));
			double rm_atom_j = fortran_vars.rm.get(atom_j);

			for (set<uint>::iterator atom_k_it = cube.nucleii.begin(); atom_k_it != cube.nucleii.end(); ++atom_k_it) {
				uint atom_k = *atom_k_it;
				if (atom_k == atom_j) continue;
				const double3& pos_atom_k(fortran_vars.atom_positions.get(atom_k));
				double u = ((point_position - pos_atom_j).length() - (point_position - pos_atom_k).length()) / fortran_vars.atom_atom_dists.get(atom_j, atom_k);

				double x;
				x = rm_atom_j / fortran_vars.rm.get(atom_k);
				x = (x - 1.0) / (x + 1.0);
				u += (x / (x * x - 1.0)) * (1.0 - u * u);

#if BECKE
				u = 1.5 * u - 0.5 * (u * u * u);
				u = 1.5 * u - 0.5 * (u * u * u);
				u = 1.5 * u - 0.5 * (u * u * u);
				u = 0.5 * (1.0 - u);
#else
				if (u <= -a) u = 1;
				else if (u >= a) u = 0;
				else {
					double ua = u / a;
					u = 2.1875 * (ua * (1 - ua) * (1 + ua)) + 1.3125 * pow(ua, 5) - 0.3125 * pow(ua, 7);
					u = 0.5 * (1.0 - u);
				}
#endif	
				//cout << u << endl;

				P_curr *= u;
				if (P_curr == 0.0) break;
			}

			if (atom_j == atom) {
				P_atom = P_curr;
				if (P_atom == 0.0) break;
			}

			P_total += P_curr;
		}

		
		double atom_weight = (P_total == 0.0 ? 0.0 : (P_atom / P_total));
		it->weight *= atom_weight;
		//cout << "peso " << P_atom << " " << P_total << " " << it->weight << endl;

		/*if (it->weight == 0.0) {
			it = cube.points.erase(it);
			cube.number_of_points--;
		}
		else*/ ++it;
	}
}
