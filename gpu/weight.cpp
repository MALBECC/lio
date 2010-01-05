#include <iostream>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "common.h"
#include "init.h"
#include "cuda/double.h"
#include "matrix.h"
#include "partition.h"
#include "weight.h"
using namespace std;
using namespace G2G;
 
void cpu_compute_group_weights(PointGroup& group)
{
	list<Point>::iterator it = group.points.begin();
	while (it != group.points.end()) {
		uint atom = it->atom;
		double atom_weight;

		const double3& point_position = it->position;

		double P_total = 0.0;
		double P_atom = 0.0;

		for (set<uint>::iterator atom_j_it = group.nucleii.begin(); atom_j_it != group.nucleii.end(); ++atom_j_it) {
			double P_curr = 1.0;
			uint atom_j = *atom_j_it;
			const double3& pos_atom_j(fortran_vars.atom_positions.get(atom_j));
			double rm_atom_j = fortran_vars.rm.get(atom_j);

			for (set<uint>::iterator atom_k_it = group.nucleii.begin(); atom_k_it != group.nucleii.end(); ++atom_k_it) {
				uint atom_k = *atom_k_it;
				if (atom_k == atom_j) continue;
				const double3& pos_atom_k(fortran_vars.atom_positions.get(atom_k));
				double u = ((point_position - pos_atom_j).length() - (point_position - pos_atom_k).length()) / fortran_vars.atom_atom_dists.get(atom_j, atom_k);

				double x;
				x = rm_atom_j / fortran_vars.rm.get(atom_k);
				x = (x - 1.0) / (x + 1.0);
				u += (x / (x * x - 1.0)) * (1.0 - u * u);

				u = 1.5 * u - 0.5 * (u * u * u);
				u = 1.5 * u - 0.5 * (u * u * u);
				u = 1.5 * u - 0.5 * (u * u * u);
				u = 0.5 * (1.0 - u);
				//cout << u << endl;

				P_curr *= u;
				if (P_curr == 0.0) { /*cout << "product" << endl;*/ break; }
			}

			if (atom_j == atom) {
				P_atom = P_curr;
				if (P_atom == 0.0) { /*cout << "curr" << endl;*/ break; }
			}

			P_total += P_curr;
		}

    // punto que no tiene a su propio atomo entre los vecinos
		if (group.nucleii.find(atom) == group.nucleii.end()) {
			P_atom = 1.0;
			uint atom_j = atom;
			const double3& pos_atom_j(fortran_vars.atom_positions.get(atom_j));
			double rm_atom_j = fortran_vars.rm.get(atom_j);
			
			for (set<uint>::iterator atom_k_it = group.nucleii.begin(); atom_k_it != group.nucleii.end(); ++atom_k_it) {
				uint atom_k = *atom_k_it;
				const double3& pos_atom_k(fortran_vars.atom_positions.get(atom_k));
				double u = ((point_position - pos_atom_j).length() - (point_position - pos_atom_k).length()) / fortran_vars.atom_atom_dists.get(atom_j, atom_k);

				double x;
				x = rm_atom_j / fortran_vars.rm.get(atom_k);
				x = (x - 1.0) / (x + 1.0);
				u += (x / (x * x - 1.0)) * (1.0 - u * u);

				u = 1.5 * u - 0.5 * (u * u * u);
				u = 1.5 * u - 0.5 * (u * u * u);
				u = 1.5 * u - 0.5 * (u * u * u);
				u = 0.5 * (1.0 - u);

				P_atom *= u;
				if (P_atom == 0.0) { /*cout << "product" << endl;*/ break; }
			}
			//cout << P_atom << " " << P_total << " " << P_atom / P_total << endl;
		}
		
		atom_weight = (P_total == 0.0 ? 0.0 : (P_atom / P_total));
		it->weight *= atom_weight;
		//cout << "peso " << P_atom << " " << P_total << " " << it->weight << endl;

		if (it->weight == 0.0) {
			it = group.points.erase(it);
			group.number_of_points--;
		}
		else ++it;
	}
}

