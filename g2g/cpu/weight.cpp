#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../common.h"
#include "../init.h"
#include "../matrix.h"
#include "../partition.h"
using namespace std;

namespace G2G {
template<class scalar_type>
void PointGroup<scalar_type>::compute_weights(void)
{
	for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it){
		uint atom = it->atom;
		double atom_weight;

		const double3& point_position = it->position;

		double P_total = 0.0;
		double P_atom = 0.0;

		for (uint j = 0; j < total_nucleii(); ++j) {
			double P_curr = 1.0;
			uint atom_j = local2global_nuc[j];
			const double3& pos_atom_j(fortran_vars.atom_positions(atom_j));
			double rm_atom_j = fortran_vars.rm(atom_j);

			for (uint k = 0; k < total_nucleii(); ++k) {
				uint atom_k = local2global_nuc[k];
				if (atom_k == atom_j) continue;
				const double3& pos_atom_k(fortran_vars.atom_positions(atom_k));
				double u = (length(point_position - pos_atom_j) - length(point_position - pos_atom_k)) / fortran_vars.atom_atom_dists(atom_j, atom_k);

				double x;
				x = rm_atom_j / fortran_vars.rm(atom_k);
        //cout << rm_atom_j << " " << fortran_vars.rm(atom_k) << endl;
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
    if (!has_nucleii(atom)) {
			P_atom = 1.0;
			uint atom_j = atom;
			const double3& pos_atom_j(fortran_vars.atom_positions(atom_j));
			double rm_atom_j = fortran_vars.rm(atom_j);

      for (uint k = 0; k < total_nucleii(); ++k) {
				uint atom_k = local2global_nuc[k];
				const double3& pos_atom_k(fortran_vars.atom_positions(atom_k));
				double u = (length(point_position - pos_atom_j) - length(point_position - pos_atom_k)) / fortran_vars.atom_atom_dists(atom_j, atom_k);

				double x;
				x = rm_atom_j / fortran_vars.rm(atom_k);
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
	}

    if (remove_zero_weights) {
        vector<Point> filteredPoints;
        for(vector<Point>::const_iterator it = points.begin(); it != points.end(); ++it){
            if(it->weight != 0.0) filteredPoints.push_back(*it);
        }
        points.swap(filteredPoints);
        number_of_points = points.size();
    }
}

template class PointGroup<double>;
template class PointGroup<float>;
}
