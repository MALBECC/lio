/* includes */
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cassert>

#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"

using namespace std;
using namespace G2G;

/************************************************************
 * Construct partition
 ************************************************************/

// Sorting the cubes in increasing order of size in bytes in GPU.
template <typename T>
bool comparison_by_size(const T& a, const T& b) {
  return a.size_in_gpu() < b.size_in_gpu();
}

template <typename T>
void sortBySize(std::vector<T>& input) {
  sort(input.begin(), input.end(), comparison_by_size<T>);
}

void load_pools(const vector<int>& elements, const vector<vector<int> >& work,
                vector<int>& pool_sizes) {
  pool_sizes.clear();
  for (uint i = 0; i < work.size(); i++) {
    int largest_pool = 0;
    for (uint j = 0; j < work[i].size(); j++) {
      largest_pool = max(largest_pool, elements[work[i][j]]);
    }
    pool_sizes.push_back(largest_pool);
  }
}

template <typename T>
long long total_costs(const vector<T*>& elements) {
  long long res = 0;
  for (uint i = 0; i < elements.size(); i++) res += elements[i]->cost();
  return res;
}

// Computes Becke a_ij factors for Becke partitioning.
int compute_becke_a(HostMatrix<double>& becke_a_in) {
  HostMatrix<double> cov_r = HostMatrix<double>(97);

  // The following are the average covalent radii according to
  // Cambridge Structural Database. All units are in Bohrs.
  // Available data is up to Cm (Z=96), if you want more
  // go measure them.
  cov_r( 0) = 0.0000; // For ghost atoms.
  /* USING COVALENT RADII
  cov_r( 1) = 0.5858; cov_r( 2) = 0.5291; cov_r( 3) = 2.4189;
  cov_r( 4) = 1.8141; cov_r( 5) = 1.5874; cov_r( 6) = 1.3795;
  cov_r( 7) = 1.3417; cov_r( 8) = 1.2472; cov_r( 9) = 1.0771;
  cov_r(10) = 1.0960; cov_r(11) = 3.1370; cov_r(12) = 2.6645;
  cov_r(13) = 2.2866; cov_r(14) = 2.0976; cov_r(15) = 2.0220;
  cov_r(16) = 1.9842; cov_r(17) = 1.9275; cov_r(18) = 2.0031;
  cov_r(19) = 3.8362; cov_r(20) = 3.3259; cov_r(21) = 3.2125;
  cov_r(22) = 3.0236; cov_r(23) = 2.8913; cov_r(24) = 2.6267;
  cov_r(25) = 2.8346; cov_r(26) = 2.6834; cov_r(27) = 2.6078;
  cov_r(28) = 2.3433; cov_r(29) = 2.4944; cov_r(30) = 2.3055;
  cov_r(31) = 2.3055; cov_r(32) = 2.2677; cov_r(33) = 2.2488;
  cov_r(34) = 2.2677; cov_r(35) = 2.2677; cov_r(36) = 2.1921;
  cov_r(37) = 4.1574; cov_r(38) = 3.6850; cov_r(39) = 3.5905;
  cov_r(40) = 3.3070; cov_r(41) = 3.0992; cov_r(42) = 2.9102;
  cov_r(43) = 2.7779; cov_r(44) = 2.7590; cov_r(45) = 2.6834;
  cov_r(46) = 2.6267; cov_r(47) = 2.7401; cov_r(48) = 2.7212;
  cov_r(49) = 2.6834; cov_r(50) = 2.6267; cov_r(51) = 2.6267;
  cov_r(52) = 2.6078; cov_r(53) = 2.6267; cov_r(54) = 2.6456;
  cov_r(55) = 4.6109; cov_r(56) = 4.0629; cov_r(57) = 3.9117;
  cov_r(58) = 3.8550; cov_r(59) = 3.8362; cov_r(60) = 3.7984;
  cov_r(61) = 3.7606; cov_r(62) = 3.7417; cov_r(63) = 3.7417;
  cov_r(64) = 3.7039; cov_r(65) = 3.6661; cov_r(66) = 3.6283;
  cov_r(67) = 3.6283; cov_r(68) = 3.5716; cov_r(69) = 3.5905;
  cov_r(70) = 3.5338; cov_r(71) = 3.3070; cov_r(72) = 3.5338;
  cov_r(73) = 3.2125; cov_r(74) = 3.0614; cov_r(75) = 2.8535;
  cov_r(76) = 2.7212; cov_r(77) = 2.6645; cov_r(78) = 2.5700;
  cov_r(79) = 2.5700; cov_r(80) = 2.4944; cov_r(81) = 2.7401;
  cov_r(82) = 2.7590; cov_r(83) = 2.7968; cov_r(84) = 2.6456;
  cov_r(85) = 2.8346; cov_r(86) = 2.8346; cov_r(87) = 4.9133;
  cov_r(88) = 4.1763; cov_r(89) = 4.0629; cov_r(90) = 3.8928;
  cov_r(91) = 3.7795; cov_r(92) = 3.7039; cov_r(93) = 3.5905;
  cov_r(94) = 3.5338; cov_r(95) = 3.4015; cov_r(96) = 3.1936;
  */
  // USING SLATER-BRAGG ATOMIC RADII
  // H is modified to be 0.35A instead of 0.25A
  cov_r( 1) = 0.6614; cov_r( 2) = 2.2677; cov_r( 3) = 2.7401;
  cov_r( 4) = 1.9842; cov_r( 5) = 1.6063; cov_r( 6) = 1.3228;
  cov_r( 7) = 1.2283; cov_r( 8) = 1.1338; cov_r( 9) = 0.9449;
  cov_r(10) = 3.0236; cov_r(11) = 3.4015; cov_r(12) = 2.8346;
  cov_r(13) = 2.3622; cov_r(14) = 2.0787; cov_r(15) = 1.8897;
  cov_r(16) = 1.8897; cov_r(17) = 1.8897; cov_r(18) = 1.3417;
  cov_r(19) = 4.1574; cov_r(20) = 3.4015; cov_r(21) = 3.0236;
  cov_r(22) = 2.6456; cov_r(23) = 2.5511; cov_r(24) = 2.6456;
  cov_r(25) = 2.6456; cov_r(26) = 2.6456; cov_r(27) = 2.5511;
  cov_r(28) = 2.5511; cov_r(29) = 2.5511; cov_r(30) = 2.5511;
  cov_r(31) = 2.4566; cov_r(32) = 2.3622; cov_r(33) = 2.1732;
  cov_r(34) = 2.1732; cov_r(35) = 2.1732; cov_r(36) = 1.6630;
  cov_r(37) = 4.4409; cov_r(38) = 3.7795; cov_r(39) = 3.4015;
  cov_r(40) = 2.9291; cov_r(41) = 2.7401; cov_r(42) = 2.7401;
  cov_r(43) = 2.5511; cov_r(44) = 2.4566; cov_r(45) = 2.5511;
  cov_r(46) = 2.6456; cov_r(47) = 3.0236; cov_r(48) = 2.9291;
  cov_r(49) = 2.9291; cov_r(50) = 2.7401; cov_r(51) = 2.7401;
  cov_r(52) = 2.6456; cov_r(53) = 2.6456; cov_r(54) = 2.0409;
  cov_r(55) = 4.9133; cov_r(56) = 4.0629; cov_r(57) = 3.6850;
  cov_r(58) = 3.4960; cov_r(59) = 3.4960; cov_r(60) = 3.4960;
  cov_r(61) = 3.4960; cov_r(62) = 3.4960; cov_r(63) = 3.4960;
  cov_r(64) = 3.4015; cov_r(65) = 3.3070; cov_r(66) = 3.3070;
  cov_r(67) = 3.3070; cov_r(68) = 3.3070; cov_r(69) = 3.3070;
  cov_r(70) = 3.3070; cov_r(71) = 3.3070; cov_r(72) = 2.9291;
  cov_r(73) = 2.7401; cov_r(74) = 2.5511; cov_r(75) = 2.5511;
  cov_r(76) = 2.4566; cov_r(77) = 2.5511; cov_r(78) = 2.5511;
  cov_r(79) = 2.5511; cov_r(80) = 2.8346; cov_r(81) = 3.5905;
  cov_r(82) = 3.4015; cov_r(83) = 3.0236; cov_r(84) = 3.5905;
  cov_r(85) = 2.4000; cov_r(86) = 2.2677; cov_r(87) = 6.5763;
  cov_r(88) = 4.0629; cov_r(89) = 3.6850; cov_r(90) = 3.4015;
  cov_r(91) = 3.4015; cov_r(92) = 3.3070; cov_r(93) = 3.3070;
  cov_r(94) = 3.3070; cov_r(95) = 3.3070; cov_r(96) = 3.1936;

  for (uint iatom = 0; iatom < fortran_vars.atoms; iatom++) {
    for (uint jatom = 0; jatom < iatom; jatom++) {
      // Uij = (Ri - Rj) / (Ri + Rj)
      // Aij = Uij / (Uij^2 - 1)
      if (fortran_vars.atom_types(iatom) == fortran_vars.atom_types(jatom)) {
        becke_a_in(iatom, jatom) = 0.0;
      } else {
        double r_i = cov_r(fortran_vars.atom_types(iatom)+1);
        double r_j = cov_r(fortran_vars.atom_types(jatom)+1);
        double u_ij = (r_i - r_j) / (r_i + r_j);
        u_ij = u_ij / (u_ij * u_ij - 1.0);        

        if (u_ij > 0.5) {
          u_ij = 0.5;
        } else if (u_ij < -0.5) {
          u_ij = -0.5;
        }

        becke_a_in(iatom, jatom) = u_ij;
      } 
      becke_a_in(jatom, iatom) = -becke_a_in(iatom,jatom);
    }
    becke_a_in(iatom, iatom) = 0.0;

  }
  return 0;
} // compute_becke_a

int compute_becke_w(HostMatrix<double>& becke_w, double3 p_pos,
                    HostMatrix<double>& becke_a){
  // P_i is the so-called "atomic cell function". It is 
  // initialized with 1.0 in order to use in the final products.
  HostMatrix<double> P_i;
  P_i.resize(fortran_vars.atoms);
  for (uint iatom = 0; iatom < fortran_vars.atoms; iatom++) {
    P_i(iatom) = 1.0f;
  }

  double P_total = 0.0f;
  for (uint iatom = 0; iatom < fortran_vars.atoms; iatom++) {
    for (uint jatom = 0; jatom < fortran_vars.atoms; jatom++) {
      if (!(iatom == jatom)) {
        double3 i_pos = fortran_vars.atom_positions(iatom);
        double3 j_pos = fortran_vars.atom_positions(jatom);
        double dist_ip = sqrt((i_pos.x - p_pos.x) * (i_pos.x - p_pos.x) +
                              (i_pos.y - p_pos.y) * (i_pos.y - p_pos.y) +
                              (i_pos.z - p_pos.z) * (i_pos.z - p_pos.z) );

        double dist_jp = sqrt((j_pos.x - p_pos.x) * (j_pos.x - p_pos.x) +
                              (j_pos.y - p_pos.y) * (j_pos.y - p_pos.y) +
                              (j_pos.z - p_pos.z) * (j_pos.z - p_pos.z) );

        double dist_ij = sqrt((i_pos.x - j_pos.x) * (i_pos.x - j_pos.x) +
                              (i_pos.y - j_pos.y) * (i_pos.y - j_pos.y) +
                              (i_pos.z - j_pos.z) * (i_pos.z - j_pos.z) );

        double mu_ij  = (dist_ip - dist_jp) / dist_ij;

        // The first f_beck is nu_ij, then it is calculated recursively
        // up to order 3. Finally, f_beck becomes S_ij
        double f_beck = mu_ij + becke_a(iatom, jatom) * (1.0f - mu_ij * mu_ij);
        for (int k = 0; k < 3; k++) {
          f_beck = 0.5f * f_beck * (3.0f - f_beck * f_beck);
        }

        f_beck = 0.5f * (1.0f - f_beck);
        P_i(iatom) = P_i(iatom) * f_beck;
      }
    }
    P_total = P_total + P_i(iatom);
  }
  becke_w(fortran_vars.atoms-1) = 1.0f;
  for (uint iatom = 0; iatom < (fortran_vars.atoms-1); iatom++) {
    becke_w(iatom) = P_i(iatom) / P_total;
    if (becke_w(iatom) < (float) 1.0E-6 ) becke_w(iatom) = 0.0f;
    becke_w(fortran_vars.atoms-1) -= becke_w(iatom);
  }


  return 0;
} // compute_becke_w

int split_bins(const vector<pair<long long, int> >& costs,
               vector<vector<int> >& workloads, long long capacity) {
  // Bin Packing heuristic
  workloads.clear();
  for (uint i = 0; i < costs.size(); i++) {
    int next_bin = -1;
    for (uint j = 0; j < workloads.size(); j++) {
      long long slack = capacity;
      for (uint k = 0; k < workloads[j].size(); k++) {
        slack -= costs[workloads[j][k]].first;
      }
      if (slack >= costs[i].first && next_bin == -1) {
        next_bin = j;
        break;
      }
    }
    if (next_bin == -1) {
      if (capacity < costs[i].second) {
        return INT_MAX;
      }
      next_bin = workloads.size();
      workloads.push_back(vector<int>());
    }
    workloads[next_bin].push_back(i);
  }
  return workloads.size();
}

void Partition::compute_work_partition() {
  if (G2G::cpu_threads == 0) return;
  vector<pair<long long, int> > costs;
  for (uint i = 0; i < cubes.size(); i++)
    if (!cubes[i]->is_big_group())
      costs.push_back(make_pair(cubes[i]->cost(), i));

  const uint ncubes = cubes.size();
  for (uint i = 0; i < spheres.size(); i++)
    if (!spheres[i]->is_big_group())
      costs.push_back(make_pair(spheres[i]->cost(), ncubes + i));

  if (costs.empty()) return;

  sort(costs.begin(), costs.end());
  reverse(costs.begin(), costs.end());

  long long min_cost = costs.front().second - 1,
            max_cost = total_costs(cubes) + total_costs(spheres) + 1;

  while (max_cost - min_cost > 1) {
    long long candidate = min_cost + (max_cost - min_cost) / 2;

    vector<vector<int> > workloads;
    int bins = split_bins(costs, workloads, candidate);
    if (bins <= G2G::cpu_threads) {
      max_cost = candidate;
    } else {
      min_cost = candidate;
    }
  }

  split_bins(costs, work, max_cost);
  for (uint i = 0; i < work.size(); i++) sort(work[i].begin(), work[i].end());

  double maxp = 0, minp = total_costs(cubes) + total_costs(spheres) + 1;
  for (uint i = 0; i < work.size(); i++) {
    long long total = 0;
    for (uint j = 0; j < work[i].size(); j++) {
      long long c = costs[work[i][j]].first;
      work[i][j] = costs[work[i][j]].second;
      total += c;
    }
    if (minp > total) minp = total;
    if (maxp < total) maxp = total;
    if (verbose > 4) printf("  Partition %d: %lld\n", i, total);
  }
  if (verbose > 4) printf("  MAX/MIN ratio: %lf\n", maxp / minp);
}

int getintenv(const char* str, int default_value) {
  char* v = getenv(str);
  if (v == NULL) return default_value;
  int ret = strtol(v, NULL, 10);
  return ret;
}

void diagnostic() {
  printf("  Threads OMP: %d - Threads CPU: %d - Threads GPU: %d\n",
         omp_get_max_threads(), G2G::cpu_threads, G2G::gpu_threads);
  printf("  Small cube correction: %d - Separation points: %d\n",
         MINCOST, SPLITPOINTS);
}

template <class T>
bool is_big_group(const T& points) {
  assert(G2G::cpu_threads > 0 || G2G::gpu_threads > 0);
  if (G2G::cpu_threads == 0) return true;
  if (G2G::gpu_threads == 0) return false;
  return (points.size() > G2G::SPLITPOINTS);
}

struct Sorter {
  template <class T>
  bool operator()(const T l, const T r) {
    return (l->cost() < r->cost());
  }
};

/* methods */
void Partition::regenerate(void) {
  Timer tweights;
  // Determina el exponente minimo para cada tipo de atomo.
  // uno por elemento de la tabla periodica.
  vector<double> min_exps(120, numeric_limits<double>::max());
  for (uint i = 0; i < fortran_vars.m; i++) {
    uint contractions = fortran_vars.contractions(i);
    uint nuc = fortran_vars.nucleii(i) - 1;
    uint nuc_type = fortran_vars.atom_types(nuc);
    for (uint j = 0; j < contractions; j++) {
      min_exps[nuc_type] = min(min_exps[nuc_type], fortran_vars.a_values(i, j));
    }
  }

  // Un exponente y un coeficiente por funcion.
  vector<double> min_exps_func(fortran_vars.m, numeric_limits<double>::max());
  vector<double> min_coeff_func(fortran_vars.m);
  for (uint i = 0; i < fortran_vars.m; i++) {
    uint contractions = fortran_vars.contractions(i);
    for (uint j = 0; j < contractions; j++) {
      if (fortran_vars.a_values(i, j) < min_exps_func[i]) {
        min_exps_func[i] = fortran_vars.a_values(i, j);
        min_coeff_func[i] = fortran_vars.c_values(i, j);
      }
    }
  }

  // Encontrando el prisma conteniendo el sistema.
  double3 x0 = make_double3(0, 0, 0);
  double3 x1 = make_double3(0, 0, 0);
  for (uint atom = 0; atom < fortran_vars.atoms; atom++) {
    double3 atom_position(fortran_vars.atom_positions(atom));
    uint atom_type = fortran_vars.atom_types(atom);
    // TODO: max_radios esta al doble de lo que deberia porque el criterio no
    // esta bien aplicado.
    double max_radius = 2 * sqrt(max_function_exponent / min_exps[atom_type]);
    _DBG(cout << "tipo: " << atom_type << " " << min_exps[atom_type]
              << " radio: " << max_radius << endl);
    double3 tuple_max_radius = make_double3(max_radius, max_radius, max_radius);
    if (atom == 0) {
      x0 = atom_position - tuple_max_radius;
      x1 = atom_position + tuple_max_radius;
    } else {
      x0.x = min(x0.x, atom_position.x - max_radius);
      x0.y = min(x0.y, atom_position.y - max_radius);
      x0.z = min(x0.z, atom_position.z - max_radius);

      x1.x = max(x1.x, atom_position.x + max_radius);
      x1.y = max(x1.y, atom_position.y + max_radius);
      x1.z = max(x1.z, atom_position.z + max_radius);
    }
  }

  // Computes covalent radii relation for Becke partition.
  HostMatrix<double> becke_a;
  if (fortran_vars.becke) {
    becke_a.resize(fortran_vars.atoms, fortran_vars.atoms);
    compute_becke_a(becke_a);
  }

  // El prisma tiene vertices (x,y), con x0 el vertice inferior, izquierdo y mas
  // lejano
  // y x1 el vertice superior, derecho y mas cercano.

  // Generamos la particion en cubos.
  uint3 prism_size = ceil_uint3((x1 - x0) / little_cube_size);

  typedef vector<Point> Group;
  vector<vector<vector<Group> > > prism(
      prism_size.x,
      vector<vector<Group> >(prism_size.y, vector<Group>(prism_size.z)));

  // Inicializamos las esferas.
  vector<Group> sphere_points;
  vector<double> sphere_radius_array;
  if (sphere_radius > 0) {
    sphere_radius_array.resize(fortran_vars.atoms);
    sphere_points.resize(fortran_vars.atoms);
    for (uint atom = 0; atom < fortran_vars.atoms; atom++) {
      uint atom_shells = fortran_vars.shells(atom);
      uint included_shells = (uint)ceil(sphere_radius * atom_shells);
      double radius;
      if (included_shells == 0) {
        radius = 0;
      } else {
        double x = cos((M_PI / (atom_shells + 1)) *
                       (atom_shells - included_shells + 1));
        double rm = fortran_vars.rm(atom);
        radius = rm * (1.0 + x) / (1.0 - x);
      }
      _DBG(cout << "esfera incluye " << included_shells << " capas de "
                << atom_shells << " (radio: " << radius << ")" << endl);
      sphere_radius_array[atom] = radius;
    }
  }

  // Precomputamos las distancias entre atomos.
  for (uint i = 0; i < fortran_vars.atoms; i++) {
    const double3& atom_i_position(fortran_vars.atom_positions(i));
    double nearest_neighbor_dist = numeric_limits<double>::max();

    for (uint j = 0; j < fortran_vars.atoms; j++) {
      const double3& atom_j_position(fortran_vars.atom_positions(j));
      double dist = length(atom_i_position - atom_j_position);
      fortran_vars.atom_atom_dists(i, j) = dist;
      if (i != j) nearest_neighbor_dist = min(nearest_neighbor_dist, dist);
    }
    fortran_vars.nearest_neighbor_dists(i) = nearest_neighbor_dist;
  }

  // Computamos los puntos y los asignamos a los cubos y esferas.
  uint puntos_totales = 0;
  uint puntos_finales = 0;
  uint funciones_finales = 0;
  uint costo = 0;

  // Limpiamos las colecciones de las esferas y cubos que tengamos guardadas.
  this->clear();

  // Computa las posiciones de los puntos (y los guarda).
  for (uint atom = 0; atom < fortran_vars.atoms; atom++) {
    uint atom_shells = fortran_vars.shells(atom);
    const double3& atom_position(fortran_vars.atom_positions(atom));

    double t0 = M_PI / (atom_shells + 1);
    double rm = fortran_vars.rm(atom);

    puntos_totales += (uint)fortran_vars.grid_size * atom_shells;
    for (uint shell = 0; shell < atom_shells; shell++) {
      double t1 = t0 * (shell + 1);
      double x = cos(t1);
      double w = t0 * abs(sin(t1));
      double r1 = rm * (1.0 + x) / (1.0 - x);
      double wrad = w * (r1 * r1) * rm * 2.0 / ((1.0 - x) * (1.0 - x));

      for (uint point = 0; point < (uint)fortran_vars.grid_size; point++) {
        double3 rel_point_position =
            make_double3(fortran_vars.e(point, 0), fortran_vars.e(point, 1),
                         fortran_vars.e(point, 2));
        double3 point_position = atom_position + rel_point_position * r1;
        bool inside_prism =
            ((x0.x <= point_position.x && point_position.x <= x1.x) &&
             (x0.y <= point_position.y && point_position.y <= x1.y) &&
             (x0.z <= point_position.z && point_position.z <= x1.z));
        if (inside_prism) {
          double point_weight =
              wrad * fortran_vars.wang(point);  // integration weight

          // Atom weights for Becke partitioning.
          HostMatrix<double> atom_weights;
          if (fortran_vars.becke) {
            atom_weights.resize(fortran_vars.atoms);
            compute_becke_w(atom_weights, point_position, becke_a);
          }
          
          Point point_object(atom, shell, point, point_position,
                             point_weight, atom_weights);
          uint included_shells = (uint)ceil(sphere_radius * atom_shells);

          // Si esta capa esta muy lejos del nucleo, la modelamos como esfera,
          // sino como cubo.
          if (shell >= (atom_shells - included_shells)) {
            // Asignamos este punto a la esfera de este atomo.
            sphere_points[atom].push_back(point_object);
          } else {
            // Insertamos este punto en el cubo correspondiente.
            uint3 cube_coord =
                floor_uint3((point_position - x0) / little_cube_size);
            if (cube_coord.x >= prism_size.x || cube_coord.y >= prism_size.y ||
                cube_coord.z >= prism_size.z)
              throw std::runtime_error("Se accedio a un cubo invalido");
            prism[cube_coord.x][cube_coord.y][cube_coord.z].push_back(
                point_object);
          }
        }
      }
    }
  }

  G2G::MINCOST = getintenv("LIO_MINCOST_OFFSET", 250000);
  G2G::THRESHOLD = getintenv("LIO_SPLIT_THRESHOLD", 80);
  G2G::SPLITPOINTS = getintenv("LIO_SPLIT_POINTS", 200);

  // La grilla computada ahora tiene |puntos_totales| puntos, y |fortran_vars.m|
  // funciones.
  uint nco_m = 0;
  uint m_m = 0;

  puntos_finales = 0;
  // Completamos los parametros de los cubos y los agregamos a la particion.
  for (uint i = 0; i < prism_size.x; i++) {
    for (uint j = 0; j < prism_size.y; j++) {
      for (uint k = 0; k < prism_size.z; k++) {
        Group& points_ijk = prism[i][j][k];

        double3 cube_coord_abs = x0 + make_uint3(i, j, k) * little_cube_size;
        PointGroup<base_scalar_type>* cube_funcs;
#if GPU_KERNELS
        if (is_big_group(points_ijk))
          cube_funcs = new PointGroupGPU<base_scalar_type>();
        else
          cube_funcs = new PointGroupCPU<base_scalar_type>();
#else
        cube_funcs = new PointGroupCPU<base_scalar_type>();
#endif

        for (uint point = 0; point < prism[i][j][k].size(); point++)
          cube_funcs->add_point(prism[i][j][k][point]);

        cube_funcs->assign_functions_as_cube(cube_coord_abs, min_exps_func,
                                             min_coeff_func);

        if ((cube_funcs->total_functions_simple() == 0) ||
            (cube_funcs->number_of_points < min_points_per_cube)) {
          // Este cubo no tiene funciones o no tiene suficientes puntos.
          delete cube_funcs;
          continue;
        }
        PointGroup<base_scalar_type>* cube;

#if GPU_KERNELS
        if (is_big_group(points_ijk))
          cube = new PointGroupGPU<base_scalar_type>(
              *(static_cast<PointGroupGPU<base_scalar_type>*>(cube_funcs)));
        else
          cube = new PointGroupCPU<base_scalar_type>(
              *(static_cast<PointGroupCPU<base_scalar_type>*>(cube_funcs)));
#else
        cube = new PointGroupCPU<base_scalar_type>(
            *(static_cast<PointGroupCPU<base_scalar_type>*>(cube_funcs)));
#endif

        delete cube_funcs;

        tweights.start();
        cube->compute_weights();
        tweights.pause();

        if (cube->number_of_points < min_points_per_cube) {
          cout << "CUBE: not enough points" << endl;
          delete cube;
          continue;
        }
        cubes.push_back(cube);

        puntos_finales += cube->number_of_points;
        funciones_finales += cube->number_of_points * cube->total_functions();
        costo += cube->number_of_points *
                 (cube->total_functions() * cube->total_functions());
        nco_m += cube->total_functions() * fortran_vars.nco;
        m_m += cube->total_functions() * cube->total_functions();
      }
    }
  }

  // Si esta habilitada la particion en esferas, entonces clasificamos y las
  // agregamos a la particion tambien.
  if (sphere_radius > 0) {
    for (uint i = 0; i < fortran_vars.atoms; i++) {
      Group& sphere_i = sphere_points[i];
      assert(sphere_i.size() != 0);

      PointGroup<base_scalar_type>* sphere_funcs;
#if GPU_KERNELS
      if (is_big_group(sphere_i))
        sphere_funcs = new PointGroupGPU<base_scalar_type>();
      else
        sphere_funcs = new PointGroupCPU<base_scalar_type>();
#else
      sphere_funcs = new PointGroupCPU<base_scalar_type>();
#endif
      for (uint point = 0; point < sphere_i.size(); point++)
        sphere_funcs->add_point(sphere_i[point]);

      sphere_funcs->assign_functions_as_sphere(i, sphere_radius_array[i],
                                               min_exps_func, min_coeff_func);
      assert(sphere_funcs->total_functions_simple() != 0);
      if (sphere_funcs->number_of_points < min_points_per_cube) {
        cout << "not enough points" << endl;
        delete sphere_funcs;
        continue;
      }

      PointGroup<base_scalar_type>* sphere;
#if GPU_KERNELS
      if (is_big_group(sphere_i)) {
        sphere = new PointGroupGPU<base_scalar_type>(
            *(static_cast<PointGroupGPU<base_scalar_type>*>(sphere_funcs)));
      } else {
        sphere = new PointGroupCPU<base_scalar_type>(
            *(static_cast<PointGroupCPU<base_scalar_type>*>(sphere_funcs)));
      }
#else
      sphere = new PointGroupCPU<base_scalar_type>(
          *(static_cast<PointGroupCPU<base_scalar_type>*>(sphere_funcs)));
#endif
      delete sphere_funcs;

      assert(sphere->number_of_points != 0);
      tweights.start();
      sphere->compute_weights();
      tweights.pause();
      if (sphere->number_of_points < min_points_per_cube) {
        cout << "not enough points" << endl;
        delete sphere;
        continue;
      }
      assert(sphere->number_of_points != 0);
      spheres.push_back(sphere);

      puntos_finales += sphere->number_of_points;
      funciones_finales += sphere->number_of_points * sphere->total_functions();
      costo += sphere->number_of_points *
               (sphere->total_functions() * sphere->total_functions());
      nco_m += sphere->total_functions() * fortran_vars.nco;
      m_m += sphere->total_functions() * sphere->total_functions();
    }
  }

  // TODO fix these sorts now that spheres and cubes are pointers
  sort(spheres.begin(), spheres.end(), Sorter());
  sort(cubes.begin(), cubes.end(), Sorter());

  // Initialize the global memory pool for CUDA, with the default safety factor
  // If it is CPU, then this doesn't matter
  GlobalMemoryPool::init(G2G::free_global_memory);

  if (timer_single) cout << "  Weights: " << tweights << endl;

  for (uint i = 0; i < cubes.size(); i++) {
    cubes[i]->compute_indexes();
  }
  for (uint i = 0; i < spheres.size(); i++) {
    spheres[i]->compute_indexes();
  }

  compute_work_partition();

  timeforgroup.resize(cubes.size() + spheres.size());
  next.resize(G2G::cpu_threads + G2G::gpu_threads);

  fort_forces_ms.resize(G2G::cpu_threads + G2G::gpu_threads);
  if (!fortran_vars.OPEN) {
    rmm_outputs.resize(G2G::cpu_threads + G2G::gpu_threads);
  } else {
    rmm_outputs_a.resize(G2G::cpu_threads + G2G::gpu_threads);
    rmm_outputs_b.resize(G2G::cpu_threads + G2G::gpu_threads);
  };

  for (int i = 0; i < G2G::cpu_threads + G2G::gpu_threads; i++) {
    fort_forces_ms[i].resize(fortran_vars.max_atoms, 3);
    if (fortran_vars.OPEN) {
      rmm_outputs_a[i].resize(fortran_vars.rmm_output_a.width,
                              fortran_vars.rmm_output_a.height);
      rmm_outputs_b[i].resize(fortran_vars.rmm_output_b.width,
                              fortran_vars.rmm_output_b.height);
    } else {
      rmm_outputs[i].resize(fortran_vars.rmm_output.width,
                            fortran_vars.rmm_output.height);
    }
  }

  if (fortran_vars.becke) {
    becke_dens.resize(G2G::cpu_threads + G2G::gpu_threads);

    for (int i = 0; i < G2G::cpu_threads + G2G::gpu_threads; i++) {
      becke_dens[i].resize(fortran_vars.atoms);
    }

    if (fortran_vars.OPEN) {
      becke_spin.resize(G2G::cpu_threads + G2G::gpu_threads);

      for (int i = 0; i < G2G::cpu_threads + G2G::gpu_threads; i++) {
        becke_spin[i].resize(fortran_vars.atoms);
      }
    }
  }

  int current_gpu = 0;
  for (int i = work.size(); i < G2G::cpu_threads + G2G::gpu_threads; i++)
    work.push_back(vector<int>());

  for (uint i = 0; i < cubes.size(); i++)
    if (cubes[i]->is_big_group()) {
      work[G2G::cpu_threads + current_gpu].push_back(i);
      current_gpu = (current_gpu + 1) % G2G::gpu_threads;
    }

  for (uint i = 0; i < spheres.size(); i++)
    if (spheres[i]->is_big_group()) {
      work[G2G::cpu_threads + current_gpu].push_back(i + cubes.size());
      current_gpu = (current_gpu + 1) % G2G::gpu_threads;
    }
  if (verbose > 4) diagnostic();
}
