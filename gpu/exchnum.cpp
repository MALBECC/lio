/* -*- mode: c -*- */
#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "common.h"
#include "init.h"
#include "cuda_extra.h"
#include "../matrix.h"
#include "exchnum.h"
#include "gpu_variables.h"
#include "../timer.h"
#include "double.h"
#include "cubes.h"

using namespace G2G;
using namespace std;

void cpu_compute_functions();

void cpu_compute_cube_functions(void) {
  cout << "<===== computing functions (CPU) ========>" << endl;
  HostMatrixFloat3 points_position_cpu;
  HostMatrixFloat2 factor_ac_cpu;
  HostMatrixUInt nuc_cpu;
  HostMatrixUInt contractions_cpu;

  Timer t1;
  t1.sync();
  t1.start();

  for (list<LittleCube>::iterator it = final_cube.begin(); it != final_cube.end(); ++it) {
    LittleCube& cube = *it;

    /** Load cube points **/
    points_position_cpu.resize(cube.number_of_points, 1);
    {
      uint i = 0;
      for (list<Point>::const_iterator p = cube.points.begin(); p != cube.points.end(); ++p, ++i) {
        points_position_cpu.get(i) = make_float3(p->position.x, p->position.y, p->position.z);
      }
    }

    /* Load cube functions */
    uint cube_m = cube.s_functions + cube.p_functions * 3 + cube.d_functions * 6;
    uint cube_spd = cube.s_functions + cube.p_functions + cube.d_functions;
    {
      factor_ac_cpu.resize(cube_spd, MAX_CONTRACTIONS);
      nuc_cpu.resize(cube_spd, 1);
      contractions_cpu.resize(cube_spd, 1);

      uint i = 0;
      for (set<uint>::const_iterator func = cube.functions.begin(); func != cube.functions.end(); ++func, ++i) {
        nuc_cpu.get(i) = fortran_vars.nucleii.get(*func) - 1;
        contractions_cpu.get(i) = fortran_vars.contractions.get(*func);
        assert(contractions_cpu.get(i) <= MAX_CONTRACTIONS);

        for (unsigned int k = 0; k < contractions_cpu.get(i); k++)
          factor_ac_cpu.get(i, k) = make_float2(fortran_vars.a_values.get(*func, k), fortran_vars.c_values.get(*func, k));
      }
    }

    /** Compute Functions **/
    cube.function_values.resize(cube_m, cube.number_of_points);
    if (fortran_vars.do_forces) cube.gradient_values.resize(cube_m, cube.number_of_points);

    //cout << "points: " << threads.x << " " << threadGrid.x << " " << threadBlock.x << endl;

    cpu_compute_functions(cube, fortran_vars.do_forces, contractions_cpu, factor_ac_cpu, nuc_cpu, cube_spd);

#if 0
    if (fortran_vars.grid_type == BIG_GRID) {
      cout << "s_funcs: " << cube.s_functions << " p_funcs " << cube.p_functions << " d_funcs " << cube.d_functions << endl;
      HostMatrixFloat functions_cpu(cube.function_values);
      HostMatrixFloat3 gradients_cpu(cube.gradient_values);
      uint i = 0;
      for (list<Point>::const_iterator p = cube.points.begin(); p != cube.points.end(); ++p, ++i) {
        uint func_idx = 0;
        for (set<uint>::const_iterator func = cube.functions.begin(); func != cube.functions.end(); ++func, ++func_idx) {
          if (fortran_vars.nucleii.get(*func) - 1 != 0) continue;
          if (func_idx < cube.s_functions)
            cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(func_idx, i).x << " " << gradients_cpu.get(func_idx, i).y << " " << gradients_cpu.get(func_idx, i).z << " F " << functions_cpu.get(func_idx, i) << " " << func_idx << endl;
          else if (func_idx < cube.p_functions + cube.s_functions) {
            uint p_idx = 3 * (func_idx - cube.s_functions) + cube.s_functions;
            for (uint j = 0; j < 3; j++)
              cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(p_idx + j, i).x << " " << gradients_cpu.get(p_idx + j, i).y << " " << gradients_cpu.get(p_idx + j, i).z << " F " << functions_cpu.get(p_idx + j, i) << " " << p_idx + j << endl;
          } else {
            uint s_idx = cube.s_functions + cube.p_functions * 3 + 6 * (func_idx - cube.s_functions - cube.p_functions);
            for (uint j = 0; j < 6; j++)
              cout << "* point (" << p->atom << "," << p->shell << "," << p->point << ") - Fg(" << *func << ")=" << gradients_cpu.get(s_idx + j, i).x << " " << gradients_cpu.get(s_idx + j, i).y << " " << gradients_cpu.get(s_idx + j, i).z << " F " << functions_cpu.get(s_idx + j, i) << " " << s_idx + j << endl;

          }
          //				cout << "* point " << p->position.x << " " << p->position.y << " " << p->position.z << " " << functions_cpu.get(p_idx, i) << endl;
        }
      }
    }
#endif
  }

  t1.sync();
  t1.stop();
  cout << "TIMER: funcs: " << t1 << endl;
}

static void cpu_compute_function(uint idx, float3 point_position, const HostMatrixUInt& contractions, const HostMatrixFloat2& factor_ac, const HostMatrixUInt& nuc,
  uint spd, float& t, float& tg, float3& v) {
  float3 atom_nuc_position = fortran_vars.atom_positions.get(nuc.get(idx));
  v = point_position - atom_nuc_position;
  float dist = length2(v);
  uint func_contractions = contractions.get(idx);

  t = 0.0f;
  if (do_forces) tg = 0.0f;

  for (uint contraction = 0; contraction < func_contractions; contraction++) {
    float2 curr_factor_ac = factor_ac.get(idx, contraction);
    float rexp = curr_factor_ac.x * dist;
    //if (rexp > 30.0f) continue;
    float t0 = expf(-rexp) * curr_factor_ac.y;
    t += t0;
    if (do_forces) tg += t0 * curr_factor_ac.x;
  }
}

void cpu_compute_functions(LittleCube& cube, bool do_forces, const HostMatrixUInt& contractions, const HostMatrixFloat2& factor_ac, const HostMatrixUInt& nuc,
  uint spd) {
  float t, tg;
  float3 v;

  uint point = 0;
  for (list<Point>::const_iterator p = cube.points.begin(); p != cube.points.end(); ++p, ++point) {
    // s functions
    for (uint i = 0; i < cube.s_functions; i++, big_i++) {
      cpu_compute_function<do_forces>(i, p->position, contractions, factor_ac, nuc, spd, t, tg, v);

      // TODO: invertir los indices
      cube.function_values.get(big_i, point) = t;
      if (do_forces) cube.gradient_values.get(big_i, point) = v * (2.0f * tg);
    }

    // p functions
    for (uint i = 0; i < cube.p_functions; i++, big_i += 3) {
      cpu_compute_function<do_forces>(cube.s_functions + i, p->position, contractions, factor_ac, nuc, spd, t, tg, v);

      cube.function_values.get(big_i + 0, point) = v.x * t;
      cube.function_values.get(big_i + 1, point) = v.y * t;
      cube.function_values.get(big_i + 2, point) = v.z * t;

      if (do_forces) {
        cube.gradient_values.get(big_i + 0, point) = v * 2.0f * v.x * tg - make_float3(t, 0, 0);
        cube.gradient_values.get(big_i + 1, point) = v * 2.0f * v.y * tg - make_float3(0, t, 0);
        cube.gradient_values.get(big_i + 2, point) = v * 2.0f * v.z * tg - make_float3(0, 0, t);
      }
    }

    // d functions
    for (uint i = 0; i < cube.d_functions; i++, big_i += 6) {
      cpu_compute_function<do_forces>(cube.s_functions + cube.d_functions + i, p->position, contractions, factor_ac, nuc, spd, t, tg, v);

      float tx = t * v.x;
      float ty = t * v.y;
      float tz = t * v.z;

      cube.function_values.get(big_i + 0, point) = tx * v.x * gpu_normalization_factor;
      cube.function_values.get(big_i + 1, point) = ty * v.x;
      cube.function_values.get(big_i + 2, point) = ty * v.y * gpu_normalization_factor;
      cube.function_values.get(big_i + 3, point) = tz * v.x;
      cube.function_values.get(big_i + 4, point) = tz * v.y;
      cube.function_values.get(big_i + 5, point) = tz * v.z * gpu_normalization_factor;

      if (do_forces) {
        float tgx = tg * v.x;
        float tgy = tg * v.y;
        float tgz = tg * v.z;

        cube.gradient_values.get(big_i + 0, point) = v * 2.0f * tgx * v.x * gpu_normalization_factor - make_float3(2 * tx * gpu_normalization_factor, 0, 0);
        cube.gradient_values.get(big_i + 1, point) = v * 2.0f * tgy * v.x - make_float3(ty, tx, 0);
        cube.gradient_values.get(big_i + 2, point) = v * 2.0f * tgy * v.y * gpu_normalization_factor - make_float3(0, 2 * ty * gpu_normalization_factor, 0);
        cube.gradient_values.get(big_i + 3, point) = v * 2.0f * tgz * v.x - make_float3(tz, 0, tx);
        cube.gradient_values.get(big_i + 4, point) = v * 2.0f * tgz * v.y - make_float3(0, tz, ty);
        cube.gradient_values.get(big_i + 5, point) = v * 2.0f * tgz * v.z * gpu_normalization_factor - make_float3(0, 0, 2 * tz * gpu_normalization_factor);
      }
    }
  }

}