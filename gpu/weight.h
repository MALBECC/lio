#ifndef _WEIGHT_H
#define _WEIGHT_H

double compute_point_weight(const double3& point_position, double wrad, uint atom, uint point);

void assign_cube_weights(LittleCube& cube);

#define BECKE_CUTOFF 0
#define BECKE 1
#define WEIGHT_CUTOFFS 1

#endif
