#ifndef _WEIGHT_H
#define _WEIGHT_H

double compute_point_weight(const double3& point_position, double wrad, uint atom, uint point);

void assign_cube_weights(LittleCube& cube);

#define WEIGHT_GPU 1

/* pesos en GPU */
#define GPU_WEIGHT_CUTOFF 0 // TODO: esta funcando?

/* pesos en CPU */
#define BECKE_CUTOFF 0    // cota minima para el peso
#define BECKE 1           // formula de Becke o de Stratman
#define WEIGHT_CUTOFFS 1  // acota vecinos por nucleo

#endif
