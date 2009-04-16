#ifndef __G2G_COMMON_H__
#define __G2G_COMMON_H__

#ifdef _DEBUG
#define _DBG(x) x
#else
#define _DBG(x)
#endif

#define MAX_CONTRACTIONS 10

// TODO: ver cuanto se pueden subir
#define RMM_BLOCK_SIZE_XY 16

#define DENSITY_BLOCK_SIZE 256
#define FUNCTIONS_BLOCK_SIZE 128
#define FORCE_BLOCK_SIZE 256
#define WEIGHT_BLOCK_SIZE 256

// used for "types" constant memory
#define MAX_ATOMS 80

#endif
