#ifndef __G2G_COMMON_H__
#define __G2G_COMMON_H__

#ifdef _DEBUG
#define _DBG(x) x
#else
#define _DBG(x)
#endif

#define MAX_CONTRACTIONS 7

/* Los block sizes deben ser multiplos de 16 */

#define FUNCTIONS_BLOCK_SIZE 128
#define WEIGHT_BLOCK_SIZE 128
#define RMM_BLOCK_SIZE_XY 16

#define NCO_BATCH_SIZE 8
#define DENSITY_BLOCK_SIZE 128

#define DENSITY_DERIV_BLOCK_SIZE 256
#define FORCE_BLOCK_SIZE 256

// used for "types" constant memory
#define MAX_ATOMS 100

#define REMOVE_ZEROS 1

#endif
