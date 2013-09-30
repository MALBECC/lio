#ifndef __G2G_COMMON_H__
#define __G2G_COMMON_H__

#include <stdexcept>
#include <sys/types.h>

#ifdef _DEBUG
#define _DBG(x) x
#define assert(x) if (!(x)) throw std::runtime_error("assertion failed: " #x);
#else
#define _DBG(x)
#define assert(x)
#endif

#define MAX_CONTRACTIONS 7 

/* Los block sizes deben ser multiplos de 16 */

#define FUNCTIONS_BLOCK_SIZE 128
#define WEIGHT_BLOCK_SIZE 128
#define RMM_BLOCK_SIZE_XY 16

#define DENSITY_BATCH_SIZE 128
#define DENSITY_ACCUM_BLOCK_SIZE 128
#define DENSITY_BLOCK_SIZE 128

/*#define DENSITY_DERIV_BLOCK_SIZE 128
#define DENSITY_DERIV_BATCH_SIZE 32
#define DENSITY_DERIV_BATCH_SIZE2 32*/

#define DENSITY_DERIV_BLOCK_SIZE 128
#define DENSITY_DERIV_BATCH_SIZE 128
#define DENSITY_DERIV_BATCH_SIZE2 128

#define FORCE_BLOCK_SIZE 256

// used for "types" constant memory
#define MAX_ATOMS 200

#define COMPUTE_RMM 					0
#define COMPUTE_ENERGY_ONLY		1
#define COMPUTE_ENERGY_FORCE	2
#define COMPUTE_FORCE_ONLY		3

#endif
