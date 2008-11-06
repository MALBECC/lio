#ifndef __G2G_COMMON_H__
#define __G2G_COMMON_H__

#ifdef _DEBUG
#define _DBG(x) x
#else
#define _DBG(x)
#endif

#define PI 3.141592654f

#define MAX_CONTRACTIONS 10

// This should correspond to the maximum number in 'layers' and 'layers2'
#define MAX_LAYERS 50

// TODO jugar bajando estos valores
#define RMM_BLOCK_SIZE_X 16
#define RMM_BLOCK_SIZE_Y 16

// used for "types" constant memory
#define MAX_ATOMS 80

#define GPU_LAYERS_1 0
#define GPU_LAYERS_2 1

#endif
