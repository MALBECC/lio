#ifndef __EXCHNUM_CONSTANTS_H__
#define __EXCHNUM_CONSTANTS_H__

#define PI 3.141592654f

#define MAX_CONTRACTIONS 10

// This should correspond to the maximum number in 'layers' and 'layers2'
#define MAX_LAYERS 50

// TODO jugar bajando estos valores
#define RMM_BLOCK_SIZE_X 16
#define RMM_BLOCK_SIZE_Y 16

#define ENERGY_FORCE_BLOCK_SIZE_X 1
#define ENERGY_FORCE_BLOCK_SIZE_Y 64

//#define ENERGY_S_BLOCK_SIZE_X 1
#define ENERGY_S_BLOCK_SIZE_X 4
#define ENERGY_M_BLOCK_SIZE_X 1
//#define ENERGY_M_BLOCK_SIZE_X 2
#define ENERGY_B_BLOCK_SIZE_X 1

//#define ENERGY_S_BLOCK_SIZE_Y 128
#define ENERGY_S_BLOCK_SIZE_Y 32
#define ENERGY_M_BLOCK_SIZE_Y 128
//#define ENERGY_M_BLOCK_SIZE_Y 64
#define ENERGY_B_BLOCK_SIZE_Y 128

#define FORCE_BLOCK_SIZE 64

// used for "types" constant memory
#define MAX_ATOMS 60

#define GPU_LAYERS_1 0
#define GPU_LAYERS_2 1

const uint cpu_layers[] = {
	30,30,35,35,35,35,35,35,35,35,
	40,40,40,40,40,40,40,40,
	45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,
	50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50
};

__device__ __constant__ uint gpu_layers_1[] = {
  30,30,35,35,35,35,35,35,35,35,
  40,40,40,40,40,40,40,40,
  45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,
  50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50
};

const uint cpu_layers2[] = {
	20,20,25,25,25,25,25,25,25,25,
	30,30,30,30,30,30,30,30,
	35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,
	40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40
};

__device__ __constant__ uint gpu_layers_2[] = {
	20,20,25,25,25,25,25,25,25,25,
	30,30,30,30,30,30,30,30,
	35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,
	40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40
};	

__device__ __constant__ float rm_factor[] = {
	0.330702203610512, 0.878722998165075, 1.37005198638641, 0.992106610831537,
	0.803133923054101, 0.661404407221024, 0.614161235276665, 0.566918063332307,
	0.472431719443589, 0.670853041609896, 3.40150837999384, 1.41729515833077,
	1.18107929860897, 1.0393497827759, 0.944863438887178, 0.944863438887178,
	0.944863438887178, 1.88972687777436, 2.07869956555179, 1.70075418999692,
	1.51178150221948, 1.32280881444205, 1.27556564249769, 1.32280881444205,
	1.32280881444205, 1.32280881444205, 1.27556564249769, 1.27556564249769,
	1.27556564249769, 1.27556564249769, 1.22832247055333, 1.18107929860897,
	1.08659295472025, 1.08659295472025, 1.08659295472025, 1.08659295472025,
	2.22042908138487, 1.88972687777436, 1.70075418999692, 1.46453833027513,
	1.37005198638641, 1.37005198638641, 1.27556564249769, 1.22832247055333,
	1.27556564249769, 1.32280881444205, 1.51178150221948, 1.46453833027513,
	1.46453833027513, 1.37005198638641, 1.37005198638641, 1.32280881444205,
	1.32280881444205, 1.2377711049422
};
		
/* grid */
__device__ __constant__ float3 gpu_point_positions[EXCHNUM_BIG_GRID_SIZE];
__device__ __constant__ float gpu_wang[EXCHNUM_BIG_GRID_SIZE];

__device__ __constant__ uint gpu_types[MAX_ATOMS];
__device__ __constant__ float3 gpu_atom_positions[MAX_ATOMS];

/* pot_kernel constants */
__device__ __constant__ float pot_alpha = -0.738558766382022447;
__device__ __constant__ float pot_gl = 0.620350490899400087;
	
__device__ __constant__ float pot_vosko_a1 = 0.03109205;
__device__ __constant__ float pot_vosko_b1 = 3.72744;
__device__ __constant__ float pot_vosko_c1 = 12.9352;
__device__ __constant__ float pot_vosko_x0 = -0.10498;
__device__ __constant__ float pot_vosko_q = 6.15199066246304849;
__device__ __constant__ float pot_vosko_a16 = 0.005182008333;
__device__ __constant__ float pot_vosko_a2 = 0.015546025;
__device__ __constant__ float pot_vosko_b2 = 7.06042;
__device__ __constant__ float pot_vosko_c2 = 18.0578;
__device__ __constant__ float pot_vosko_x02 = -0.32500;
__device__ __constant__ float pot_vosko_q2 = 4.7309269;
__device__ __constant__ float pot_vosko_a26 = 0.0025910042;

__device__ __constant__ uint Iexch_gpu = 0;

#endif
