#define BLOCK_DIM 16

template <class input_type>
__global__ void transpose(input_type *odata, input_type *idata, int width,
                          int height) {
  __shared__ input_type block[BLOCK_DIM][BLOCK_DIM + 1];

  // read the matrix tile into shared memory
  unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
  unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
  if ((xIndex < width) && (yIndex < height)) {
    unsigned int index_in = yIndex * width + xIndex;
    block[threadIdx.y][threadIdx.x] = idata[index_in];
  }

  __syncthreads();

  // write the transposed matrix tile to global memory
  xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
  yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
  if ((xIndex < height) && (yIndex < width)) {
    unsigned int index_out = yIndex * height + xIndex;
    odata[index_out] = block[threadIdx.x][threadIdx.y];
  }
}
