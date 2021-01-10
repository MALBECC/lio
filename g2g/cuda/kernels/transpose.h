#define BLOCK_DIM 16

template <class input_type>
__global__ void transpose(input_type *odata, input_type *idata, int width,
                          int height) {
  __shared__ input_type block[BLOCK_DIM][BLOCK_DIM + 1];

  // Efficiency trick: diagonal reordering.
  unsigned int blockIdx_x, blockIdx_y;
  if (width == height) {
    blockIdx_y = blockIdx.x;
    blockIdx_x = (blockIdx.x + blockIdx.y) % gridDim.x;
  } else {
    int bid = blockIdx.x + blockIdx.y * gridDim.x;
    blockIdx_y = bid % gridDim.y;
    blockIdx_x = ((bid / gridDim.y) + blockIdx_y) % gridDim.x;
  }

  // read the matrix tile into shared memory
  unsigned int xIndex = blockIdx_x * BLOCK_DIM + threadIdx.x;
  unsigned int yIndex = blockIdx_y * BLOCK_DIM + threadIdx.y;
  unsigned int index_in = yIndex * width + xIndex;
  if ((xIndex < width) && (yIndex < height)) {
    block[threadIdx.y][threadIdx.x] = idata[index_in];
  }

  __syncthreads();

  // write the transposed matrix tile to global memory
  xIndex = blockIdx_y * BLOCK_DIM + threadIdx.x;
  yIndex = blockIdx_x * BLOCK_DIM + threadIdx.y;
  unsigned int index_out = yIndex * height + xIndex;
  if ((xIndex < height) && (yIndex < width)) {
    odata[index_out] = block[threadIdx.x][threadIdx.y];
  }
}
