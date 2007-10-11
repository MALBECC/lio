/*-*- mode: c -*-*/

__global__ void dists(float* coords, float* dists, int atoms, int pitch);

void run_test(float* coords_cpu, float* dists_cpu, int atoms) {
	float* coords;
	float* dists;
	
	cudaMalloc(&coords, atoms * sizeof(float));
	cudaMemcpy(coords, coords_cpu, atoms * sizeof(float), cudaMemcpyHostToDevice);
	
	int pitch;	
	cudaMallocPitch(&dists, &pitch, atoms * atoms * sizeof(float), atoms);
	
	dim3 dimBlock(16, 16);
	dim3 dimSize(atoms / dimBlock.x, atoms / dimBlock.y);
	dists<<<dimSize, dimBlock>>>(coords, dists, atoms, pitch);
	
	cudaMemcpy(dists_cpu, dists, atoms * atoms * sizeof(float), cudaMemcpyDeviceToHost);	
		
	cudaFree(coords);
	cudaFree(dists);	
}

__global__ void dists(float* coords, float* dists, int atoms, int pitch)
{
	for (unsigned int i = 0; i < atoms; i++) {
		for (unsigned int j = 0; j < atoms; j++) {
			dists[i][j] = dist(coords[i], coords[j]);
		}
	}
}
