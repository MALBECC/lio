void run_test(float* coords, float* dists, int atoms);

int main(void) {
	int atoms = 16;
	float* coords = new float[atoms];
	float* dists = new float[atoms * atoms];
	
	for (unsigned int i = 0; i < atoms; i++) coords[i] = rand();
	
	run_test(coords, dists, atoms);
	
	for (unsigned int i = 0; i < atoms; i++) {
		for (unsigned int j = 0; j < atoms; j++) {
			cout << dists[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	
	delete[] coords;
	delete[] dists;
}
