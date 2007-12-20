#include <cmath>
#include <iostream>
#include "timer.h"
#include "cuda/dist.h"
using namespace std;
using namespace G2G;

int main(void) {
	unsigned int atoms = 4096;
	
	float* coords = new float[atoms];
	float* dists = new float[atoms * atoms];
	
	for (unsigned int i = 0; i < atoms; i++) coords[i] = 1.0 * i;
	
	for (unsigned int i = 0; i < atoms; i++) {
		for (unsigned int j = 0; j < atoms; j++) {
			dists[i * atoms + j] = -57.0;
		}
	}	

	Timer total, actual;
	for (unsigned int it = 0; it < 100; it++) {
		actual.start();
		calc_dists(coords, dists, atoms);
		actual.stop();
//		cout << actual << endl;
	}
	
	/*for (unsigned int i = 0; i < atoms; i++) {
		for (unsigned int j = 0; j < atoms; j++) {
			cout << dists[i * atoms + j] << " ";
		}
		cout << endl;
	}
	cout << endl;*/
	
	delete[] coords;
	delete[] dists;
}
