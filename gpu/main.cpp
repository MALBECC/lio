#include <cmath>
#include <iostream>
#include "timer.h"
using namespace std;
using namespace G2G;

extern "C" void run_test(float* coords, float* dists, unsigned int atoms);

int main(void) {
	unsigned int atoms = 16;
	float* coords = new float[atoms];
	float* dists = new float[atoms * atoms];
	
	for (unsigned int i = 0; i < atoms; i++) coords[i] = 1.0 * i;

	Timer total, actual;
	for (unsigned int it = 0; it < 100; it++) {
		actual.start();
		run_test(coords, dists, atoms);
		actual.stop();
		cout << actual << endl;
	}
	
	for (unsigned int i = 0; i < atoms; i++) {
		for (unsigned int j = 0; j < atoms; j++) {
			cout << dists[i * atoms + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	
	delete[] coords;
	delete[] dists;
}
