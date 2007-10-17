#include <cmath>
#include <iostream>
#include "timer.h"
#include "matrix.h"
#include "cuda/accum.h"
using namespace std;
using namespace G2G;

int main(void) {
	unsigned int n = 32;
	Matrix input(n, n);
	Matrix output(-57.0f);
	
	for (unsigned int i = 0; i < n * n; i++) input.data[i] = 1.0;
	
	Timer total, actual;
	for (unsigned int it = 0; it < 100; it++) {
		actual.start();
		calc_accum(input, output);
		actual.stop();
		cout << actual << endl;
	}
	
	cout << "resultado: " << *output.data << endl;
}
