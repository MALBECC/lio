#include "common.h"
using namespace G2G;

bool isPowerOfTwo(unsigned long n) {
	return (n != 0) && ((n & (n - 1)) == 0);
}
