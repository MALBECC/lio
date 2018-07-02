#include <iostream>
using namespace std;

#include "testLibxcProxy.h"

using std::cout;
using std::endl;


void runAllProxyCPUTests () {
    libxcProxyTest::testProxy02<double> ();
    libxcProxyTest::testProxy03<double> (false);
    libxcProxyTest::testProxy03<double> (true);
    libxcProxyTest::testProxy05<double> ();
    libxcProxyTest::testProxy05<float> ();
    libxcProxyTest::testProxy06 ();
    libxcProxyTest::testProxy07 ();
    libxcProxyTest::testProxy08 ();
    libxcProxyTest::testProxy09 ();
    libxcProxyTest::testProxy10 ();
    libxcProxyTest::printFunctionalsInformationTest ();
}

void runAllProxyGPUTest () {
}

void runRandomGeneratorTest (int numberCount) {
    int i, n;
    time_t t;

    n = numberCount;

    /* Intializes random number generator */
    srand((unsigned) time(&t));

    /* Print 5 random numbers from 0 to 49 */
    for( i = 0 ; i < n ; i++ ) {
        printf("%lf \n", (rand() % 10) * 0.1f);
    }
}


int main()
{
    cout << "Test: Lio met Libxc - BEGIN" << endl;
    try {
        runAllProxyCPUTests();
    } catch (int e) {
	cout << "An exception occurred. Exception Nr. " << e << endl;
	exit (EXIT_FAILURE);
    }
    cout << "Test: Lio mets Libxc - END" << endl;
    return 0;
}

