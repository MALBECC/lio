#ifndef TESTLIBXCPROXY_H
#define TESTLIBXCPROXY_H

#include <stdio.h>
#include <xc.h>
#include "../../../g2g/fix_compile.h"
#include "../../../g2g/pointxc/calc_ggaCS.h"
#include "../../../g2g/libxc/libxcproxy.h"

namespace libxcProxyTest
{
#define PEPITO

///////////////////////////////////////////
// TestProxy01
// Simple constructor test
//
//
void testProxy01 ()
{
    printf ("=================================== \n");
    printf ("testProxy01 \n");
    printf ("Test del constructor del proxy hacia libxc \n");
    printf ("=================================== \n");

    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<double,0> aProxy();
    LibxcProxy<double,0> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    printf ("=================================== \n");
}

///////////////////////////////////////////
// TestProxy02
// This test simulates a call to Libxc
// like if we where inside the iteration.cpp class
//
template<class scalar_type>
void testProxy02 ()
{
    printf ("=================================== \n");
    printf ("testProxy02 \n");
    printf ("Test del proxy para libx con parametros que usa lio \n");
    printf ("=================================== \n");

    // Con este parametro determinamos si es capa
    // abierta o capa cerrada (OS, CS).
    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<scalar_type,3> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    // Aca simulamos una llamada desde iteration.cpp
    scalar_type exc = 0, corr = 0, y2a = 0;
    const G2G::vec_type <scalar_type, 3> dxyz(0,0,0.447213595);
    const G2G::vec_type <scalar_type, 3> dd1(0,0,0);
    const G2G::vec_type <scalar_type, 3> dd2(0,0,0);
    scalar_type densidad = 0.1;

    aProxy2.doGGA(densidad, dxyz, dd1, dd2, exc, corr, y2a);

    fprintf (stdout, "Termino el test \n");
    printf ("=================================== \n");
}

///////////////////////////////////////////
// TestProxy03
// Test the Proxy to libxc with the
// input in the Lio-like format
//
template<class scalar_type>
void testProxy03 (bool overrideForLio)
{
    printf ("=================================== \n");
    printf ("testProxy03 \n");
    printf ("Simulacro de como seria la llamada en lio \n");
    printf ("=================================== \n");

    // Con este parametro determinamos si es capa
    // abierta o capa cerrada (OS, CS).
    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<scalar_type,3> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    // Aca simulamos una llamada desde iteration.cpp
    scalar_type exc = 0, corr = 0, y2a = 0;
    const G2G::vec_type <scalar_type, 3> dxyz(0,0,0.447213595);
    const G2G::vec_type <scalar_type, 3> dd1(0,0,0);
    const G2G::vec_type <scalar_type, 3> dd2(0,0,0);
    scalar_type densidad = 0.1;
    int iexch = 9;

#if !USE_LIBXC
    if (overrideForLio) {
        fprintf (stdout, "LIBXC configurado, pero usa LIO \n");
        G2G::calc_ggaCS_in<double, 3>(densidad, dxyz, dd1, dd2, exc, corr, y2a, iexch);
    } else {
        fprintf (stdout, "Usa LIBXC \n");
        try {
            aProxy2.doGGA(densidad, dxyz, dd1, dd2, exc, corr, y2a);
        } catch (int exception) {
            fprintf (stderr, "Exception ocurred calling doGGA '%d' \n", exception);
        }
    }
#else
    fprintf (stdout, "Usa LIO \n");
    G2G::calc_ggaCS_in<double, 3>(densidad, dxyz, dd1, dd2, exc, corr, y2a, iexch);
#endif


    fprintf (stdout, "Exchange '%lf' \n", exc);
    fprintf (stdout, "Correlation '%lf' \n", corr);
    fprintf (stdout, "Y2A '%lf' \n", y2a);
    fprintf (stdout, "Termino el test \n");
    printf ("=================================== \n");
}

///////////////////////////////////////////
// TestProxy04
// This test simulates a call to Libxc
// like if we where inside the iteration.cpp class
//
template<class scalar_type>
void testProxy04 (bool overrideForLio)
{
    printf ("=================================== \n");
    printf ("testProxy03 \n");
    printf ("Simulacro de como seria la llamada en lio \n");
    printf ("=================================== \n");

    // Con este parametro determinamos si es capa
    // abierta o capa cerrada (OS, CS).
    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<float,3> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    // Aca simulamos una llamada desde iteration.cpp
    float exc = 0, corr = 0, y2a = 0;
    const G2G::vec_type <float, 3> dxyz(0,0,0.447213595);
    const G2G::vec_type <float, 3> dd1(0,0,0);
    const G2G::vec_type <float, 3> dd2(0,0,0);
    float densidad = 0.1;
    int iexch = 9;

    // This is how is gonna to be implemented in iteration.cpp
#if !USE_LIBXC
    if (overrideForLio) {
        fprintf (stdout, "LIBXC configurado, pero usa LIO \n");
        G2G::calc_ggaCS_in<float, 3>(densidad, dxyz, dd1, dd2, exc, corr, y2a, iexch);
    } else {
        fprintf (stdout, "Usa LIBXC \n");
        try {
            aProxy2.doGGA(densidad, dxyz, dd1, dd2, exc, corr, y2a);
        } catch (int exception) {
            fprintf (stderr, "Exception ocurred calling doGGA '%d' \n", exception);
        }
    }
#else
    fprintf (stdout, "Usa LIO \n");
    G2G::calc_ggaCS_in<float, 3>(densidad, dxyz, dd1, dd2, exc, corr, y2a, iexch);
#endif


    fprintf (stdout, "Exchange '%lf' \n", exc);
    fprintf (stdout, "Correlation '%lf' \n", corr);
    fprintf (stdout, "Y2A '%lf' \n", y2a);

    // Ahora veamos los resultados.
    fprintf (stdout, "Termino el test \n");
    printf ("=================================== \n");
}

////////////////////////////////////
// Test Proxy 05
// Several call to the functionals
// through the Proxy
//
template<class scalar_type>
void testProxy05 ()
{
    printf ("=================================== \n");
    printf ("testProxy05 \n");
    printf ("Varias llamadas a los funcionales   \n");
    printf ("con el mismo proxy                  \n");
    printf ("=================================== \n");

    // Con este parametro determinamos si es capa
    // abierta o capa cerrada (OS, CS).
    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<scalar_type,3> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    // Aca simulamos una llamada desde iteration.cpp
    scalar_type exc = 0, corr = 0, y2a = 0;
    const G2G::vec_type <scalar_type, 3> dxyz(0,0,0.447213595);
    const G2G::vec_type <scalar_type, 3> dd1(0,0,0);
    const G2G::vec_type <scalar_type, 3> dd2(0,0,0);
    scalar_type densidad = 0.1;

    int calls = 10;
    for (int i=0; i<calls; i++) {
        aProxy2.doGGA(densidad, dxyz, dd1, dd2, exc, corr, y2a);
        fprintf (stdout, "Call '%i' \n", i);
        fprintf (stdout, "Density '%lf' \n", densidad);
        fprintf (stdout, "Exchange '%lf' \n", exc);
        fprintf (stdout, "Correlation '%lf' \n", corr);
        fprintf (stdout, "Y2A '%lf' \n", y2a);
        printf  ("-------------------- \n");
        densidad += 0.01;
    }

    printf ("Termino el test \n");
    printf ("=================================== \n");
}

double getRandomNumber (int numberCount) {
    int i, n;
    time_t t;

    n = numberCount;

    /* Intializes random number generator */
    srand((unsigned) time(&t));

    double number = 0;
    for( i = 0 ; i < n ; i++ ) {
        //printf("%lf \n", (rand() % 10) * 0.1f);
	number = rand() * 0.1f;
    }
    return number;
}

//////////////////////////////////////////
// TestProxy06
// Simple call to libxc though the Proxy
//
void testProxy06 ()
{
    printf ("=================================== \n");
    printf ("TestProxy06 \n");
    printf ("=================================== \n");

    // Con este parametro determinamos si es capa
    // abierta o capa cerrada (OS, CS).
    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;
    int number_of_points = 10;

    // Create the proxy.
    LibxcProxy<double,3> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    // Parameters
    double* exc;
    double* corr;
    double* y2a;
    double* densidad;
    G2G::vec_type<double,3>* dxyz;
    G2G::vec_type<double,3>* dd1;
    G2G::vec_type<double,3>* dd2;

    // Now alloc memory for the data
    exc 	= (double*)malloc(sizeof(double)*number_of_points);
    corr 	= (double*)malloc(sizeof(double)*number_of_points);
    y2a 	= (double*)malloc(sizeof(double)*number_of_points);
    densidad 	= (double*)malloc(sizeof(double)*number_of_points);

    dxyz = (G2G::vec_type<double,3>*)malloc(sizeof(G2G::vec_type<double,3>)*number_of_points);
    dd1  = (G2G::vec_type<double,3>*)malloc(sizeof(G2G::vec_type<double,3>)*number_of_points);
    dd2  = (G2G::vec_type<double,3>*)malloc(sizeof(G2G::vec_type<double,3>)*number_of_points);

    // Fill the input parameters.
    for (int i=0; i<number_of_points; i++) {
	exc[i] = 0;
	corr[i] = 0;
	y2a[i] = 0;
	densidad[i] = 0.1*i;
	dxyz[i] = G2G::vec_type<double,3>(0,0,0.447213595);
	dd1[i] = G2G::vec_type<double,3>(0,0,0);
	dd2[i] = G2G::vec_type<double,3>(0,0,0);
    }

    // Make the call to the proxy.
    aProxy2.doGGA (densidad, number_of_points, dxyz, dd1, dd2, exc, corr, y2a);
    for (int i=0; i<number_of_points; i++) {
        printf ("Point '%i' \n", i);
        printf ("Density '%lf' \n", densidad[i]);
        printf ("Exchange '%lf' \n", exc[i]);
        printf ("Correlation '%lf' \n", corr[i]);
        printf ("Y2A '%lf' \n", y2a[i]);
        printf  ("-------------------- \n");
    }

    // Free the memory.
    free(exc);
    free(corr);
    free(y2a);
    free(densidad);
    free(dxyz);
    free(dd1);
    free(dd2);

    printf ("Termino el test \n");
    printf ("=================================== \n");
}

/////////////////////////////
// TestProxy07
// Constructor test
//
void testProxy07 ()
{
    printf ("=================================== \n");
    printf ("testProxy07 \n");
    printf ("Test del init del proxy hacia libxc \n");
    printf ("=================================== \n");

    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<float,3> aProxy;
    LibxcProxy<float,3> anotherProxy (functionalExchange, functionalCorrelation, nspin);

    aProxy.closeProxy ();
    anotherProxy.closeProxy ();

    printf ("=================================== \n");
}

/////////////////////////////
// TestProxy08
// Constructor test
//
void testProxy08 ()
{
    printf ("=================================== \n");
    printf ("testProxy08 \n");
    printf ("Test del init del proxy hacia libxc \n");
    printf ("=================================== \n");

    int nspin = XC_UNPOLARIZED;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<double,0> aProxy;
    aProxy.init (functionalExchange, functionalCorrelation, nspin);

    aProxy.closeProxy ();

    printf ("=================================== \n");
}

/////////////////////////////
// Test: TestProxy09
//
// Calling xc_gga_xfc test 01
// We use the cpu version LibxcProxy component to make calls
// to the Libxc xc_gga_fxc function.

void testProxy09 ()
{
    printf ("=================================== \n");
    printf ("testProxy09 \n");
    printf (" \n");
    printf ("=================================== \n");

    int nspin = 1;
    int functionalExchange = 101;
    int functionalCorrelation = 130;
    int number_of_points = 10;

    LibxcProxy<double,3> aProxy;
    aProxy.init (functionalExchange, functionalCorrelation, nspin);

    // Parameters
    double* rho;
    double* sigma;
    double* v2rho2;
    double* v2rhosigma;
    double* v2sigma2;

    // Now alloc memory for the data
    rho = (double*)malloc(sizeof(double)*number_of_points);
    sigma = (double*)malloc(sizeof(double)*number_of_points);
    v2rho2 = (double*)malloc(sizeof(double)*number_of_points);
    v2rhosigma = (double*)malloc(sizeof(double)*number_of_points);
    v2sigma2 = (double*)malloc(sizeof(double)*number_of_points);


    // Fill the input parameters.
    for (int i=0; i<number_of_points; i++) {
	rho[i] = 0.1*i;
	sigma[i] = 0.1;
	v2rho2[i] = 0;
	v2rhosigma[i] = 0;
	v2sigma2[i] = 0;
    }

    // Make the call to the proxy.
    aProxy.doGGA (rho, sigma, number_of_points, v2rho2, v2rhosigma, v2sigma2);

    for (int i=0; i<number_of_points; i++) {
        printf ("Point '%i' \n", i);
        printf ("rho '%lf' \n", rho[i]);
        printf ("sigma '%lf' \n", sigma[i]);
        printf ("v2rho '%lf' \n", v2rho2[i]);
        printf ("v2rhosigma '%lf' \n", v2rhosigma[i]);
        printf ("v2sigma2 '%lf' \n", v2sigma2[i]);
        printf  ("-------------------- \n");
    }

    // Free the memory.
    free(rho);
    free(sigma);
    free(v2rho2);
    free(v2rhosigma);
    free(v2sigma2);

    aProxy.closeProxy ();

    printf ("=================================== \n");
}


/////////////////////////////
// Test: TestProxy10
//
// Calling xc_gga_xfc test 02
// We use the cpu version LibxcProxy component to make calls
// to the Libxc xc_gga_fxc function.
//
void testProxy10 ()
{
    printf ("=================================== \n");
    printf ("testProxy10 \n");
    printf (" \n");
    printf ("=================================== \n");

    int nspin = 1;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<double,3> aProxy;
    aProxy.init (functionalExchange, functionalCorrelation, nspin);

    // Parameters
    double rho;
    double sigma;
    double v2rho2;
    double v2rhosigma;
    double v2sigma2;

    // Set the data value
    rho = 0.1;
    sigma = 0.1;
    v2rho2 = 0;
    v2rhosigma = 0;
    v2sigma2 = 0;


    // Make the call to the proxy.
    aProxy.doGGA (&rho, &sigma, 1, &v2rho2, &v2rhosigma, &v2sigma2);

    printf ("rho '%lf' \n", rho);
    printf ("sigma '%lf' \n", sigma);
    printf ("v2rho '%lf' \n", v2rho2);
    printf ("v2rhosigma '%lf' \n", v2rhosigma);
    printf ("v2sigma2 '%lf' \n", v2sigma2);
    printf  ("-------------------- \n");

    aProxy.closeProxy ();

    printf ("=================================== \n");
}

void printFunctionalsInformationTest () 
{
    int nspin = 1;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<double,3> aProxy;

    aProxy.printFunctionalsInformation (functionalExchange, functionalCorrelation);

}


}

#endif // TESTLIBXCPROXY_H
