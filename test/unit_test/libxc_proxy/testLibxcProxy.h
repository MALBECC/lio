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
void testProxy01 ()
{
    printf ("=================================== \n");
    printf ("testProxy01 \n");
    printf ("Test del constructor del proxy hacia libxc \n");
    printf ("=================================== \n");

    int nspin = XC_UNPOLARIZED;
    int functionalExchange = XC_GGA_X_PBE;
    int functionalCorrelation = XC_GGA_C_PBE;

    LibxcProxy<double,0> aProxy();
    LibxcProxy<double,0> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    printf ("=================================== \n");
}

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
    int functionalExchange = XC_GGA_X_PBE;
    int functionalCorrelation = XC_GGA_C_PBE;

    LibxcProxy<scalar_type,3> aProxy2(functionalExchange,
                                 functionalCorrelation,
                                 nspin);

    // Aca simulamos una llamada desde iteration.cpp
    scalar_type exc = 0, corr = 0, y2a = 0;
    const G2G::vec_type <scalar_type, 3> dxyz(0,0,0.447213595);
    const G2G::vec_type <scalar_type, 3> dd1(0,0,0);
    const G2G::vec_type <scalar_type, 3> dd2(0,0,0);
    scalar_type densidad = 0.1;

    try {
        aProxy2.doGGA(densidad, dxyz, dd1, dd2, exc, corr, y2a);
    } catch (int exception) {
        fprintf (stderr, "Exception ocurred calling doGGA '%d' \n", exception);
    }

//    aProxy2.closeProxy();

    // Ahora veamos los resultados.
    fprintf (stdout, "Termino el test \n");
    printf ("=================================== \n");
}

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
    int functionalExchange = XC_GGA_X_PBE;
    int functionalCorrelation = XC_GGA_C_PBE;

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
    // TODO: aca van las llamadas a lio normalmente.
    fprintf (stdout, "Usa LIO \n");
    G2G::calc_ggaCS_in<double, 3>(densidad, dxyz, dd1, dd2, exc, corr, y2a, iexch);
#endif


    fprintf (stdout, "Exchange '%lf' \n", exc);
    fprintf (stdout, "Correlation '%lf' \n", corr);
    fprintf (stdout, "Y2A '%lf' \n", y2a);

    // Ahora veamos los resultados.
//    aProxy2.closeProxy();

    fprintf (stdout, "Termino el test \n");
    printf ("=================================== \n");
}

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
    int functionalExchange = XC_GGA_X_PBE;
    int functionalCorrelation = XC_GGA_C_PBE;

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
    // TODO: aca van las llamadas a lio normalmente.
    fprintf (stdout, "Usa LIO \n");
    G2G::calc_ggaCS_in<float, 3>(densidad, dxyz, dd1, dd2, exc, corr, y2a, iexch);
#endif


    fprintf (stdout, "Exchange '%lf' \n", exc);
    fprintf (stdout, "Correlation '%lf' \n", corr);
    fprintf (stdout, "Y2A '%lf' \n", y2a);

//    aProxy2.closeProxy();

    // Ahora veamos los resultados.
    fprintf (stdout, "Termino el test \n");
    printf ("=================================== \n");
}


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
    int functionalExchange = XC_GGA_X_PBE;
    int functionalCorrelation = XC_GGA_C_PBE;

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

    aProxy2.doGGA (densidad, number_of_points, dxyz, dd1, dd2, exc, corr, y2a);
    for (int i=0; i<number_of_points; i++) {
        printf ("Point '%i' \n", i);
        printf ("Density '%lf' \n", densidad[i]);
        printf ("Exchange '%lf' \n", exc[i]);
        printf ("Correlation '%lf' \n", corr[i]);
        printf ("Y2A '%lf' \n", y2a[i]);
        printf  ("-------------------- \n");
    }

    // Make the call to the proxy.
    //int calls = 10;
    //for (int i=0; i<calls; i++) {
    //    aProxy2.doGGA(densidad, dxyz, dd1, dd2, exc, corr, y2a);
    //    printf ("Call '%i' \n", i);
    //    printf ("Density '%lf' \n", densidad);
    //    printf ("Exchange '%lf' \n", exc);
    //    printf ("Correlation '%lf' \n", corr);
    //    printf ("Y2A '%lf' \n", y2a);
    //    printf  ("-------------------- \n");
    //    densidad += 0.01;
    //}

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

void testProxy07 ()
{
    printf ("=================================== \n");
    printf ("testProxy07 \n");
    printf ("Test del init del proxy hacia libxc \n");
    printf ("=================================== \n");

    int nspin = XC_UNPOLARIZED;
    int functionalExchange = XC_GGA_X_PBE;
    int functionalCorrelation = XC_GGA_C_PBE;

    LibxcProxy<double,0> aProxy;
    aProxy.init (functionalExchange, functionalCorrelation, nspin);

    aProxy.close ();

    printf ("=================================== \n");
}


}

#endif // TESTLIBXCPROXY_H
