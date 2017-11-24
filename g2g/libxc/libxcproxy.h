#ifndef LIBXCPROXY_H
#define LIBXCPROXY_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <xc.h>
#include <vector>
#include <scalar_vector_types.h>


template <class scalar_type, int width>
class LibxcProxy
{
private:

    // La clase que contiene los funcionales.
    xc_func_type funcForExchange;
    xc_func_type funcForCorrelation;

    // El id del funcional.
    int funcIdForExchange;
    int funcIdForCorrelation;
    int nspin;

public:
    LibxcProxy ();
    LibxcProxy (int exchangeFunctionId, int correlationFuncionalId, int nspin);
    ~LibxcProxy ();

    void doGGA (scalar_type dens,
                const G2G::vec_type<scalar_type,width>& grad,
                const G2G::vec_type<scalar_type,width>& hess1,
                const G2G::vec_type<scalar_type,width>& hess2,
                scalar_type& ex,
                scalar_type& ec,
                scalar_type& y2a);

    void doLDA (scalar_type dens,
                const G2G::vec_type<scalar_type,width>& grad,
                const G2G::vec_type<scalar_type,width>& hess1,
                const G2G::vec_type<scalar_type,width>& hess2,
                scalar_type& ex,
                scalar_type& ec,
                scalar_type& y2a);

};

template <class scalar_type, int width>
LibxcProxy <scalar_type, width>::LibxcProxy()
{
    funcIdForExchange = 0;
    funcIdForCorrelation = 0;
    nspin = 0;
}

template <class scalar_type, int width>
LibxcProxy <scalar_type, width>::LibxcProxy (int exchangeFunctionId, int correlationFuncionalId, int nspin)
{
    printf("LibxcProxy::LibxcProxy (%d,%d,%d) \n", exchangeFunctionId, correlationFuncionalId, nspin);

    funcIdForExchange = exchangeFunctionId;
    funcIdForCorrelation = correlationFuncionalId;
    nspin = nspin;

    printf("LibxcProxy::initFunctional () \n");
    if (xc_func_init (&funcForExchange, funcIdForExchange, nspin) != 0){
        fprintf (stderr, "Functional '%d' not found\n", funcIdForExchange);
    } else {
        fprintf (stdout, "Functional '%d' inited\n", funcIdForExchange);
    }

    printf("LibxcProxy::initFunctional () \n");
    if (xc_func_init (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
        fprintf (stderr, "Functional '%d' not found\n", funcIdForCorrelation);
    } else {
        fprintf (stdout, "Functional '%d' inited\n", funcIdForCorrelation);
    }

}

template <class scalar_type, int width>
LibxcProxy <scalar_type, width>::~LibxcProxy ()
{
    xc_func_end (&funcForExchange);
    xc_func_end (&funcForCorrelation);
}

template <class scalar_type, int width>
void LibxcProxy <scalar_type, width>::doGGA(scalar_type dens,
    const G2G::vec_type<scalar_type, width> &grad,
    const G2G::vec_type<scalar_type, width> &hess1,
    const G2G::vec_type<scalar_type, width> &hess2,
    scalar_type &ex, scalar_type &ec, scalar_type &y2a)
{
    printf("LibxcProxy::doGGA () \n");

    scalar_type rho[1] = {dens};
    // Libxc necesita el valor 'contraido' del gradiente.
    scalar_type sigma[1] = {(grad.x * grad.x) + (grad.y * grad.y) + (grad.z * grad.z)};
    scalar_type exchange[1];
    scalar_type correlation[1];

    //fprintf (stdout, "rho '%lf' \n", rho[0]);
    //fprintf (stdout, "sigma'%lf' \n", sigma [0]);

    try {
        xc_gga (&funcForExchange, 1,
                rho,
                sigma,
                exchange,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL, NULL, NULL, NULL);
    } catch (int exception) {
        fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
    }

    try {
        // Ya en este punto solo necesitamos el coeficiente
        // de correlacio, el resto de las variables las
        // calculó en la llamada de más arriab.
        xc_gga (&funcForCorrelation, 1,
                rho,
                sigma,
                correlation,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL, NULL, NULL, NULL);
    } catch (int exception) {
        fprintf (stderr, "Exception ocurred calling xc_gga for Correlation '%d' \n", exception);
    }

    ex = exchange [0];
    ec = correlation [0];
    y2a = 0;
}

template <class scalar_type, int width>
void LibxcProxy <scalar_type, width>::doLDA(scalar_type dens, const G2G::vec_type<scalar_type, width> &grad, const G2G::vec_type<scalar_type, width> &hess1, const G2G::vec_type<scalar_type, width> &hess2, scalar_type &ex, scalar_type &ec, scalar_type &y2a)
{
    //TODO: aun no implementada
    return;
}



#endif // LIBXCPROXY_H
