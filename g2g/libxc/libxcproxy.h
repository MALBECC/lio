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

    void closeProxy ();

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
    //printf("LibxcProxy::LibxcProxy (%d,%d,%d) \n", exchangeFunctionId, correlationFuncionalId, nspin);

    funcIdForExchange = exchangeFunctionId;
    funcIdForCorrelation = correlationFuncionalId;
    nspin = nspin;

    //printf("LibxcProxy::initFunctional () \n");
    if (xc_func_init (&funcForExchange, funcIdForExchange, nspin) != 0) {
        //fprintf (stderr, "Functional '%d' not found\n", funcIdForExchange);
    } else {
        //fprintf (stdout, "Functional '%d' inited\n", funcIdForExchange);
    }

    //printf("LibxcProxy::initFunctional () \n");
    if (xc_func_init (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
        //fprintf (stderr, "Functional '%d' not found\n", funcIdForCorrelation);
    } else {
        //fprintf (stdout, "Functional '%d' inited\n", funcIdForCorrelation);
    }

}

template <class scalar_type, int width>
LibxcProxy <scalar_type, width>::~LibxcProxy ()
{
    xc_func_end (&funcForExchange);
    xc_func_end (&funcForCorrelation);
}

template <class scalar_type, int width>
void LibxcProxy <scalar_type, width>::closeProxy ()
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
    //printf("LibxcProxy::doGGA () \n");

    double rho[1] = {dens};
    // Libxc needs the 'contracted gradient'
    double sigma[1] = {(grad.x * grad.x) + (grad.y * grad.y) + (grad.z * grad.z)};
    double exchange[1];
    double correlation[1];

    // The outputs for exchange
    double vrho [1];
    double vsigma [1];
    double v2rho [1];
    double v2rhosigma[1];
    double v2sigma [1];

    // The outputs for correlation
    double vrhoC [1];
    double vsigmaC [1];
    double v2rhoC [1];
    double v2rhosigmaC [1];
    double v2sigmaC [1];

    try {
        xc_gga (&funcForExchange, 1,
                rho,
                sigma,
                exchange,
                vrho,
                vsigma,
                v2rho,
                v2rhosigma,
                v2sigma,
                NULL, NULL, NULL, NULL);
    } catch (int exception) {
        //fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
        return;
    }

    try {
        // Now the correlation value.
        xc_gga (&funcForCorrelation, 1,
                rho,
                sigma,
                correlation,
                vrhoC,
                vsigmaC,
                v2rhoC,
                v2rhosigmaC,
                v2sigmaC,
                NULL, NULL, NULL, NULL);
    } catch (int exception) {
        //fprintf (stderr, "Exception ocurred calling xc_gga for Correlation '%d' \n", exception);
        return;
    }

    // TODO: tener el cuenta el tema del FULL_DOUBLE
    // a la vuelta del calculo en libxc. Si FULL_DOUBLE=1 -> son doubles
    // sino son floats y hay que castear antes de hacer las cuentas
    // sino da todo cero.
    ex = exchange[0];
    ec = correlation[0];
    // Merge the results for the derivatives.
    vrho[0] += vrhoC[0];
    vsigma[0] += vsigmaC[0];
    v2rho[0] += v2rhoC[0];
    v2rhosigma[0] += v2rhosigmaC[0];
    v2sigma[0] += v2sigmaC[0];

    // Now, compute y2a value.
    y2a = vrho[0] - (2 * sigma[0] * v2rhosigma[0]
            + 2 * (hess1.x + hess1.y + hess1.z) * vsigma[0]
            + 4 * v2sigma[0] * (grad.x * grad.x * hess1.x + grad.y * grad.y * hess1.y + grad.z * grad.z * hess1.z + 2 * grad.x * grad.y * hess2.x + 2 * grad.x * grad.z * hess2.y + 2 * grad.y * grad.z * hess2.z));

    return;
}

template <class scalar_type, int width>
void LibxcProxy <scalar_type, width>::doLDA(scalar_type dens, const G2G::vec_type<scalar_type, width> &grad, const G2G::vec_type<scalar_type, width> &hess1, const G2G::vec_type<scalar_type, width> &hess2, scalar_type &ex, scalar_type &ec, scalar_type &y2a)
{
    //TODO: not implemented yet!
    return;
}



#endif // LIBXCPROXY_H
