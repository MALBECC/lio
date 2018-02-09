#ifndef LIBXCPROXY_H
#define LIBXCPROXY_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <xc.h>
#include <vector>
//#include <scalar_vector_types.h>
#include "../scalar_vector_types.h"

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

    void doGGA (scalar_type* dens,
                const int number_of_points,
		const G2G::vec_type<scalar_type,width>* grad,
                const G2G::vec_type<scalar_type,width>* hess1,
                const G2G::vec_type<scalar_type,width>* hess2,
                scalar_type* ex,
                scalar_type* ec,
                scalar_type* y2a);

    void doGGA (scalar_type* dens,
                const int number_of_points,
		const scalar_type* contracted_grad,
		const G2G::vec_type<scalar_type,width>* grad,
                const G2G::vec_type<scalar_type,width>* hess1,
                const G2G::vec_type<scalar_type,width>* hess2,
                scalar_type* ex,
                scalar_type* ec,
                scalar_type* y2a);

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
LibxcProxy <scalar_type, width>::LibxcProxy (int exchangeFunctionalId, int correlationFuncionalId, int nSpin)
{
#ifdef _DEBUG
    printf("LibxcProxy::LibxcProxy (%i,%i,%i) \n", exchangeFunctionalId, correlationFuncionalId, nSpin);
#endif


    funcIdForExchange = exchangeFunctionalId;
    funcIdForCorrelation = correlationFuncionalId;
    nspin = nSpin;

    if (xc_func_init (&funcForExchange, funcIdForExchange, nspin) != 0) {
        fprintf (stderr, "Functional '%d' not found\n", funcIdForExchange);
	exit(-1);
    }

    if (xc_func_init (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
	fprintf (stderr, "Functional '%d' not found\n", funcIdForCorrelation);
	exit(-1);
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
#ifdef _DEBUG
    printf("LibxcProxy::doGGA (...) \n");
#endif

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
__global__ void joinResults(scalar_type* ex, scalar_type* exchange,
		    scalar_type* ec, scalar_type* correlation,
		    scalar_type* vrho, scalar_type* vrhoC,
		    scalar_type* vsigma, scalar_type* vsigmaC,
		    scalar_type* v2rho, scalar_type* v2rhoC,
		    scalar_type* v2rhosigma, scalar_type* v2rhosigmaC,
		    scalar_type* v2sigma, scalar_type* v2sigmaC,
		    scalar_type* y2a,
		    const scalar_type* sigma,
		    const G2G::vec_type<scalar_type, width>* grad,
		    const G2G::vec_type<scalar_type, width>* hess1,
		    const G2G::vec_type<scalar_type, width>* hess2,
		    int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
	ex[i] = exchange[i];
	ec[i] = correlation[i];
	// Merge the results for the derivatives.
	vrho[i] += vrhoC[i];
        vsigma[i] += vsigmaC[i];
        v2rho[i] += v2rhoC[i];
        v2rhosigma[i] += v2rhosigmaC[i];
        v2sigma[i] += v2sigmaC[i];
        // Now, compute y2a value.
	y2a[i] = vrho[i] - (2 * sigma[i] * v2rhosigma[i]
            + 2 * (hess1[i].x + hess1[i].y + hess1[i].z) * vsigma[i]
            + 4 * v2sigma[i] * (grad[i].x * grad[i].x * hess1[i].x + grad[i].y * grad[i].y * hess1[i].y + grad[i].z * grad[i].z * hess1[i].z + 2 * grad[i].x * grad[i].y * hess2[i].x + 2 * grad[i].x * grad[i].z * hess2[i].y + 2 * grad[i].y * grad[i].z * hess2[i].z));
    }
}


template <class scalar_type, int width>
void LibxcProxy <scalar_type, width>::doGGA(scalar_type* dens,
    const int number_of_points,
    const scalar_type* contracted_grad,
    const G2G::vec_type<scalar_type, width>* grad,
    const G2G::vec_type<scalar_type, width>* hess1,
    const G2G::vec_type<scalar_type, width>* hess2,
    scalar_type* ex,
    scalar_type* ec,
    scalar_type* y2a)
{

    printf("LibxcProxy::doGGA cuda (...) \n");

    double* rho = dens;
    const double* sigma = contracted_grad;

    cudaError_t err = cudaSuccess;

    int array_size = sizeof(double) * number_of_points;
    double* exchange = NULL;
    err = cudaMalloc((void **)&exchange, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device exchange!\n");
        exit(EXIT_FAILURE);
    }

    double* correlation = NULL;
    err = cudaMalloc((void **)&correlation, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device correlation!\n");
        exit(EXIT_FAILURE);
    }

    // The outputs for exchange
    double* vrho = NULL;
    double* vsigma = NULL;
    double* v2rho = NULL;
    double* v2rhosigma = NULL;
    double* v2sigma = NULL;

    err = cudaMalloc((void **)&vrho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigma!\n");
        exit(EXIT_FAILURE);
    }

    // The outputs for correlation
    double* vrhoC = NULL;
    double* vsigmaC = NULL;
    double* v2rhoC = NULL;
    double* v2rhosigmaC = NULL;
    double* v2sigmaC = NULL;

    err = cudaMalloc((void **)&vrhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigmaC!\n");
        exit(EXIT_FAILURE);
    }

    try {
        xc_gga (&funcForExchange, number_of_points,
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
        fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
        return;
    }

    try {
        // Now the correlation value.
        xc_gga (&funcForCorrelation, number_of_points,
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
        fprintf (stderr, "Exception ocurred calling xc_gga for Correlation '%d' \n", exception);
        return;
    }

    // Gather the results.
    // Variables for the Kernels
    int threadsPerBlock = 256;
    int blocksPerGrid = (number_of_points + threadsPerBlock - 1) / threadsPerBlock;
    joinResults<scalar_type, width><<<blocksPerGrid, threadsPerBlock>>>(ex, exchange,
	ec, correlation,
	vrho, vrhoC,
	vsigma, vsigmaC,
	v2rho, v2rhoC,
	v2rhosigma, v2rhosigmaC,
	v2sigma, v2sigmaC,
	y2a,
	sigma,
	grad,
	hess1,
	hess2, 
	number_of_points);

    // Free device memory.
    if (exchange != NULL) {
	cudaFree(exchange);
    }
    if (correlation != NULL) {
	cudaFree(correlation);
    }
    if (vrho != NULL) {
        cudaFree(vrho);
    }
    if (vsigma != NULL) {
	cudaFree(vsigma);
    }
    if (v2rho != NULL) {
	cudaFree(v2rho);
    }
    if (v2rhosigma != NULL) {
	cudaFree(v2rhosigma);
    }
    if (v2sigma != NULL) {
	cudaFree(v2sigma);
    }
    if (vrhoC != NULL) {
        cudaFree(vrhoC);
    }
    if (vsigmaC != NULL) {
	cudaFree(vsigmaC);
    }
    if (v2rhoC != NULL) {
	cudaFree(v2rhoC);
    }
    if (v2rhosigmaC != NULL) {
	cudaFree(v2rhosigmaC);
    }
    if (v2sigmaC != NULL) {
	cudaFree(v2sigmaC);
    }

    return;
}


template <class scalar_type, int width>
void LibxcProxy <scalar_type, width>::doLDA(scalar_type dens, const G2G::vec_type<scalar_type, width> &grad, const G2G::vec_type<scalar_type, width> &hess1, const G2G::vec_type<scalar_type, width> &hess2, scalar_type &ex, scalar_type &ec, scalar_type &y2a)
{
    //TODO: not implemented yet!
    return;
}

#endif // LIBXCPROXY_H
