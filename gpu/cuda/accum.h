#ifndef __G2G_ACCUM_H__
#define __G2G_ACCUM_H__

extern "C" void calc_accum(const G2G::HostMatrixFloat& input, G2G::HostMatrixFloat& output);
extern "C" void calc_accum_cuda(const G2G::CudaMatrixFloat& input, G2G::CudaMatrixFloat& output);		// TODO: realmente es necesario el extern?

#endif
