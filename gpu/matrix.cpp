#include <iostream>
#include <stdexcept>
#include <cuda_runtime.h>
#include "matrix.h"
using namespace G2G;
using namespace std;

/*
 * Matrix
 */

Matrix::Matrix(void) : data(NULL), width(0), height(0) /*, components(0)*/ {}

Matrix::~Matrix(void) { }

unsigned int Matrix::bytes(void) const {
	return elements() * sizeof(float);	
}

unsigned int Matrix::elements(void) const {
	return width * height /* * components */;
}

/*
 * HostMatrix
 */

HostMatrix::HostMatrix(void) : Matrix() { }

HostMatrix::HostMatrix(float elem) : Matrix() {
	width = height = 1;
	data = new float[1];
	data[0] = elem;
}

HostMatrix::HostMatrix(unsigned int _width, unsigned _height) : Matrix() {
	width = _width; height = _height;
	data = new float[elements()];
}

HostMatrix::~HostMatrix(void) {
	delete[] data;	
}

void HostMatrix::copy_submatrix(const CudaMatrix& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? bytes() : _elements * sizeof(float));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes()) throw runtime_error("Can't copy more elements than what operator has");
	
	cudaMemcpy(data, c.data, _bytes, cudaMemcpyDeviceToHost);
}

/*
 * CudaMatrix
 */

CudaMatrix::CudaMatrix(void) : Matrix() { }

CudaMatrix::CudaMatrix(unsigned int _width, unsigned int _height) : Matrix() {
	width = _width; height = _height;
 	cudaMalloc((void**)&data, bytes());
}

CudaMatrix::CudaMatrix(const CudaMatrix& c) : Matrix() {
	*this = c;
}

CudaMatrix::CudaMatrix(const HostMatrix& c) : Matrix() {
	*this = c;
}

CudaMatrix::~CudaMatrix(void) {
	if (data) cudaFree(data);	
}

void CudaMatrix::copy_submatrix(const HostMatrix& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? bytes() : _elements * sizeof(float));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes()) throw runtime_error("CudaMatrix: Can't copy more elements than what operator has");

	cudaMemcpy(data, c.data, _bytes, cudaMemcpyHostToDevice);
}

CudaMatrix& CudaMatrix::operator=(const HostMatrix& c) {
	if (!c.data) {
		if (data) { cudaFree(data); width = height = 0; }
	}
	else {
		if (data) {
			if (bytes() != c.bytes()) {
				cudaFree(data);
				width = c.width; height = c.height;
				cudaMalloc((void**)&data, bytes());
			}			
		}
		else {
			width = c.width; height = c.height;
			cudaMalloc((void**)&data, bytes());
		}
		
		copy_submatrix(c);
	}
	
	return *this;
}

CudaMatrix& CudaMatrix::operator=(const CudaMatrix& c) {
	// copies data from c, only if necessary (always frees this's data, if any)
	if (!c.data) {
		if (data) { cudaFree(data); width = height = 0; }
	}
	else {
		if (data) {
			if (bytes() != c.bytes()) {
				cudaFree(data);
				width = c.width; height = c.height;
				cudaMalloc((void**)&data, bytes());
			}
		}
		else {
			width = c.width; height = c.height;
			cudaMalloc((void**)&data, bytes());
		}
		
		cudaMemcpy(data, c.data, bytes(), cudaMemcpyDeviceToDevice);		
	}
	
	return *this;
}
