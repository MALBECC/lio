#include <iostream>
#include <stdexcept>
#include <cuda_runtime.h>
#include <cassert>
#include "matrix.h"
using namespace G2G;
using namespace std;

/*
 * Matrix
 */

template<class T> Matrix<T>::Matrix(void) : data(NULL), width(0), height(0) /*, components(0)*/ {}

template<class T> Matrix<T>::~Matrix(void) { } 

template<class T> unsigned int Matrix<T>::bytes(void) const {
	return elements() * sizeof(T);
}

template<class T> unsigned int Matrix<T>::elements(void) const {
	return width * height /* * components */;
}

/*
 * HostMatrix
 */

template<class T> HostMatrix<T>::HostMatrix(void) : Matrix<T>() { }

/*HostMatrix::HostMatrix(const HostMatrix& c) : Matrix() {
	*this = c;
}*/

template<class T> HostMatrix<T>::HostMatrix(unsigned int _width, unsigned _height) : Matrix<T>() {
	resize(_width, _height);
}

template<class T> HostMatrix<T>::HostMatrix(const CudaMatrix<T>& c) : Matrix<T>() {
	*this = c;
}

template<class T> HostMatrix<T>::~HostMatrix(void) {
	delete[] this->data;	
}

template<class T> HostMatrix<T>& HostMatrix<T>::resize(unsigned int _width, unsigned _height) {
	if (this->data) delete[] this->data;
	this->width = _width; this->height = _height;
	this->data = new T[this->elements()];
	return *this;
}

template<class T> HostMatrix<T>& HostMatrix<T>::fill(const T& value) {
	for (uint i = 0; i < this->elements(); i++) this->data[i] = value;
	return *this;
}

/*const HostMatrix& HostMatrix::operator=(const HostMatrix& c) {
	if (!c.data) delete[] data;
	else {
		if (data) {
			if (bytes != c.bytes()) { delete[] data; width = c.width; height = c.height; data = new float[this->elements()]; }
		}
		else { width = c.width; height = c.height; data = new float[this->elements()]; }
		
		copy_submatrix(c);
	}
	return *this;
}*/

template <class T> HostMatrix<T>& HostMatrix<T>::operator=(const CudaMatrix<T>& c) {
	if (!c.data) {
		if (this->data) { delete[] this->data; this->width = this->height = 0; this->data = NULL; }
	}
	else {
		if (this->data) {
			if (this->bytes() != c.bytes()) {
				delete[] this->data;
				this->width = c.width; this->height = c.height;
				this->data = new T[c.elements()];
			}			
		}
		else {
			this->width = c.width; this->height = c.height;
			this->data = new T[c.elements()];
		}

		copy_submatrix(c);
	}

	return *this;		
}

template<class T> void HostMatrix<T>::copy_submatrix(const CudaMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes()) throw runtime_error("Can't copy more elements than what operator has");
	
	cudaMemcpy(this->data, c.data, _bytes, cudaMemcpyDeviceToHost);
}

/*
 * CudaMatrix
 */

template<class T> CudaMatrix<T>::CudaMatrix(void) : Matrix<T>() {
	assert(this->data == NULL);
}

template<class T> CudaMatrix<T>::CudaMatrix(unsigned int _width, unsigned int _height) : Matrix<T>() {
	assert(this->data == NULL);
	resize(_width, _height);
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::resize(unsigned int _width, unsigned int _height) {
	if (this->data) cudaFree(this->data);
	this->width = _width; this->height = _height;
	cudaError_t error_status = cudaMalloc((void**)&this->data, this->bytes());
	assert(error_status != cudaErrorMemoryAllocation);
	return *this;		
}

template<class T> CudaMatrix<T>::CudaMatrix(const CudaMatrix<T>& c) : Matrix<T>() {
	assert(this->data == NULL);	
	*this = c;
}

template<class T> CudaMatrix<T>::CudaMatrix(const HostMatrix<T>& c) : Matrix<T>() {
	assert(this->data == NULL);	
	*this = c;
}

template<class T> CudaMatrix<T>::~CudaMatrix(void) {
	if (this->data) cudaFree(this->data);
}

template<class T> void CudaMatrix<T>::copy_submatrix(const HostMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes()) throw runtime_error("CudaMatrix: Can't copy more elements than what operand has");

	cudaMemcpy(this->data, c.data, _bytes, cudaMemcpyHostToDevice);
}

template<class T> void CudaMatrix<T>::copy_submatrix(const CudaMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	if (_bytes > c.bytes()) throw runtime_error("CudaMatrix: Can't copy more elements than what operand has");
	cudaMemcpy(c.data, this->data, _bytes, cudaMemcpyDeviceToDevice);
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::operator=(const HostMatrix<T>& c) {
	if (!c.data) {
		if (this->data) { cudaFree(this->data); this->width = this->height = 0; this->data = NULL; }
	}
	else {
		if (this->data) {
			if (this->bytes() != c.bytes()) {
				cudaFree(this->data);
				this->width = c.width; this->height = c.height;
				cudaMalloc((void**)&this->data, this->bytes());
			}			
		}
		else {
			this->width = c.width; this->height = c.height;
			cudaMalloc((void**)&this->data, this->bytes());
		}
		
		copy_submatrix(c);
	}

	return *this;
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::operator=(const CudaMatrix<T>& c) {
	// copies data from c, only if necessary (always frees this's data, if any)
	if (!c.data) {
		if (this->data) { cudaFree(this->data); this->width = this->height = 0; this->data = NULL; }
	}
	else {
		if (this->data) {
			if (this->bytes() != c.bytes()) {
				cudaFree(this->data);
				this->width = c.width; this->height = c.height;
				cudaMalloc((void**)&this->data, this->bytes());
			}
		}
		else {
			this->width = c.width; this->height = c.height;
			cudaMalloc((void**)&this->data, this->bytes());
		}
		
		cudaMemcpy(this->data, c.data, this->bytes(), cudaMemcpyDeviceToDevice);		
	}
	
	return *this;
}

/**
 * Instantiations
 */
template class Matrix<double>;
template class Matrix<float>;
template class Matrix<float1>;
template class Matrix<float2>;
template class Matrix<float3>;
template class Matrix<float4>;
template class Matrix<uint1>;
template class Matrix<uint>;

template class HostMatrix<double>;
template class HostMatrix<float>;
template class HostMatrix<float1>;
template class HostMatrix<float2>;
template class HostMatrix<float3>;
template class HostMatrix<float4>;
template class HostMatrix<uint1>;
template class HostMatrix<uint>;

template class CudaMatrix<float>;
template class CudaMatrix<float1>;
template class CudaMatrix<float2>;
template class CudaMatrix<float3>;
template class CudaMatrix<float4>;
template class CudaMatrix<uint>;
