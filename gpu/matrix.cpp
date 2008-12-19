#include <iostream>
#include <stdexcept>
#include <cassert>
#include <cuda_runtime.h>
#include "matrix.h"
using namespace G2G;
using namespace std;

/***************************
 * Matrix
 ***************************/

template<class T> Matrix<T>::Matrix(void) : data(NULL), width(0), height(0) /*, components(0)*/ {}

template<class T> Matrix<T>::~Matrix(void) { } 

template<class T> unsigned int Matrix<T>::bytes(void) const {
	return elements() * sizeof(T);
}

template<class T> unsigned int Matrix<T>::elements(void) const {
	return width * height /* * components */;
}

template<class T> bool Matrix<T>::is_allocated(void) const {
	return data;
}

/***************************
 * HostMatrix
 ***************************/
template<class T> void HostMatrix<T>::alloc_data(void) {
	if (pinned) {
		cudaError_t error_status = cudaMallocHost((void**)&this->data, this->bytes());
		assert(error_status != cudaErrorMemoryAllocation);
	}	
	else this->data = new T[this->elements()];
	
	assert(this->data);
}

template<class T> void HostMatrix<T>::dealloc_data(void) {
	if (pinned) cudaFreeHost(this->data);
	else delete[] this->data;
}

template<class T> void HostMatrix<T>::deallocate(void) {
	dealloc_data();
	this->data = NULL;	
}

template<class T> HostMatrix<T>::HostMatrix(void) : Matrix<T>(), pinned(false) { }

//template<class T> HostMatrix<T>::HostMatrix(bool _pinned) : Matrix<T>(), pinned(_pinned) { }

/*HostMatrix::HostMatrix(const HostMatrix& c) : Matrix() {
	*this = c;
}*/

template<class T> HostMatrix<T>::HostMatrix(unsigned int _width, unsigned _height) : Matrix<T>(), pinned(false) {
	resize(_width, _height);
}

template<class T> HostMatrix<T>::HostMatrix(const CudaMatrix<T>& c) : Matrix<T>(), pinned(false) {
	*this = c;
}

template<class T> HostMatrix<T>::HostMatrix(const HostMatrix<T>& m) : Matrix<T>(), pinned(false) {
	*this = m;
}

template<class T> HostMatrix<T>::~HostMatrix(void) {
	deallocate();
}

template<class T> HostMatrix<T>& HostMatrix<T>::resize(unsigned int _width, unsigned _height) {
	if (this->data) dealloc_data();
	this->width = _width; this->height = _height;	
	alloc_data();
	
	return *this;
}

template<class T> HostMatrix<T>& HostMatrix<T>::fill(const T& value) {
	for (uint i = 0; i < this->elements(); i++) this->data[i] = value;
	return *this;
}

template<class T> HostMatrix<T>& HostMatrix<T>::operator=(const HostMatrix<T>& c) {
	assert(!this->pinned);
	
	if (!c.data) {
		if (this->data) { dealloc_data(); this->width = this->height = 0; this->data = NULL; }
	}
	else {
		if (this->data) {
			if (this->bytes() != c.bytes()) {
				dealloc_data();
				this->width = c.width; this->height = c.height; 
				alloc_data();
			}
		}
		else {
			this->width = c.width; this->height = c.height;
			alloc_data();
		}
		
		copy_submatrix(c);
	}
	
	return *this;
}

template <class T> HostMatrix<T>& HostMatrix<T>::operator=(const CudaMatrix<T>& c) {
	if (!c.data) {
		if (this->data) { dealloc_data(); this->width = this->height = 0; this->data = NULL; }
	}
	else {
		if (this->data) {
			if (this->bytes() != c.bytes()) {
				dealloc_data();
				this->width = c.width; this->height = c.height;
				alloc_data();
			}			
		}
		else {
			this->width = c.width; this->height = c.height;
			alloc_data();
		}

		copy_submatrix(c);
	}

	return *this;		
}

template<class T> void HostMatrix<T>::copy_submatrix(const HostMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes()) throw runtime_error("Can't copy more elements than what operator has");
	memcpy(this->data, c.data, _bytes);
}

template<class T> void HostMatrix<T>::copy_submatrix(const CudaMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes()) throw runtime_error("Can't copy more elements than what operator has");
	
	cudaError_t ret = cudaMemcpy(this->data, c.data, _bytes, cudaMemcpyDeviceToHost);
	assert(ret == cudaSuccess);
}

template<class T> void HostMatrix<T>::to_constant(const char* symbol) {
	cudaMemcpyToSymbol(symbol, this->data, this->bytes(), 0, cudaMemcpyHostToDevice);	
}

/*template<class T, class S> void HostMatrix<T>::to_constant<S>(const char* constant, const S& value) {
	cudaMemcpyToSymbol(constant, &value, sizeof(S), 0, cudaMemcpyHostToDevice);	
}*/

/******************************
 * CudaMatrix
 ******************************/

template<class T> CudaMatrix<T>::CudaMatrix(void) : Matrix<T>() { }

template<class T> CudaMatrix<T>::CudaMatrix(unsigned int _width, unsigned int _height) : Matrix<T>() {
	resize(_width, _height);
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::resize(unsigned int _width, unsigned int _height) {
	if (this->data) cudaFree(this->data);
	this->width = _width; this->height = _height;
	cudaError_t error_status = cudaMalloc((void**)&this->data, this->bytes());
	assert(error_status != cudaErrorMemoryAllocation);
	return *this;		
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::fill(int value) {
	assert(this->data);
	cudaMemset(this->data, value, this->bytes());	
	return *this;
}

template<class T> CudaMatrix<T>::CudaMatrix(const CudaMatrix<T>& c) : Matrix<T>() {
	*this = c;
}

template<class T> CudaMatrix<T>::CudaMatrix(const HostMatrix<T>& c) : Matrix<T>() {
	*this = c;
}

template<class T> CudaMatrix<T>::~CudaMatrix(void) {
	if (this->data) cudaFree(this->data);
}

template<class T> void CudaMatrix<T>::deallocate(void) {
	if (this->data) cudaFree(this->data);	
	this->data = NULL;
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
				cudaError_t error_status = cudaMalloc((void**)&this->data, this->bytes());
				assert(error_status != cudaErrorMemoryAllocation);
			}			
		}
		else {
			this->width = c.width; this->height = c.height;
			cudaError_t error_status = cudaMalloc((void**)&this->data, this->bytes());
			assert(error_status != cudaErrorMemoryAllocation);			
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
				cudaError_t error_status = cudaMalloc((void**)&this->data, this->bytes());
				assert(error_status != cudaErrorMemoryAllocation);				
			}
		}
		else {
			this->width = c.width; this->height = c.height;
			cudaError_t error_status = cudaMalloc((void**)&this->data, this->bytes());
			assert(error_status != cudaErrorMemoryAllocation);			
		}
		
		cudaMemcpy(this->data, c.data, this->bytes(), cudaMemcpyDeviceToDevice);		
	}
	
	return *this;
}

/*************************************
 * FortranMatrix
 *************************************/
template<class T> FortranMatrix<T>::FortranMatrix(void)
	: Matrix<T>(), fortran_width(0)
{ }

template<class T> FortranMatrix<T>::FortranMatrix(T* _data, unsigned int _width, unsigned int _height, unsigned int _fortran_width)
	: Matrix<T>(), fortran_width(_fortran_width)
{
	this->data = _data;
	this->width = _width; this->height = _height;
	assert(this->data);
}

/**
 * Instantiations
 */
template class Matrix<double3>;
template class Matrix<double>;
template class Matrix<float>;
template class Matrix<float1>;
template class Matrix<float2>;
template class Matrix<float3>;
template class Matrix<float4>;
template class Matrix<uint1>;
template class Matrix<uint>;

template class HostMatrix<double3>;
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
template class CudaMatrix<double>;

template class FortranMatrix<double>;
template class FortranMatrix<uint>;
