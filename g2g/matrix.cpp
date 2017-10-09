#include <iostream>
#include <stdexcept>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "common.h"
#include "matrix.h"
#include "scalar_vector_types.h"
using namespace std;

namespace G2G {

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
  assert(this->bytes() != 0);
  int bytes = this->bytes();
  int posix_return = 0;

  if (pinned) {
#if GPU_KERNELS
    cudaError_t error_status = cudaMallocHost((void**)&this->data, this->bytes());
    assert(error_status == cudaSuccess);
#else
    assert(false);
#endif
  }
  else
  {
    posix_return = posix_memalign((void **) &this->data, 64, this->bytes());
    if ( posix_return != 0) 
    { 
       std::cout <<"HostMatrix: Error in posix_memalign.\n"; 
       exit(1);
    };
  };

  assert(this->data);
}

template<class T> void HostMatrix<T>::dealloc_data(void) {
	if (pinned) {
    #if GPU_KERNELS
    cudaFreeHost(this->data);
    #else
    assert(false);
    #endif
  }
	else free(this->data); //mkl_free(this->data);
}

template<class T> void HostMatrix<T>::copy_to_tmp(T * dst) const {
    memcpy(dst, this->data, this->bytes());
}

template<class T> void HostMatrix<T>::deallocate(void) {
	dealloc_data();
	this->data = NULL;
  this->width = this->height = 0;
}

template<class T> HostMatrix<T>::HostMatrix(PinnedFlag _pinned) : Matrix<T>() {
  pinned = (_pinned == Pinned);
}

template<class T> HostMatrix<T>::HostMatrix(unsigned int _width, unsigned _height, PinnedFlag _pinned) : Matrix<T>() {
  pinned = (_pinned == Pinned);
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
  if (_width == 0 ) throw std::runtime_error("El ancho no puede ser 0");
  if (_height == 0 ) throw std::runtime_error("La altura no puede ser 0");
  if (_width != this->width || _height != this->height) {
    if (this->data) dealloc_data();
    this->width = _width; this->height = _height;
    alloc_data();
  }

	return *this;
}

template<class T> HostMatrix<T>& HostMatrix<T>::shrink(unsigned int _width, unsigned int _height) {
  if (_width == 0 || _height == 0) throw std::runtime_error("La dimension no puede ser 0");
  if (_width != this->width || _height != this->height) {
    HostMatrix<T> temp_matrix(_width, _height);
    temp_matrix.copy_submatrix(temp_matrix, _width * _height);
    resize(_width, _height);
    *this = temp_matrix;
  }

	return *this;
}

template<class T> HostMatrix<T>& HostMatrix<T>::zero(void) {
  memset(this->data, 0, this->bytes());
	return *this;
}

template<class T> HostMatrix<T>& HostMatrix<T>::fill(T value) {
  for (uint i = 0; i < this->elements(); i++) { this->data[i] = value; }
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
	if (_bytes > c.bytes())
    throw runtime_error("Can't copy more elements than what operator has");
	memcpy(this->data, c.data, _bytes);
}

template<class T> void HostMatrix<T>::copy_submatrix(const CudaMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes())
    throw runtime_error("Can't copy more elements than what operator has");

  #if GPU_KERNELS
	cudaMemcpy(this->data, c.data, _bytes, cudaMemcpyDeviceToHost);
  cudaAssertNoError("HostMatrix::copy_submatrix");
  #else
  assert(false);
  #endif
}

template<class T> void HostMatrix<T>::to_constant(const char* symbol) {
  #if GPU_KERNELS
	cudaMemcpyToSymbol(symbol, this->data, this->bytes(), 0, cudaMemcpyHostToDevice);
  cudaAssertNoError("to_constant");
  #endif
}

template<class T> void HostMatrix<T>::transpose(HostMatrix<T>& out) const {
  out.resize(this->height, this->width);
  out.zero();
  for (uint i = 0; i < this->width; i++) {
    for (uint j = 0; j < this->height; j++) {
      out(j, i) = (*this)(i, j);
    }
  }
}
template<class T> void HostMatrix<T>::copy_transpose(const CudaMatrix<T>& cuda_matrix) {
  if (cuda_matrix.width != this->height || cuda_matrix.height != this->width) throw runtime_error("Matrix dimensions for copy_transpose don't agree");
  HostMatrix<T> cuda_matrix_copy(cuda_matrix);
  for (uint i = 0; i < cuda_matrix.width; i++) {
    for (uint j = 0; j < cuda_matrix.height; j++) {
      (*this)(j, i) = cuda_matrix_copy(i, j);
    }
  }
}

template<class T> void to_constant(const char* constant, const T& value) {
  #if GPU_KERNELS
	cudaMemcpyToSymbol(constant, &value, sizeof(T), 0, cudaMemcpyHostToDevice);
  cudaAssertNoError("to_constant(value)");
  #endif
}

template void to_constant<uint>(const char* constant, const uint& value);
template void to_constant<float>(const char* constant, const float& value);
template void to_constant<double>(const char* constant, const double& value);

/******************************
 * CudaMatrix
 ******************************/

template<class T> CudaMatrix<T>::CudaMatrix(void) : Matrix<T>() { }

template<class T> CudaMatrix<T>::CudaMatrix(unsigned int _width, unsigned int _height) : Matrix<T>() {
	resize(_width, _height);
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::resize(unsigned int _width, unsigned int _height) {
  assert(_width * _height != 0);

  #if GPU_KERNELS
  if (_width != this->width || _height != this->height) {
    if (this->data) cudaFree(this->data);
    this->width = _width; this->height = _height;
    cudaMalloc((void**)&this->data, this->bytes());
    cudaAssertNoError("CudaMatrix::resize");
  }
  #endif
	return *this;
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::zero(void) {
  #if GPU_KERNELS
	assert(this->data);
	cudaMemset(this->data, 0, this->bytes());
  cudaAssertNoError("CudaMatrix::zero");
  #endif
	return *this;
}

template<class T> CudaMatrix<T>::CudaMatrix(const CudaMatrix<T>& c) : Matrix<T>() {
	*this = c;
}

template<class T> CudaMatrix<T>::CudaMatrix(const HostMatrix<T>& c) : Matrix<T>() {
	*this = c;
}

template<class T> CudaMatrix<T>::CudaMatrix(const std::vector<T>& v) : Matrix<T>() {
	*this = v;
}

template<class T> CudaMatrix<T>::~CudaMatrix(void) {
  deallocate();
}

template<class T> void CudaMatrix<T>::deallocate(void) {
  #if GPU_KERNELS
	if (this->data) cudaFree(this->data);
	this->data = NULL;
  this->width = this->height = 0;
  #endif
}

template<class T> void CudaMatrix<T>::copy_submatrix(const HostMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	//cout << "bytes: " << _bytes << ", c.bytes: " << c.bytes() << endl;
	if (_bytes > c.bytes()) throw runtime_error("CudaMatrix: Can't copy more elements than what operand has");

  #if GPU_KERNELS
	cudaMemcpy(this->data, c.data, _bytes, cudaMemcpyHostToDevice);
  cudaAssertNoError("CudaMatrix::copy_submatrix");
  #endif
}

template<class T> void CudaMatrix<T>::copy_submatrix(const CudaMatrix<T>& c, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	if (_bytes > c.bytes()) throw runtime_error("CudaMatrix: Can't copy more elements than what operand has");

  #if GPU_KERNELS
	cudaMemcpy(c.data, this->data, _bytes, cudaMemcpyDeviceToDevice);
  cudaAssertNoError("CudaMatrix::copy_submatrix");
  #endif
}

template<class T> void CudaMatrix<T>::copy_submatrix(const std::vector<T>& v, unsigned int _elements) {
	unsigned int _bytes = (_elements == 0 ? this->bytes() : _elements * sizeof(T));
	if (_bytes > v.size() * sizeof(T)) throw runtime_error("CudaMatrix: Can't copy more elements than what operand has");

  #if GPU_KERNELS
	cudaMemcpy(this->data, (T*)&v[0], _bytes, cudaMemcpyHostToDevice);
  cudaAssertNoError("CudaMatrix::copy_submatrix");
  #endif
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::operator=(const HostMatrix<T>& c) {
  #if GPU_KERNELS
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
		cudaAssertNoError("CudaMatrix::operator=");
		copy_submatrix(c);
	}
  #endif
	return *this;
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::operator=(const std::vector<T>& v) {
  #if GPU_KERNELS
	if (v.empty()) {
		if (this->data) { cudaFree(this->data); this->width = this->height = 0; this->data = NULL; }
	}
	else {
		if (this->data) {
			if (this->elements() != v.size()) {
				cudaFree(this->data);
				this->width = v.size(); this->height = 1;
				cudaMalloc((void**)&this->data, this->bytes());
			}
		}
		else {
			this->width = v.size(); this->height = 1;
			cudaMalloc((void**)&this->data, this->bytes());
		}
    cudaAssertNoError("CudaMatrix::operator=");
		copy_submatrix(v);
	}
  #endif
	return *this;
}

template<class T> CudaMatrix<T>& CudaMatrix<T>::operator=(const CudaMatrix<T>& c) {
  #if GPU_KERNELS
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
  cudaAssertNoError("CudaMatrix::operator=");
	#endif
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
template class Matrix<double>;
template class Matrix<double3>;
template class Matrix<float>;
template class Matrix<float3>;
template class Matrix<uint>;

template class Matrix< vec_type<float, 2> >;
template class Matrix< vec_type<float, 3> >;
template class Matrix< vec_type<double, 2> >;
template class Matrix< vec_type<double, 3> >;

template class HostMatrix< vec_type<float, 2> >;
template class HostMatrix< vec_type<float, 3> >;
template class HostMatrix< vec_type<float, 4> >;
template class HostMatrix< vec_type<double, 2> >;
template class HostMatrix< vec_type<double, 3> >;
template class HostMatrix< vec_type<double, 4> >;
template class CudaMatrix< vec_type<float, 2> >;
template class CudaMatrix< vec_type<float, 3> >;
template class CudaMatrix< vec_type<float, 4> >;
template class CudaMatrix< vec_type<double, 2> >;
template class CudaMatrix< vec_type<double, 3> >;
template class CudaMatrix< vec_type<double, 4> >;

template class HostMatrix<double>;
template class HostMatrix<float>;

template class HostMatrix<double3>;
template class HostMatrix<float3>;
template class HostMatrix<uint>;

template class CudaMatrix<float>;
template class CudaMatrix<uint>;
template class CudaMatrix<double>;

template class FortranMatrix<double>;
template class FortranMatrix<unsigned int>;

}
