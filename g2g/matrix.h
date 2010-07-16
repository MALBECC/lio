#ifndef __G2G_MATRIX_H__
#define __G2G_MATRIX_H__

#include <vector>
#include "cpu/cpu_vector_types.h"
#define COALESCED_DIMENSION(d) (d + 16 - (d % 16))

#include "cuda/cuda_extra.h"

namespace G2G {
	template<class T> class Matrix {
		public:
			Matrix(void);
			virtual ~Matrix(void);				

			T* data;
			unsigned int width;
			unsigned int height;
			// unsigned int components;
	
			unsigned int bytes(void) const;
			unsigned int elements(void) const;
		
			bool is_allocated(void) const;

      virtual void deallocate(void) = 0;
	};
	
	template<class T> class CudaMatrix;
	
	template<class T> class HostMatrix : public Matrix<T> {
		public:
      enum PinnedFlag { Pinned, NonPinned };

      HostMatrix(PinnedFlag pinned = NonPinned);
			HostMatrix(unsigned int width, unsigned int height = 1, PinnedFlag pinned = NonPinned);
			HostMatrix(const CudaMatrix<T>& c);
			HostMatrix(const HostMatrix<T>& c);
			~HostMatrix(void);
		
			HostMatrix<T>& operator=(const CudaMatrix<T>& c);
			HostMatrix<T>& operator=(const HostMatrix<T>& c);

		  // TODO: mejorar esto
			inline const T& operator()(unsigned int i = 0, unsigned int j = 0) const {
				assert(i < this->width);
				assert(j < this->height);
        return this->data[j * this->width + i];
			}
			inline T& operator()(unsigned int i = 0, unsigned int j = 0) {
				assert(i < this->width);
				assert(j < this->height);
				return this->data[j * this->width + i];
			}
			
			void copy_submatrix(const CudaMatrix<T>& c, unsigned int elements = 0);
			void copy_submatrix(const HostMatrix<T>& c, unsigned int elements = 0);

      void copy_transpose(const CudaMatrix<T>& cuda_matrix);
      void transpose(HostMatrix<T>& out);
		
			HostMatrix<T>& resize(unsigned int width, unsigned int height = 1);
			HostMatrix<T>& zero(void);
      HostMatrix<T>& fill(T value);

			void deallocate(void);

      /* BLAS methods */
      enum UpperLowerTriangle { UpperTriangle, LowerTriangle };
      static void blas_ssyr(UpperLowerTriangle triangle, float alpha, const HostMatrix<T>& x, const HostMatrix<T>& A, unsigned int x_row);

			void to_constant(const char* constant);
		
		private:
			bool pinned;
			void alloc_data(void);
			void dealloc_data(void);
	};

  template<class T> void to_constant(const char* constant, const T& value);
	
	template <class T> class CudaMatrix : public Matrix<T> {
		public:
			CudaMatrix(void);
			CudaMatrix(const CudaMatrix<T>& c);
			CudaMatrix(const HostMatrix<T>& c);
      CudaMatrix(const std::vector<T>& v);

			CudaMatrix(unsigned int width, unsigned int height = 1);
			~CudaMatrix(void);
		
			CudaMatrix<T>& resize(unsigned int width, unsigned int height = 1);
			CudaMatrix<T>& zero(void);

			CudaMatrix& operator=(const HostMatrix<T>& c);
			CudaMatrix& operator=(const CudaMatrix<T>& c);
      CudaMatrix& operator=(const std::vector<T>& v);
		
			void copy_submatrix(const HostMatrix<T>& c, unsigned int elements = 0);
			void copy_submatrix(const CudaMatrix<T>& c, unsigned int elements = 0);
      void copy_submatrix(const std::vector<T>& v, unsigned int elements = 0);
		
			void deallocate(void);		
	};
	
	template <class T> class FortranMatrix : public Matrix<T> {
		public:
			FortranMatrix(void);					
			FortranMatrix(T* ptr, unsigned int width, unsigned int height = 1, unsigned int fortran_width = 1);
		
			inline T& operator()(unsigned int x = 0, unsigned int y = 0) {
				assert(x < this->width);
				assert(y < this->height);
				return this->data[y * fortran_width + x];
			}
			inline const T& operator()(unsigned int x = 0, unsigned int y = 0) const {
				assert(x < this->width);
				assert(y < this->height);
				return this->data[y * fortran_width + x];
			}

      void deallocate(void) {};
		
		private:
			unsigned int fortran_width;
	};
	
	typedef HostMatrix<double> HostMatrixDouble;
	typedef HostMatrix<double3> HostMatrixDouble3;
	typedef HostMatrix<float> HostMatrixFloat;
	typedef HostMatrix<float1> HostMatrixFloat1;
	typedef HostMatrix<float2> HostMatrixFloat2;
	typedef HostMatrix<float3> HostMatrixFloat3;
  typedef HostMatrix<float4> HostMatrixFloat4;
	typedef HostMatrix<uint> HostMatrixUInt;
	typedef HostMatrix<uint1> HostMatrixUInt1;
  typedef HostMatrix<uint2> HostMatrixUInt2;

  #ifndef __CUDACC__
  typedef HostMatrix<cfloat3> HostMatrixCFloat3;
	typedef HostMatrix<cfloat4> HostMatrixCFloat4;
  #endif
	
	typedef CudaMatrix<float> CudaMatrixFloat;
	typedef CudaMatrix<float1> CudaMatrixFloat1;
	typedef CudaMatrix<float2> CudaMatrixFloat2;
	typedef CudaMatrix<float3> CudaMatrixFloat3;
	typedef CudaMatrix<float4> CudaMatrixFloat4;	
	typedef CudaMatrix<uint> CudaMatrixUInt;
	typedef CudaMatrix<uint1> CudaMatrixUInt1;
  typedef CudaMatrix<uint2> CudaMatrixUInt2;
}

#endif
