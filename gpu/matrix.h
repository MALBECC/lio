#ifndef __G2G_MATRIX_H__
#define __G2G_MATRIX_H__

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
	};
	
	template<class T> class CudaMatrix;
	
	template<class T> class HostMatrix : public Matrix<T> {
		public:
			HostMatrix(void);
			HostMatrix(unsigned int width, unsigned int height);
			HostMatrix(const CudaMatrix<T>& c);
			~HostMatrix(void);
		
			HostMatrix<T>& operator=(const CudaMatrix<T>& c);
			
			void copy_submatrix(const CudaMatrix<T>& c, unsigned int elements = 0);
		
			HostMatrix<T>& resize(unsigned int width, unsigned int height = 1);
			HostMatrix<T>& fill(const T& value);
		
		// TODO: implement these
		private:
			HostMatrix(bool pinned); // TODO: habilitar
			HostMatrix(const HostMatrix<T>& c);
			const HostMatrix& operator=(const HostMatrix<T>& c);
			void copy_submatrix(const HostMatrix<T>& c, unsigned int elements = 0);
		
			bool pinned;
			void alloc_data(void);
			void dealloc_data(void);
	};
	
	template <class T> class CudaMatrix : public Matrix<T> {
		public:
			CudaMatrix(void);
			CudaMatrix(const CudaMatrix<T>& c);
			CudaMatrix(const HostMatrix<T>& c);		
			CudaMatrix(unsigned int width, unsigned int height = 1);
			~CudaMatrix(void);
		
			CudaMatrix<T>& resize(unsigned int width, unsigned int height = 1);

			CudaMatrix& operator=(const HostMatrix<T>& c);
			CudaMatrix& operator=(const CudaMatrix<T>& c);
		
			void copy_submatrix(const HostMatrix<T>& c, unsigned int elements = 0);
			void copy_submatrix(const CudaMatrix<T>& c, unsigned int elements = 0);
	};
	
	typedef HostMatrix<double> HostMatrixDouble;
	typedef HostMatrix<float> HostMatrixFloat;
	typedef HostMatrix<float1> HostMatrixFloat1;
	typedef HostMatrix<float2> HostMatrixFloat2;
	typedef HostMatrix<float3> HostMatrixFloat3;
	typedef HostMatrix<float4> HostMatrixFloat4;
	typedef HostMatrix<uint> HostMatrixUInt;
	typedef HostMatrix<uint1> HostMatrixUInt1;	
	
	typedef CudaMatrix<float> CudaMatrixFloat;
	typedef CudaMatrix<float1> CudaMatrixFloat1;
	typedef CudaMatrix<float2> CudaMatrixFloat2;
	typedef CudaMatrix<float3> CudaMatrixFloat3;
	typedef CudaMatrix<float4> CudaMatrixFloat4;	
	typedef CudaMatrix<uint> CudaMatrixUInt;
	typedef CudaMatrix<uint1> CudaMatrixUInt1;
}

#endif
