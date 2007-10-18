#ifndef __G2G_MATRIX_H__
#define __G2G_MATRIX_H__

namespace G2G {
	class Matrix {
		public:
			Matrix(void);
			virtual ~Matrix(void);				

			float* data;
			unsigned int width;
			unsigned int height;
			// unsigned int components;
	
			unsigned int bytes(void) const;
			unsigned int elements(void) const;
	};
	
	class CudaMatrix;
	
	class HostMatrix : public Matrix {
		public:
			HostMatrix(void);
			HostMatrix(float elem);
			HostMatrix(unsigned int width, unsigned height = 1);
			~HostMatrix(void);
			
			void copy_submatrix(const CudaMatrix& c, unsigned int elements = 0);
	};
	
	class CudaMatrix : public Matrix {
		public:
			CudaMatrix(void);
			CudaMatrix(const CudaMatrix& c);
			CudaMatrix(const HostMatrix& c);		
			CudaMatrix(unsigned int width, unsigned height = 1);
			~CudaMatrix(void);		

			CudaMatrix& operator=(const HostMatrix& c);
			CudaMatrix& operator=(const CudaMatrix& c);
		
			void copy_submatrix(const HostMatrix& c, unsigned int elements = 0);
	};
}

#endif
