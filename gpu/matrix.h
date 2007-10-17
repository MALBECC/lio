#ifndef __G2G_MATRIX_H__
#define __G2G_MATRIX_H__

namespace G2G {
	class Matrix {
		public:
			Matrix(void);
			Matrix(float elem);
			Matrix(unsigned int width, unsigned height = 1);
			~Matrix(void);		
		

			float* data;
			unsigned int width;
			unsigned int height;
			// unsigned int components;
	
			unsigned int bytes(void) const;
			unsigned int elements(void) const;		
	};
}

#endif
