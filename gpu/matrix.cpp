#include <iostream>
#include "matrix.h"
using namespace G2G;

Matrix::Matrix(void) : data(NULL), width(0), height(0) /*, components(0)*/ {}

Matrix::Matrix(float elem) : width(1), height(1) {
	data = new float[1];
	data[0] = elem;
}

Matrix::Matrix(unsigned int _width, unsigned _height) : width(_width), height(_height) {
	data = new float[width * height];	
}

Matrix::~Matrix(void) {
	delete[] data;	
}

unsigned int Matrix::bytes(void) const {
	return elements() * sizeof(float);	
}

unsigned int Matrix::elements(void) const {
	return width * height /* * components */;
}
