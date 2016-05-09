#include "LA_Matrix.h"

LA_Matrix::LA_Matrix() {}

LA_Matrix::LA_Matrix(int size) : size(size) 
{
	element.resize(size, vector<double> (size, 0));
}

LA_Matrix LA_Matrix::operator + (const LA_Matrix& M) const 
{
	LA_Matrix result = LA_Matrix(size);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			result.element[i][j] = element[i][j] + M.element[i][j];
		}
	}

	return result;
}

void LA_Matrix::operator *= (double d) 
{
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			element[i][j] = d * element[i][j];
		}
	}
}