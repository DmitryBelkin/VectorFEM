#include "LA_Vector.h"

LA_Vector::LA_Vector() {}

LA_Vector::LA_Vector(int size) : size(size) 
{
	element.resize(size, 0);
}

LA_Vector::LA_Vector(int size, double init_value) : size(size) 
{
	element.resize(size, init_value);
}

LA_Vector::LA_Vector(vector<double> init_vector) : size(init_vector.size()) 
{
	element.resize(size);

	for (int i = 0; i < size; i++) {
		element[i] = init_vector[i];
	}
}

void LA_Vector::operator *= (double d) 
{
	for (int i = 0; i < size; i++) {
		element[i] = d * element[i];
	}
}