#include "MatrixVectorOP.h"

Mtrx::Mtrx() {}

Mtrx::Mtrx(uint32_t size)
	: size(size) 
{
	element.resize(size, std::vector<double> (size, 0));
}

Mtrx Mtrx::operator + (const Mtrx& M) const 
{
	Mtrx result = Mtrx(size);

	for (uint32_t i = 0; i < size; ++i) 
		for (uint32_t j = 0; j < size; ++j)
			result.element[i][j] = element[i][j] + M.element[i][j];

	return result;
}

void Mtrx::operator *= (double d) 
{
	for (uint32_t i = 0; i < size; i++)
		for (uint32_t j = 0; j < size; j++)
			element[i][j] = d * element[i][j];
}


Vctr::Vctr() {}

Vctr::Vctr(uint32_t size)
	: size(size) 
{
	element.resize(size, 0);
}

Vctr::Vctr(std::vector<double> init_vector) 
	: size(init_vector.size()) 
{
	element.resize(size);

	for (uint32_t i = 0; i < size; ++i)
		element[i] = init_vector[i];
}

void Vctr::operator *= (double d) 
{
	for (uint32_t i = 0; i < size; ++i)
		element[i] = d * element[i];
}
