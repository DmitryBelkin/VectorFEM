#ifndef MATRIX_VECTOR_OP_H
#define MATRIX_VECTOR_OP_H

#include <vector>
#include <cstdint>

class Mtrx 
{
public:
	uint32_t size;
	std::vector<std::vector<double>> element;

	Mtrx();
	Mtrx(uint32_t);
	
	Mtrx operator + (const Mtrx&) const;
	void operator *= (double);
};

class Vctr 
{
public:
	uint32_t size;
	std::vector<double> element;

	Vctr();
	Vctr(uint32_t);
	Vctr(std::vector<double>);
	void operator *= (double d);
};

#endif MATRIX_VECTOR_OP_H
