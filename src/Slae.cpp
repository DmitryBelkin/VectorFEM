#include "Slae.h"
#include <iomanip>
#include <iostream>

uint32_t LOC::solveLOC() 
{
	Aq.resize(size);
	r .resize(size);
	z .resize(size);
	p .resize(size);
	dd.resize(size);
	const double normOfRightSide = norm(rightSide);
	double relativeDiscrepancy;
	uint32_t numberOfIterations = 0;

	multiplyMatrixByVector(solution, Aq);
	for (uint32_t i = 0; i < size; ++i)
		r[i] = z[i] = (rightSide[i] - Aq[i]);

	multiplyMatrixByVector(z, p);
	do {
		const double pr = dotProduct(p, r);
		const double pp = dotProduct(p, p);
		const double alpha = pr / pp;

		for (uint32_t i = 0; i < size; ++i) 
		{
			solution[i] += alpha * z[i];
			r[i] -= alpha * p[i];
		}
		multiplyMatrixByVector(r, Aq);

		const double betta = -dotProduct(p, Aq) / pp;
		for (uint32_t i = 0; i < size; ++i) 
		{
			z[i] = r[i] + betta * z[i];
			p[i] = Aq[i] + betta * p[i];
		}
		relativeDiscrepancy = norm(r) / normOfRightSide;
		++numberOfIterations;
	} while (relativeDiscrepancy > precison && numberOfIterations < maxIter);

	multiplyMatrixByVector(solution, r);
	for (uint32_t i = 0; i < size; ++i)
		r[i] -= rightSide[i];

	relativeDiscrepancy = norm(r);
	cout << numberOfIterations << '\t' << setprecision(10) << relativeDiscrepancy << endl;

	return numberOfIterations;
}

void LOC::multiplyMatrixByVector(vector<double> &v, vector<double> &r) 
{
	for (uint32_t i = 0; i < size; ++i)
		r[i] = diagonal[i] * v[i];

	for (uint32_t i = 0; i < size; ++i) 
		for (auto j : matrix[i]) 
		{
			r[i      ] += j.second * v[j.first];
			r[j.first] += j.second * v[i];
		}
}

void LOC::multiplyPreconditioner(vector<double> &r) 
{
	for (uint32_t i = 0; i < size; ++i)
		r[i] /= diagonal[i];
}

void LOC::multiplyPreconditionerByVector(vector<double> &v, vector<double> &r) 
{
	for (uint32_t i = 0; i < size; ++i)
		r[i] = v[i] / diagonal[i];
}

double LOC::norm(vector<double> &v) 
{
	double r = 0;

	for (uint32_t i = 0; i < size; ++i)
		r += v[i] * v[i];

	return sqrt(r);
}

double LOC::dotProduct(vector<double> &u, vector<double> &v) 
{
	double r = 0;

	for (uint32_t i = 0; i < size; ++i)
		r += u[i] * v[i];

	return r;
}
