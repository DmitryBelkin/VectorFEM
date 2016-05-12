#ifndef SLAE_H
#define SLAE_H

#include "Grid.h"
#include <map>
#include <cstdint>

using namespace std;

const double   precison =  1e-16;
const uint32_t maxIter  =  10000;

class LOC
{
public:
	uint32_t size;
	vector<map<uint32_t, double>> matrix;
	vector<double> diagonal;
	vector<double> rightSide;
	vector<double> solution;

	uint32_t solveLOC();
	void multiplyMatrixByVector (vector<double> &, vector<double> &);

private:
	vector<double> Aq;
	vector<double> r ;
	vector<double> z ;
	vector<double> p ;
	vector<double> dd;

	void   multiplyPreconditioner (vector<double> &);
	void   multiplyPreconditionerByVector (vector<double> &, vector<double> &);
	double norm (vector<double> &);
	double dotProduct (vector<double> &, vector<double> &);
};

#endif SLAE_H
