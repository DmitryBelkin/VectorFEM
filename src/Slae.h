#ifndef FILE_SLAE_H
#define FILE_SLAE_H

#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include "Grid.h"

using namespace std;

#define PRECISION 1e-16
#define MAX_ITERATIONS 10000

class Slae
{
public:
	unsigned int size;
	vector<map<int, double>> matrix;
	vector<double> diagonal;
	vector<double> rightSide;
	vector<double> solution;

	int solve();
	int solveUsingLocallyOptimalScheme();
	void multiplyMatrixByVector (vector<double> &, vector<double> &);

private:
	vector<double> Aq;
	vector<double> r;
	vector<double> z;
	vector<double> p;
	vector<double> dd;

	void multiplyPreconditioner (vector<double> &);
	void multiplyPreconditionerByVector (vector<double> &, vector<double> &);
	double norm (vector<double> &);
	double dotProduct (vector<double> &, vector<double> &);
};

class SlaeBuilder 
{
protected:
	Slae* slae;

public:
	SlaeBuilder();
	virtual ~SlaeBuilder();
	virtual void createSlae();
	virtual void fillUsingGrid(Grid*);
	virtual Slae* getSlae();
};

class SlaeDirector 
{
public:
	Slae* build(SlaeBuilder&, Grid*);
};

#endif