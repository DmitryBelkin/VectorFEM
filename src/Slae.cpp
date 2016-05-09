#include "Slae.h"

int Slae::solve() 
{
	return solveUsingLocallyOptimalScheme();
}

SlaeBuilder::SlaeBuilder() {}
SlaeBuilder::~SlaeBuilder() {}
void SlaeBuilder::createSlae() {}
void SlaeBuilder::fillUsingGrid(Grid*) {}
Slae* SlaeBuilder::getSlae() 
{
	return slae;
}

Slae* SlaeDirector::build(SlaeBuilder& builder, Grid *grid) 
{
	builder.createSlae();
	builder.fillUsingGrid(grid);
	return builder.getSlae();
}

int Slae::solveUsingLocallyOptimalScheme() 
{
	Aq.resize(size);
	r.resize(size);
	z.resize(size);
	p.resize(size);
	dd.resize(size);
	double normOfRightSide = norm(rightSide);
	double relativeDiscrepancy;
	int numberOfIterations = 0;

	multiplyMatrixByVector(solution, Aq);
	for (unsigned int i = 0; i < size; i++) {
		r[i] = z[i] = (rightSide[i] - Aq[i]);
	}
	multiplyMatrixByVector(z, p);
	do {
		double pr = dotProduct(p, r);
		double pp = dotProduct(p, p);
		double alpha = pr / pp;

		for (unsigned int i = 0; i < size; i++) {
			solution[i] += alpha * z[i];
			r[i] -= alpha * p[i];
		}
		multiplyMatrixByVector(r, Aq);

		double betta = -dotProduct(p, Aq) / pp;
		for (unsigned int i = 0; i < size; i++) {
			z[i] = r[i] + betta * z[i];
			p[i] = Aq[i] + betta * p[i];
		}
		relativeDiscrepancy = norm(r) / normOfRightSide;
		numberOfIterations++;
	} while (relativeDiscrepancy > PRECISION && numberOfIterations < MAX_ITERATIONS);

	multiplyMatrixByVector(solution, r);
	for (unsigned int i = 0; i < size; i++) {
		r[i] -= rightSide[i];
	}
	relativeDiscrepancy = norm(r);
	cout << numberOfIterations << '\t' << setprecision(10) << relativeDiscrepancy << endl;

	return numberOfIterations;
}

void Slae::multiplyMatrixByVector (vector<double> &v, vector<double> &r) 
{
	for (unsigned int i = 0; i < size; i++) {
		r[i] = diagonal[i] * v[i];
	}
	for (unsigned int i = 0; i < size; i++) {
		for (auto j : matrix[i]) {
			r[i] += j.second * v[j.first];
			r[j.first] += j.second * v[i];
		}
	}
}

void Slae::multiplyPreconditioner (vector<double> &r) 
{
	for (unsigned int i = 0; i < size; i++) {
		r[i] /= diagonal[i];
	}
}

void Slae::multiplyPreconditionerByVector (vector<double> &v, vector<double> &r) 
{
	for (unsigned int i = 0; i < size; i++) {
		r[i] = v[i] / diagonal[i];
	}
}

double Slae::norm (vector<double> &v) 
{
	double r = 0;

	for (unsigned int i = 0; i < size; i++) {
		r += v[i] * v[i];
	}
	return sqrt(r);
}

double Slae::dotProduct (vector<double> &u, vector<double> &v) 
{
	double r = 0;

	for (unsigned int i = 0; i < size; i++) {
		r += u[i] * v[i];
	}
	return r;
}
