#ifndef FILE_BORDER_H
#define FILE_BORDER_H

#include "Comparator.h"
#include "LA_Matrix.h"
#include "LA_Vector.h"
#include "Problem.h"
#include "Quadrature.h"

class Border 
{
public:
	int type;
	double sign;
	double x0, x1, y0, y1;
	unsigned int bf[NUMBER_OF_BORDER_BASIS_FUNCTIONS];

	Border();
	Border(int, double, double, double, double, double, unsigned int[NUMBER_OF_BORDER_BASIS_FUNCTIONS]);
	pair<double, double> tao() const;
	double psitao(int, double) const;
};


class BorderMethod 
{
public:
	static void calculateDirichletVector(const Border& base, LA_Vector& vect) 
	{
		LA_Matrix matrix(NUMBER_OF_BORDER_BASIS_FUNCTIONS);
		LA_Vector right_side(NUMBER_OF_BORDER_BASIS_FUNCTIONS);

		for (int i = 0; i < NUMBER_OF_BORDER_BASIS_FUNCTIONS; i++) {
			right_side.element[i] += calculateRightSideElement(i, base);
			for (int j = 0; j < NUMBER_OF_BORDER_BASIS_FUNCTIONS; j++) {
				matrix.element[i][j] = calculateMassMatrixElement(i, j, base);
			}
		}
		solveLU(matrix, right_side, vect);
	}

	static void calculateNeumannVector(const Border& base, LA_Vector& vect) 
	{
		for (int i = 0; i < NUMBER_OF_BORDER_BASIS_FUNCTIONS; i++) {
			vect.element[i] += calculateBoundaryElement(i, base);
		}
	}

private:
	static double calculateMassMatrixElement(int i, int j, const Border& base) 
	{
		double hx = base.x1 - base.x0;
		double hy = base.y1 - base.y0;
		double h = sqrt(hx * hx + hy * hy);
		double res = Quadrature::gauss5_1d(
			0, 
			1.0, 
			[&] (double t) {
				double x = (1.0 - t) * base.x0 + t * base.x1;
				double y = (1.0 - t) * base.y0 + t * base.y1;

				return h * base.psitao(i, t) * base.psitao(j, t);}
		);

		return res;
	}

	static double calculateRightSideElement(int i, const Border& base) 
	{
		double hx = base.x1 - base.x0;
		double hy = base.y1 - base.y0;
		double h = sqrt(hx * hx + hy * hy);
		double res = Quadrature::gauss5_1d(
			0, 
			1.0, 
			[&] (double t) {
				double x = (1.0 - t) * base.x0 + t * base.x1;
				double y = (1.0 - t) * base.y0 + t * base.y1;

				return h * base.psitao(i, t) * Problem::firstBoundary(x, y, base.tao());}
		);

		return res;
	}

	static double calculateBoundaryElement(int i, const Border& base) 
	{
		double hx = base.x1 - base.x0;
		double hy = base.y1 - base.y0;
		double h = sqrt(hx * hx + hy * hy);
		double res = Quadrature::gauss5_1d(
			0, 
			1.0, 
			[&] (double t) {
				double x = (1.0 - t) * base.x0 + t * base.x1;
				double y = (1.0 - t) * base.y0 + t * base.y1;

				return h * Problem::etta(x, y) * base.psitao(i, t) / Problem::mu(x, y);}
		);

		return res;
	}

	static void solveLU(LA_Matrix A, LA_Vector b, LA_Vector& x) 
	{
		int N = NUMBER_OF_BORDER_BASIS_FUNCTIONS;
		vector<vector<double>> LU(N, vector<double>(N));
		vector<double> y(N);

		for (int i = 0; i < N; i++) {		
			for (int j = i; j < N; j++) {
				LU[i][j] = A.element[i][j];
				for (int k = 0; k < i; k++) {
					LU[i][j] -= LU[i][k] * LU[k][j];
				}
			}
			for (int j = i + 1; j < N; j++) {
				LU[j][i] = A.element[j][i];
				for (int k = 0; k < i; k++) {
					LU[j][i] -= LU[j][k] * LU[k][i];
				}
				LU[j][i] /= LU[i][i];
			}
		}
		for (int i = 0; i < N; i++) {
			y[i] = b.element[i];
			for (int k = 0; k < i; k++) {
				y[i] -= y[k] * LU[i][k];
			}
		}
		for (int i = 0; i < N; i++) {
			x.element[N - i - 1] = y[N - i - 1];
			for (int k = N - i; k < N; k++) {
				x.element[N - i - 1] -= x.element[k] * LU[N - i - 1][k];
			}
			x.element[N - i - 1] /= LU[N - i - 1][N - i - 1];
		}
	}
};


#endif