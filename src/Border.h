#ifndef BORDER_H
#define BORDER_H

#include "MatrixVectorOP.h"
#include "Problem.h"
#include "Integrate.h"
#include <cstdint>

class Border 
{
public:
	uint32_t m_type;
	double   m_sign;
	double   m_x0, m_x1, m_y0, m_y1;
	uint32_t m_bf[amountBorderBf];

	Border();
	Border(int, double, double, double, double, double, uint32_t[amountBorderBf]);
	pair<double, double> tao() const;
	double psitao(uint32_t, double) const;
};

class BorderMethod 
{
public:
	static void calculateDirichletVector(const Border& base, Vctr& vect) 
	{
		Mtrx matrix(amountBorderBf);
		Vctr right_side(amountBorderBf);

		for (uint32_t i = 0; i < amountBorderBf; ++i) 
		{
			right_side.element[i] += calculateRightSideElement(i, base);
			for (uint32_t j = 0; j < amountBorderBf; ++j) 
				matrix.element[i][j] = calculateMassMatrixElement(i, j, base);
		}
		solveLU(matrix, right_side, vect);
	}

	static void calculateNeumannVector(const Border& base, Vctr& vect) 
	{
		for (uint32_t i = 0; i < amountBorderBf; ++i)
			vect.element[i] += calculateBoundaryElement(i, base);
	}

private:
	static double calculateMassMatrixElement(int i, int j, const Border& base) 
	{
		const double hx = base.m_x1 - base.m_x0;
		const double hy = base.m_y1 - base.m_y0;
		const double h = sqrt(hx * hx + hy * hy);
		const double res = Integrate::gauss5_1d(0, 1.0, [&] (double t) 
		{
			return h * base.psitao(i, t) * base.psitao(j, t);
		});

		return res;
	}

	static double calculateRightSideElement(int i, const Border& base) 
	{
		const double hx = base.m_x1 - base.m_x0;
		const double hy = base.m_y1 - base.m_y0;
		const double h = sqrt(hx * hx + hy * hy);
		const double res = Integrate::gauss5_1d(0, 1.0, [&] (double t)
		{
			const double x = (1.0 - t) * base.m_x0 + t * base.m_x1;
			const double y = (1.0 - t) * base.m_y0 + t * base.m_y1;
			return h * base.psitao(i, t) * Problem::firstBoundary(x, y, base.tao());
		});

		return res;
	}

	static double calculateBoundaryElement(uint32_t i, const Border& base) 
	{
		const double hx = base.m_x1 - base.m_x0;
		const double hy = base.m_y1 - base.m_y0;
		const double h = sqrt(hx * hx + hy * hy);
		const double res = Integrate::gauss5_1d(0, 1.0, [&] (double t) 
		{
			const double x = (1.0 - t) * base.m_x0 + t * base.m_x1;
			const double y = (1.0 - t) * base.m_y0 + t * base.m_y1;
			return h * Problem::etta(x, y) * base.psitao(i, t) / Problem::mu(x, y);
		});

		return res;
	}

	static void solveLU(Mtrx A, Vctr b, Vctr& x) 
	{
		const uint32_t N = amountBorderBf;
		vector<vector<double>> LU(N, vector<double>(N));
		vector<double> y(N);

		for (uint32_t i = 0; i < N; ++i) 
		{
			for (uint32_t j = i; j < N; ++j) 
			{
				LU[i][j] = A.element[i][j];
				for (uint32_t k = 0; k < i; ++k) 
					LU[i][j] -= LU[i][k] * LU[k][j];
			}
			for (uint32_t j = i + 1; j < N; ++j) 
			{
				LU[j][i] = A.element[j][i];
				for (uint32_t k = 0; k < i; ++k) 
					LU[j][i] -= LU[j][k] * LU[k][i];
				LU[j][i] /= LU[i][i];
			}
		}
		for (uint32_t i = 0; i < N; ++i) 
		{
			y[i] = b.element[i];
			for (uint32_t k = 0; k < i; ++k) 
				y[i] -= y[k] * LU[i][k];
		}
		for (uint32_t i = 0; i < N; ++i) 
		{
			x.element[N - i - 1] = y[N - i - 1];
			for (uint32_t k = N - i; k < N; ++k)
				x.element[N - i - 1] -= x.element[k] * LU[N - i - 1][k];
			x.element[N - i - 1] /= LU[N - i - 1][N - i - 1];
		}
	}
};

#endif BORDER_H
