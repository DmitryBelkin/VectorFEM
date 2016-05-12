#ifndef FE_H
#define FE_H

#include "MatrixVectorOP.h"
#include "Problem.h"
#include "Integrate.h"
#include <cstdint>

class FE 
{
public:
	double   m_y0, m_x1, m_y1, m_x0;
	uint32_t m_bf[amountLocalBf];

	FE();
	FE(double, double, double, double, uint32_t[amountLocalBf]);

	double rotpsi(uint32_t, double, double        ) const;
	double dphi  (uint32_t, double, double, double) const;
	double phi   (uint32_t, double, double, double) const;
	double phim  (double  , double, double        ) const;
	pair<double, double> psi(uint32_t, double, double) const;
};

class ElementMethod 
{
public:
	static void calculateMassMatrix(const FE& base, Mtrx& matrix) 
	{
		calculateMatrix(base, matrix, [&] (uint32_t i, uint32_t j, const FE& e) { return calculateMassElement(i, j, e);} );
	}

	static void calculateStiffnessMatrix(const FE& base, Mtrx& matrix) 
	{
		calculateMatrix(base, matrix, [&] (uint32_t i, uint32_t j, const FE& e) { return calculateStiffnessElement(i, j, e);} );
	}

	static void calculateRightSideVector(const FE& base, Vctr& vect) 
	{
		for (uint32_t i = 0; i < amountLocalBf; ++i) { vect.element[i] = calculateRightSideElement(i, base); }
	}

private:
	static void calculateMatrix(const FE& base, Mtrx& matrix, function<double (int, int, const FE&)> f) 
	{
		for (uint32_t i = 0; i < amountLocalBf; ++i) 
		{
			for (uint32_t j = 0; j < amountLocalBf; ++j) 
			{
				matrix.element[i][j] = f(i, j, base);
			}
		}
	}

	static double calculateElement(const FE& base, function<double (double, double)> f) 
	{
		return Integrate::gauss5_2d(base.m_x0, base.m_x1, base.m_y0, base.m_y1, [&] (double x, double y) { return f(x, y); } );
	}

	static double calculateMassElement(uint32_t i, uint32_t j, const FE& base) 
	{
		return calculateElement(base, [&] (double x, double y) 
		{
			const double koef = base.psi(i, x, y).first * base.psi(j, x, y).first + base.psi(i, x, y).second * base.psi(j, x, y).second;
			return Problem::gamma(x, y) * koef;
		});
	}

	static double calculateStiffnessElement(uint32_t i, uint32_t j, const FE& base) 
	{
		return calculateElement(base, [&] (double x, double y) { return base.rotpsi(i, x, y) * base.rotpsi(j, x, y) / Problem::mu(x, y); } );
	}

	static double calculateRightSideElement(uint32_t i, const FE& base) 
	{
		return Integrate::gauss5_2d(base.m_x0, base.m_x1, base.m_y0, base.m_y1, [&] (double x, double y) 
		{
			const double koeff = Problem::rightSide(x, y).first * base.psi(i, x, y).first + Problem::rightSide(x, y).second * base.psi(i, x, y).second;
			return Problem::gamma(x, y) * koeff;
		});
	}
};

#endif FE_H
