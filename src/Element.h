#ifndef FILE_ELEMENT_H
#define FILE_ELEMENT_H

#include "LA_Matrix.h"
#include "LA_Vector.h"
#include "Problem.h"
#include "Quadrature.h"

class Element 
{
public:
	double y0, x1, y1, x0;
	unsigned int bf[NUMBER_OF_LOCAL_BASIS_FUNCTIONS];

	Element();
	Element(double, double, double, double, unsigned int[NUMBER_OF_LOCAL_BASIS_FUNCTIONS]);
	double rotpsi(int, double, double) const;
	pair<double, double> psi(int, double, double) const;
	double dphi(int, double, double, double) const;
	double phi(int, double, double, double) const;
	double phim(double, double, double) const;
};


class ElementMethod 
{
public:
	static void calculateMassMatrix(const Element& base, LA_Matrix& matrix) 
	{
		calculateMatrix(
			base, 
			matrix, 
			[&] (int i, int j, const Element& e) { 
				return calculateMassElement(i, j, e);}
		);
	}

	static void calculateStiffnessMatrix(const Element& base, LA_Matrix& matrix) 
	{
		calculateMatrix(
			base, 
			matrix, 
			[&] (int i, int j, const Element& e) { 
				return calculateStiffnessElement(i, j, e);}
		);
	}

	static void calculateRightSideVector(const Element& base, LA_Vector& vect) 
	{
		for (int i = 0; i < NUMBER_OF_LOCAL_BASIS_FUNCTIONS; i++) {
			vect.element[i] = calculateRightSideElement(i, base);
		}
	}

private:
	static void calculateMatrix(
		const Element& base, 
		LA_Matrix& matrix, 
		function<double (int, int, const Element&)> f) 
	{
		for (int i = 0; i < NUMBER_OF_LOCAL_BASIS_FUNCTIONS; i++) {
			for (int j = 0; j < NUMBER_OF_LOCAL_BASIS_FUNCTIONS; j++) {
				matrix.element[i][j] = f(i, j, base);
			}
		}
	}

	static double calculateElement(
		const Element& base, 
		function<double (double, double)> f) 
	{
		return Quadrature::gauss5_2d(
			base.x0, 
			base.x1, 
			base.y0,
			base.y1,
			[&] (double x, double y) {return f(x, y);});
	}

	static double calculateMassElement(int i, int j, const Element& base) 
	{
		return calculateElement(
			base, 
			[&] (double x, double y) {
				return Problem::gamma(x, y) * Comparator::scalar(base.psi(i, x, y), base.psi(j, x, y));}
		);
	}

	static double calculateStiffnessElement(int i, int j, const Element& base) {
		return calculateElement(
			base, 
			[&] (double x, double y) {
				return base.rotpsi(i, x, y) * base.rotpsi(j, x, y) / Problem::mu(x, y);}
		);
	}

	static double calculateRightSideElement(int i, const Element& base) 
	{
		return Quadrature::gauss5_2d(
			base.x0, 
			base.x1, 
			base.y0,
			base.y1,
			[&] (double x, double y) {
				return Problem::gamma(x, y) * Comparator::scalar(Problem::rightSide(x, y), base.psi(i, x, y));}
		);
	}
};


#endif