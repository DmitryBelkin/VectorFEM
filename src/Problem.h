#ifndef FILE_PROBLEM_H
#define FILE_PROBLEM_H

#include <vector>
using namespace std;

#include "Comparator.h"
#include "Tester.h"

#ifdef LINEAR_NONFULL
#define NUMBER_OF_LOCAL_BASIS_FUNCTIONS 4
#define NUMBER_OF_BORDER_BASIS_FUNCTIONS 1
#endif

#ifdef LINEAR_FULL
#define NUMBER_OF_LOCAL_BASIS_FUNCTIONS 8
#define NUMBER_OF_BORDER_BASIS_FUNCTIONS 2
#endif

#ifdef QUADRATIC
#define NUMBER_OF_LOCAL_BASIS_FUNCTIONS 12
#define NUMBER_OF_BORDER_BASIS_FUNCTIONS 2
#endif

class Problem
{
public:
	static pair<double, double> rightSide(double x, double y) 
	{
		switch(TEST_ID) {
		case 1:	return make_pair(1 + y, 0);
		case 2:	return make_pair(0, 2 * x);
		case 3: return make_pair(y, x);
		case 4: return make_pair(x, 0);
		case 5: return make_pair(0, y);
		case 6: return make_pair(x, y);
		case 7: return make_pair(x*y + 2, 2*x*y + 1);
		case 8: return make_pair(y*y - 2, x*x - 2); 
		case 9: return make_pair(x*y*y, x*x*y); 
		case 10: return make_pair(x*x*y*y + 4*x*y - 2*x*x, x*x*y*y + 4*x*y - 2*y*y); 
		case 11: return make_pair(y*y*y - 6*y, x*x*x - 6*x);
		case 12: return make_pair(x*x*x, y*y*y);
		case 13: return make_pair(exp(x) * sin(y), exp(x) * cos(y));
		}
		return make_pair(1 + y, 0);
	}

	static pair<double, double> A(double x, double y) 
	{
		switch(TEST_ID) {
		case 1:	return make_pair(1 + y, 0);
		case 2:	return make_pair(0, 2*x);
		case 3: return make_pair(y, x);
		case 4: return make_pair(x, 0);
		case 5: return make_pair(0, y);
		case 6: return make_pair(x, y);
		case 7: return make_pair(x*y, 2*x*y);
		case 8: return make_pair(y*y, x*x); 
		case 9: return make_pair(x*y*y, x*x*y); 
		case 10: return make_pair(x*x*y*y, x*x*y*y); 
		case 11: return make_pair(y*y*y, x*x*x);
		case 12: return make_pair(x*x*x, y*y*y);
		case 13: return make_pair(exp(x) * sin(y), exp(x) * cos(y));
		}
		return make_pair(1 + y, 0);
	}

	static double mu(double x, double y)
	{
		switch(TEST_ID) {
		case 1:	
		case 2: 
		case 3: 
		case 4: 
		case 5: 
		case 6: 
		case 7: 
		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
			return 1;
		}
		return 1;
	}

	static double gamma(double x, double y)
	{
		switch(TEST_ID) {
		case 1:	
		case 2: 
		case 3: 
		case 4:
		case 5:
		case 6: 
		case 7: 
		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
			return 1;
		}
		return 1;
	}

	static double etta(double x, double y)
	{
		switch(TEST_ID) {
		case 1:	return -1;
		case 2:	return 2;
		case 3: return 0;
		case 4: return 0;
		case 5: return 0;
		case 6: return 0;
		case 7: return 2*y - x;
		case 8: return 2*x - 2*y;
		case 9: return 0;
		case 10: return 2*x*y*y - 2*x*x*y;
		case 11: return 3*x*x - 3*y*y;
		case 12: return 0;
		case 13: return 0;
		}
		return -1;
	}

	static double firstBoundary(double x, double y, pair<double, double> tao)
	{
		switch(TEST_ID) {
		case 1:	
		case 2: 
		case 3: 
		case 4:
		case 5:
		case 6: 
		case 7: 
		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
			return Comparator::scalar(A(x, y), tao);
		}
		return Comparator::scalar(A(x, y), tao);
	}
};

#endif