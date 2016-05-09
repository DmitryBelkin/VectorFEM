#ifndef FILE_COMPORATOR_H
#define FILE_COMPORATOR_H

#include <vector>
#include <math.h>
using namespace std;

#define ZERO 1e-20
#define INF 1e9

class Comparator 
{
public:
	static pair<double, double> abs2(pair<double, double> a)
	{
		return make_pair(abs(a.first), abs(a.second));
	}

	static double scalar(pair<double, double> a, pair<double, double> b)
	{
		return a.first * b.first + a.second * b.second;
	}

	static bool isZero(double a, double eps) 
	{
		if (abs(a) < eps) {
			return true;
		}
		return false;
	}

	static bool isZero(double a) 
	{
		return isZero(a, ZERO);
	}

	static bool isEqual(double a, double b, double eps) 
	{
		return isZero(a - b, eps);
	}

	static bool isEqual(double a, double b) 
	{
		return isEqual(a, b, ZERO);
	}

	static bool nonDecreasingSequence(double a, double b, double c) 
	{
		if (a <= b && b <= c) {
			return true;
		}
		return false;
	}
};

#endif