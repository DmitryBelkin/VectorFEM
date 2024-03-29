#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <array>
#include <functional>

using namespace std;

class Integrate 
{
public:
	static double gauss5_1d(double x0, double x1, function<double (double) > f) 
	{
		array<double, 5> weight     = {0.5688888888888889,  0.4786286704993665, 0.4786286704993665,  0.2369268850561891, 0.2369268850561891};
		array<double, 5> coordinate = {0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640};
		const double scale_x = (x1 - x0) / 2;
		const double shift_x = (x1 + x0) / 2;
		double result = 0;

		for (int i = 0; i < 5; ++i)
		{
			const double x = coordinate[i] * scale_x + shift_x;
			const double ff = f(x);
			result += weight[i] * ff * scale_x;
		}

		return result;
	}

	static double gauss5_2d(double x0, double x1, double y0, double y1, function<double (double, double)> f) 
	{
		array<double, 5> weight     = {0.5688888888888889,  0.4786286704993665, 0.4786286704993665,  0.2369268850561891, 0.2369268850561891};
		array<double, 5> coordinate = {0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640};
		const double scale_x = (x1 - x0) / 2.0;
		const double shift_x = (x1 + x0) / 2.0;
		const double scale_y = (y1 - y0) / 2.0;
		const double shift_y = (y1 + y0) / 2.0;
		double result = 0;

		for (int i = 0; i < 5; ++i)
		{
			const double x = coordinate[i] * scale_x + shift_x;
			for (int j = 0; j < 5; ++j)
			{
				const double y = coordinate[j] * scale_y + shift_y;
				result += weight[i] * weight[j] * f(x, y) * scale_x * scale_y;
			}
		}

		return result;
	}

	static double gauss5_2d(double x0, double x1, function<double (double)> y0, function<double (double)> y1, function<double (double, double)> f) 
	{
		array<double, 5> weight     = {0.5688888888888889,  0.4786286704993665, 0.4786286704993665,  0.2369268850561891, 0.2369268850561891};
		array<double, 5> coordinate = {0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640};
		double scale_x = (x1 - x0) / 2;
		double shift_x = (x1 + x0) / 2;
		double scale_y;
		double shift_y;
		double result = 0;

		for (int i = 0; i < 5; i++) 
		{
			const double x = coordinate[i] * scale_x + shift_x;
			const double ymin = min(y0(x), y1(x));
			const double ymax = max(y0(x), y1(x));
			scale_y = (ymax - ymin) / 2;
			shift_y = (ymax + ymin) / 2;

			for (int j = 0; j < 5; ++j) 
			{
				const double y = coordinate[j] * scale_y + shift_y;
				const double ff = f(x, y);
				result += weight[i] * weight[j] * ff * scale_x * scale_y;
			}
		}

		return result;
	}
};

#endif INTEGRATE_H
