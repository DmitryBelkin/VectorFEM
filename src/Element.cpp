#include "Element.h"

Element::Element() {}

Element::Element(
	double x0, 
	double x1, 
	double y0, 
	double y1, 
	unsigned int b[4]
) : x0(x0), x1(x1), y0(y0), y1(y1)
{
	for (int i = 0; i < NUMBER_OF_LOCAL_BASIS_FUNCTIONS; i++)
		bf[i] = b[i];
}

double Element::rotpsi(int n, double x, double y) const
{
	switch (n) {
	case 0:	return -dphi(0, y0, y1, y);
	case 1:	return -dphi(1, y0, y1, y);
	case 2:	return  dphi(0, x0, x1, x);
	case 3:	return  dphi(1, x0, x1, x);
	case 4:	return -phim(x0, x1, x) * dphi(0, y0, y1, y);
	case 5:	return -phim(x0, x1, x) * dphi(1, y0, y1, y);
	case 6:	return  phim(y0, y1, y) * dphi(0, x0, x1, x);
	case 7:	return  phim(y0, y1, y) * dphi(1, x0, x1, x);
	case 8: return -dphi(2, y0, y1, y);
	case 9:	return  dphi(2, x0, x1, x);
	case 10: return -phim(x0, x1, x) * dphi(2, y0, y1, y);
	case 11: return  phim(y0, y1, y) * dphi(2, x0, x1, x);
	}
	return 0;
}

pair<double, double> Element::psi(int n, double x, double y) const
{
	if (x < x0 || x > x1 || y < y0 || y > y1) 
		return make_pair(0, 0);

	switch (n) {
	case 0:	return make_pair(phi(0, y0, y1, y), 0);
	case 1:	return make_pair(phi(1, y0, y1, y), 0);
	case 2:	return make_pair(0, phi(0, x0, x1, x));
	case 3:	return make_pair(0, phi(1, x0, x1, x));
	case 4:	return make_pair(phim(x0, x1, x) * phi(0, y0, y1, y), 0);
	case 5:	return make_pair(phim(x0, x1, x) * phi(1, y0, y1, y), 0);
	case 6:	return make_pair(0, phim(y0, y1, y) * phi(0, x0, x1, x));
	case 7:	return make_pair(0, phim(y0, y1, y) * phi(1, x0, x1, x));
	case 8: return make_pair(phi(2, y0, y1, y), 0);
	case 9:	return make_pair(0, phi(2, x0, x1, x));
	case 10: return make_pair(phim(x0, x1, x) * phi(2, y0, y1, y), 0);
	case 11: return make_pair(0, phim(y0, y1, y) * phi(2, x0, x1, x));
	}
	return make_pair(0, 0);
}

double Element::dphi(int n, double t0, double t1, double t) const
{
	if (t < t0 || t > t1) 
		return 0;

	switch (n) {
	case 0: return -1. / (t1 - t0);
	case 1:	return 1. / (t1 - t0);
	case 2: return 4. * (t1 + t0 - 2 * t) / (t1 - t0) * (t1 - t0);
	}
	return 0;
}

double Element::phim(double t0, double t1, double t) const
{
	if (t < t0 || t > t1) 
		return 0;

	return (2. * t - (t0 + t1)) / (t1 - t0);
}

double Element::phi(int n, double t0, double t1, double t) const
{
	if (t < t0 || t > t1) 
		return 0;

	switch (n) {
	case 0: return (t1 - t) / (t1 - t0);
	case 1:	return (t - t0) / (t1 - t0);
	case 2: return 4 * (t - t0) * (t1 - t) / (t1 - t0) * (t1 - t0);
	}
	return 0;
}
