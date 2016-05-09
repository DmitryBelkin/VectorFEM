#include "Border.h"

Border::Border() {}

Border::Border(
	int type,
	double sign,
	double x0, 
	double x1, 
	double y0, 
	double y1, 
	unsigned int b[4]
) : type(type), sign(sign), x0(x0), x1(x1), y0(y0), y1(y1)
{
	for (int i = 0; i < NUMBER_OF_BORDER_BASIS_FUNCTIONS; i++)
		bf[i] = b[i];
}

pair<double, double> Border::tao() const {
	double hx = x1 - x0;
	double hy = y1 - y0;
	double h = sqrt(hx * hx + hy * hy);

	return make_pair(hx / h, hy / h);
}


double Border::psitao(int i, double t) const 
{
	switch(i) {
	case 0: return sign * 1.0;
	case 1: return (2.0 * t - 1.0);
	}
	return 0;
}
