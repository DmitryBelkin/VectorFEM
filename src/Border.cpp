#include "Border.h"

Border::Border() {}

Border::Border(int type,double sign,double x0, double x1, double y0, double y1, uint32_t b[4]) :
	m_type(type),
	m_sign(sign),
	m_x0  (x0)  ,
	m_x1  (x1)  ,
	m_y0  (y0)  ,
	m_y1  (y1)
{
	for (uint32_t i = 0; i < amountBorderBf; ++i) 
		m_bf[i] = b[i];
}

pair<double, double> Border::tao() const 
{
	const double hx = m_x1 - m_x0;
	const double hy = m_y1 - m_y0;
	const double h  = sqrt(hx * hx + hy * hy);

	return make_pair(hx / h, hy / h);
}

double Border::psitao(uint32_t i, double t) const 
{
	switch(i) 
	{
	case 0 : return m_sign * 1.0   ;
	case 1 : return (2.0 * t - 1.0);
	default: return 0              ;
	}
}
