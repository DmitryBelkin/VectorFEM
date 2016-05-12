#include "FE.h"

FE::FE() {}

FE::FE(double x0, double x1, double y0, double y1, uint32_t b[4]) :
	m_x0(x0),
	m_x1(x1),
	m_y0(y0),
	m_y1(y1)
{
	for (uint32_t i = 0; i < amountLocalBf; ++i)
		m_bf[i] = b[i];
}

double FE::rotpsi(uint32_t n, double x, double y) const
{
	switch (n) 
	{
	case 0 : return -dphi( 0   , m_y0, m_y1, y );
	case 1 : return -dphi( 1   , m_y0, m_y1, y );
	case 2 : return  dphi( 0   , m_x0, m_x1, x );
	case 3 : return  dphi( 1   , m_x0, m_x1, x );
	case 4 : return -phim( m_x0, m_x1, x       ) * dphi( 0, m_y0, m_y1, y );
	case 5 : return -phim( m_x0, m_x1, x       ) * dphi( 1, m_y0, m_y1, y );
	case 6 : return  phim( m_y0, m_y1, y       ) * dphi( 0, m_x0, m_x1, x );
	case 7 : return  phim( m_y0, m_y1, y       ) * dphi( 1, m_x0, m_x1, x );
	case 8 : return -dphi( 2   , m_y0, m_y1, y );
	case 9 : return  dphi( 2   , m_x0, m_x1, x );
	case 10: return -phim( m_x0, m_x1, x       ) * dphi( 2, m_y0, m_y1, y );
	case 11: return  phim( m_y0, m_y1, y       ) * dphi( 2, m_x0, m_x1, x );
	}
	return 0;
}

pair<double, double> FE::psi(uint32_t n, double x, double y) const
{
	if (x < m_x0 || x > m_x1 || y < m_y0 || y > m_y1) 
		return make_pair(0, 0);

	switch (n)
	{
	case 0 : return make_pair(phi(0, m_y0, m_y1, y)                      , 0                                          );
	case 1 : return make_pair(phi(1, m_y0, m_y1, y)                      , 0                                          );
	case 2 : return make_pair(0                                          , phi(0, m_x0, m_x1, x)                      );
	case 3 : return make_pair(0                                          , phi(1, m_x0, m_x1, x)                      );
	case 4 : return make_pair(phim(m_x0, m_x1, x) * phi(0, m_y0, m_y1, y), 0                                          );
	case 5 : return make_pair(phim(m_x0, m_x1, x) * phi(1, m_y0, m_y1, y), 0                                          );
	case 6 : return make_pair(0                                          , phim(m_y0, m_y1, y) * phi(0, m_x0, m_x1, x));
	case 7 : return make_pair(0                                          , phim(m_y0, m_y1, y) * phi(1, m_x0, m_x1, x));
	case 8 : return make_pair(phi(2, m_y0, m_y1, y)                      , 0                                          );
	case 9 : return make_pair(0                                          , phi(2, m_x0, m_x1, x)                      );
	case 10: return make_pair(phim(m_x0, m_x1, x) * phi(2, m_y0, m_y1, y), 0                                          );
	case 11: return make_pair(0                                          , phim(m_y0, m_y1, y) * phi(2, m_x0, m_x1, x));
	default: return make_pair(0                                          , 0                                          );
	}
}

double FE::dphi(uint32_t n, double t0, double t1, double t) const
{
	if (t < t0 || t > t1) 
		return 0;

	switch (n) 
	{
	case 0 : return -1.                     / (t1 - t0)            ;
	case 1 : return  1.                     / (t1 - t0)            ;
	case 2 : return  4. * (t1 + t0 - 2 * t) / (t1 - t0) * (t1 - t0);
	default: return  0                                             ;
	}
}

double FE::phim(double t0, double t1, double t) const
{
	if (t < t0 || t > t1) 
		return 0;

	return (2. * t - (t0 + t1)) / (t1 - t0);
}

double FE::phi(uint32_t n, double t0, double t1, double t) const
{
	if (t < t0 || t > t1) 
		return 0;

	switch (n) 
	{
	case 0 : return (t1 - t)                / (t1 - t0)            ;
	case 1 : return (t - t0)                / (t1 - t0)            ;
	case 2 : return 4 * (t - t0) * (t1 - t) / (t1 - t0) * (t1 - t0);
	default: return 0;
	}
}
