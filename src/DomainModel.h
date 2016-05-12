#ifndef DOMAIN_MODEL_H
#define DOMAIN_MODEL_H

#define NUM 12

#include <cstdint>

class DomainModel
{
public:
	uint32_t   m_n;
	double   * m_x;
	double   * m_y;
	uint32_t   m_nx, m_ny;

	void inputSplitNumbers();
	void createGrid();
	void createCheckPoints();

private:
	double m_x0, m_hx, m_y0, m_hy;
};

#endif DOMAIN_MODEL_H
