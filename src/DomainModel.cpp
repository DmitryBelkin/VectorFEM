#include "DomainModel.h"
#include <fstream>

void DomainModel::inputSplitNumbers()
{
	std::fstream in("../resources/split.txt");
	in >> m_nx >> m_ny;
	in.close();
}

void DomainModel::createGrid()
{
	m_x0 = 1.0; m_hx = 9.0;
	m_y0 = 1.0; m_hy = 9.0;

	std::fstream out("../resources/grid.txt");
	out << m_nx << '\t' << m_ny << std::endl;
	out << m_x0 << ' ' << m_x0 + m_hx << std::endl;
	out << m_y0 << ' ' << m_y0 + m_hy << std::endl;
	out.close();
	out.open("../resources/boundary.txt");
	out << "1 2 1 1\n";
	out.close();
}

void DomainModel::createCheckPoints()
{
	const uint32_t nx = 4;
	const uint32_t ny = 5;
	m_n = nx * ny;
	m_x = new double[m_n];
	m_y = new double[m_n];
	for (uint32_t i = 0; i < nx * ny; ++i) 
	{
		m_x[i] = (static_cast<double>(i % nx) * m_hx / static_cast<double>(nx)) + m_x0 + 0.2;
		m_y[i] = (static_cast<double>(i / nx) * m_hy / static_cast<double>(ny)) + m_y0 + 0.3;
	}
}
