#include "Grid.h"
#include <fstream>
#include <cstdint>

void RectangleGridGenerator::generateFromFiles(Grid& grid, const string file_grid, const string file_boundary)
{
	ifstream in(file_grid);
	in >> m_nx >> m_ny;
	in >> m_x0 >> m_x1 >> m_y0 >> m_y1;
	in.close();

	double *x = new double[m_nx + 1];
	double *y = new double[m_ny + 1];

	for (uint32_t i = 0; i < m_nx + 1; ++i) x[i] = m_x0 * (m_nx - i) / m_nx + m_x1 * i / m_nx;
	for (uint32_t i = 0; i < m_ny + 1; ++i) y[i] = m_y0 * (m_ny - i) / m_ny + m_y1 * i / m_ny;

	resizeGrid(grid);
	buildElements(grid, x, y);
	buildBorders(grid, file_boundary, x, y);
}

void RectangleGridGenerator::buildElements(Grid& grid, double* x, double* y)
{
	for (uint32_t j = 0; j < m_ny; ++j) 
	{
		for (uint32_t i = 0; i < m_nx; ++i) 
		{
			uint32_t bf_id[amountLocalBf];
			getBfidForElement(bf_id, i, j);
			grid.element[j * m_nx + i] = FE(x[i], x[i + 1], y[j], y[j + 1], bf_id);
		}
	}
}

void RectangleGridGenerator::buildBorders(Grid& grid, const string file_boundary, double* x, double* y)
{
	int bot, top, left, right;
	ifstream in(file_boundary);
	in >> bot >> top >> left >> right;
	in.close();

	for (uint32_t i = 0; i < m_nx; ++i) 
	{
		uint32_t bf[amountBorderBf];
		getBfidForBorder(bf, i);
		grid.border.push_back(Border(bot, 1.0, x[i], x[i + 1], m_y0, m_y0, bf));
	}
	for (uint32_t i = 0; i < m_nx; ++i) 
	{
		uint32_t bf[amountBorderBf];
		getBfidForBorder(bf, m_ny * m_nx + i);
		grid.border.push_back(Border(top, -1.0, x[i + 1], x[i], m_y1, m_y1, bf));
	}
	for (uint32_t j = 0; j < m_ny; ++j) 
	{
		uint32_t bf[amountBorderBf];
		getBfidForBorder(bf, (m_ny + 1) * m_nx + j * (m_nx + 1));
		grid.border.push_back(Border(left, -1.0, m_x0, m_x0, y[j + 1], y[j], bf));
	}
	for (uint32_t j = 0; j < m_ny; ++j) 
	{
		uint32_t bf[amountBorderBf];
		getBfidForBorder(bf, (m_ny + 1) * m_nx + j * (m_nx + 1) + m_nx);
		grid.border.push_back(Border(right, 1.0, m_x1, m_x1, y[j], y[j + 1], bf));
	}
}

void RectangleGridGenerator::resizeGrid(Grid& grid)
{
	const uint32_t nln = (m_nx + 1) * m_ny + (m_ny + 1) * m_nx;
	const uint32_t nlf = (m_nx + 1) * m_ny + (m_ny + 1) * m_nx;
	const uint32_t nq  = 4 * m_nx * m_ny;

	grid.element.resize(m_nx * m_ny);
	grid.number_of_bf = nln + nlf + nq;
}

void RectangleGridGenerator::getBfidForElement(unsigned int* bf_id, const int i, const int j) const
{
	const uint32_t nln = (m_nx + 1) * m_ny + (m_ny + 1) * m_nx;
	const uint32_t nlf = (m_nx + 1) * m_ny + (m_ny + 1) * m_nx;
	const uint32_t b[amountLocalBf] = 
	{
		j * m_nx + i                                    ,
		(j + 1) * m_nx + i                              ,
		(m_ny + 1) * m_nx + j * (m_nx + 1) + i          ,
		(m_ny + 1) * m_nx + j * (m_nx + 1) + i + 1      ,
		nln + j * m_nx + i                              ,
		nln + (j + 1) * m_nx + i                        ,
		nln + (m_ny + 1) * m_nx + j * (m_nx + 1) + i    ,
		nln + (m_ny + 1) * m_nx + j * (m_nx + 1) + i + 1,
		nln + nlf + 4 * (j * m_nx + i)                  ,
		nln + nlf + 4 * (j * m_nx + i) + 1              ,
		nln + nlf + 4 * (j * m_nx + i) + 2              ,
		nln + nlf + 4 * (j * m_nx + i) + 3
	};

	for (uint32_t j = 0; j < amountLocalBf; ++j)
		bf_id[j] = b[j];
}

void RectangleGridGenerator::getBfidForBorder(unsigned int* bf_id, const int i) const
{
	const uint32_t nln = (m_nx + 1) * m_ny + (m_ny + 1) * m_nx;
	const uint32_t b[amountBorderBf] = 
	{
		i,
		nln + i
	};

	for (uint32_t j = 0; j < amountBorderBf; ++j) 
		bf_id[j] = b[j];
}
