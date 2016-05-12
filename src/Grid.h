#ifndef GRID_H
#define GRID_H

#include "Border.h"
#include "FE.h"
#include <string>
#include <vector>

class Grid
{
public:
	Grid() {};

	int number_of_bf;
	std::vector<Border > border;
	std::vector<FE> element;
};

class RectangleGridGenerator
{
public:
	void generateFromFiles(Grid& grid, const std::string, const std::string);
	void buildElements    (Grid& grid, double*          , double*);
	void buildBorders     (Grid& grid, const std::string, double*, double*);
	void resizeGrid       (Grid& grid);
	void getBfidForElement(unsigned int*, const int, const int) const;
	void getBfidForBorder (unsigned int*, const int           ) const;
	
private:
	uint32_t m_nx, m_ny;
	double   m_x0, m_x1, m_y0, m_y1;
};

#endif GRID_H
