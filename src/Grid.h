#ifndef FILE_GRID_H
#define FILE_GRID_H

#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "Border.h"
#include "Element.h"

#define EPS 1e-20

class Grid
{
public:
	int number_of_bf;
	vector<Border> border;
	vector<Element> element;

	Grid();
};

class RectangleGridGenerator
{
public:
	void generateFromFiles(Grid& grid, string, string);
	void buildElements(Grid& grid, double*, double*);
	void buildBorders(Grid& grid, string, double*, double*);
	void resizeGrid(Grid& grid);
	void getBfidForElement(unsigned int*, int, int);
	void getBfidForBorder(unsigned int*, int);
	
private:
	int nx, ny;
	double x0, x1, y0, y1;
};

#endif