#include "Grid.h"

Grid::Grid() {}

void RectangleGridGenerator::generateFromFiles(
	Grid& grid, 
	string file_grid,
	string file_boundary)
{
	ifstream in(file_grid);
	in >> nx >> ny;
	in >> x0 >> x1 >> y0 >> y1;
	double *x = new double[nx + 1];
	double *y = new double[ny + 1];

	for (int i = 0; i < nx + 1; i++) 
		x[i] = x0 * (nx - i) / nx + x1 * i / nx;
	for (int i = 0; i < ny + 1; i++) 
		y[i] = y0 * (ny - i) / ny + y1 * i / ny;

	resizeGrid(grid);
	buildElements(grid, x, y);
	buildBorders(grid, file_boundary, x, y);
}

void RectangleGridGenerator::buildElements(
	Grid& grid, 
	double* x, 
	double* y)
{
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			unsigned int bf_id[NUMBER_OF_LOCAL_BASIS_FUNCTIONS];
			getBfidForElement(bf_id, i, j);
			grid.element[j * nx + i] = Element(x[i], x[i + 1], y[j], y[j + 1], bf_id);
		}
	}
}

void RectangleGridGenerator::buildBorders(
	Grid& grid, 
	string file_boundary,
	double* x,
	double* y)
{
	int nln = (nx + 1) * ny + (ny + 1) * nx;
	int bot, top, left, right;
	ifstream in(file_boundary);
	in >> bot >> top >> left >> right;

	for (int i = 0; i < nx; i++) {
		unsigned int bf[NUMBER_OF_BORDER_BASIS_FUNCTIONS];
		getBfidForBorder(bf, i);
		grid.border.push_back(Border(bot, 1.0, x[i], x[i + 1], y0, y0, bf));
	}
	for (int i = 0; i < nx; i++) {
		unsigned int bf[NUMBER_OF_BORDER_BASIS_FUNCTIONS];
		getBfidForBorder(bf, ny * nx + i);
		grid.border.push_back(Border(top, -1.0, x[i + 1], x[i], y1, y1, bf));
	}
	for (int j = 0; j < ny; j++) {
		unsigned int bf[NUMBER_OF_BORDER_BASIS_FUNCTIONS];
		getBfidForBorder(bf, (ny + 1) * nx + j * (nx + 1));
		grid.border.push_back(Border(left, -1.0, x0, x0, y[j + 1], y[j], bf));
	}
	for (int j = 0; j < ny; j++) {
		unsigned int bf[NUMBER_OF_BORDER_BASIS_FUNCTIONS];
		getBfidForBorder(bf, (ny + 1) * nx + j * (nx + 1) + nx);
		grid.border.push_back(Border(right, 1.0, x1, x1, y[j], y[j + 1], bf));
	}
}

void RectangleGridGenerator::resizeGrid(Grid& grid)
{
	int nln = (nx + 1) * ny + (ny + 1) * nx;
	int nlf = (nx + 1) * ny + (ny + 1) * nx;
	int nq = 4 * nx * ny;
#ifdef LINEAR_NONFULL
	grid.number_of_bf = nln;
#endif
#ifdef LINEAR_FULL
	grid.number_of_bf = nln + nlf;
#endif
#ifdef QUADRATIC
	grid.number_of_bf = nln + nlf + nq;
#endif
	grid.element.resize(nx * ny);
}

void RectangleGridGenerator::getBfidForElement(unsigned int* bf_id, int i, int j)
{
	int nln = (nx + 1) * ny + (ny + 1) * nx;
	int nlf = (nx + 1) * ny + (ny + 1) * nx;
	unsigned int b[NUMBER_OF_LOCAL_BASIS_FUNCTIONS] = {
#ifdef LINEAR_NONFULL
		j * nx + i,		
		(j + 1) * nx + i,		
		(ny + 1) * nx + j * (nx + 1) + i,		
		(ny + 1) * nx + j * (nx + 1) + i + 1
#endif
#ifdef LINEAR_FULL
		j * nx + i,		
		(j + 1) * nx + i,		
		(ny + 1) * nx + j * (nx + 1) + i,		
		(ny + 1) * nx + j * (nx + 1) + i + 1, 
		nln + j * nx + i,		
		nln + (j + 1) * nx + i,		
		nln + (ny + 1) * nx + j * (nx + 1) + i,		
		nln + (ny + 1) * nx + j * (nx + 1) + i + 1
#endif
#ifdef QUADRATIC
		j * nx + i,		
		(j + 1) * nx + i,		
		(ny + 1) * nx + j * (nx + 1) + i,		
		(ny + 1) * nx + j * (nx + 1) + i + 1, 
		nln + j * nx + i,		
		nln + (j + 1) * nx + i,		
		nln + (ny + 1) * nx + j * (nx + 1) + i,		
		nln + (ny + 1) * nx + j * (nx + 1) + i + 1,
		nln + nlf + 4 * (j * nx + i),
		nln + nlf + 4 * (j * nx + i) + 1,
		nln + nlf + 4 * (j * nx + i) + 2,
		nln + nlf + 4 * (j * nx + i) + 3
#endif
	};

	for (int j = 0; j < NUMBER_OF_LOCAL_BASIS_FUNCTIONS; j++) {
		bf_id[j] = b[j];
	}
}

void RectangleGridGenerator::getBfidForBorder(unsigned int* bf_id, int i)
{
	int nln = (nx + 1) * ny + (ny + 1) * nx;
	unsigned int b[NUMBER_OF_BORDER_BASIS_FUNCTIONS] = {
#ifdef LINEAR_NONFULL
		i
#endif
#ifdef LINEAR_FULL
		i,
		nln + i
#endif
#ifdef QUADRATIC
		i,
		nln + i
#endif
	};
	
	for (int j = 0; j < NUMBER_OF_BORDER_BASIS_FUNCTIONS; j++) {
		bf_id[j] = b[j];
	}
}