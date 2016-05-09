#include "Tester.h"

void Tester::generate(int doubleMesh)
{
	switch(TEST_ID) {
	case 1: 
	case 2: 
	case 3: 
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
	case 9:
	case 10:
	case 11:
	case 12:
	case 13:
		generateGrid(doubleMesh);
		generateCheckPoints();
	}
}

void Tester::generateGrid(int doubleMesh)
{
	//int nx = 1 * doubleMesh, ny = 1 * doubleMesh;
	int nx = 8 * doubleMesh, ny = 8 * doubleMesh;
	x0 = 1.0; hx = 5.0;
	y0 = 1.0; hy = 5.0;

	fstream out("../resources/grid.txt");
	out << nx << '\t' << ny << endl;
	out << x0 << ' ' << x0 + hx << endl;
	out << y0 << ' ' << y0 + hy << endl;
	out.close();
	out.open("../resources/boundary.txt");
	out << "1 2 1 1\n";
	out.close();
}

void Tester::generateCheckPoints()
{
	int nx = 5;
	int ny = 5;
	n = nx * ny;
	x = new double[n];
	y = new double[n];
	for (int i = 0; i < nx * ny; i++) {
		x[i] = (static_cast<double>(i % nx) * hx / static_cast<double>(nx)) + x0 + 0.2;
		y[i] = (static_cast<double>(i / nx) * hy / static_cast<double>(ny)) + y0 + 0.3;
	}
}
