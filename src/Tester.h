#ifndef FILE_TESTER_H
#define FILE_TESTER_H

#include <fstream>
using namespace std;

/* LINEAR_NONFULL / LINEAR_FULL / QUADRATIC */
#define QUADRATIC
#define TEST_ID 13

class Tester
{
public:
	int n;
	double *x;
	double *y;

	void generate(int);

private:
	double x0, hx;
	double y0, hy;

	void generateGrid(int);
	void generateCheckPoints();
};

#endif