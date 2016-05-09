#ifndef CLASS_LA_MATRIX
#define CLASS_LA_MATRIX

#include <vector>
using namespace std;

class LA_Matrix 
{
public:
	unsigned int size;
	vector<vector<double>> element;

	LA_Matrix();
	LA_Matrix(int);
	
	LA_Matrix operator + (const LA_Matrix&) const;
	void operator *= (double);
};

#endif