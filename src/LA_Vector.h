#ifndef CLASS_LA_VECTOR
#define CLASS_LA_VECTOR

#include <vector>
using namespace std;

class LA_Vector 
{
public:
	int size;
	vector<double> element;

	LA_Vector();
	LA_Vector(int);
	LA_Vector(int, double);
	LA_Vector(vector<double>);
	void operator *= (double d);
};

#endif