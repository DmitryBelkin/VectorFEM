#include <fstream>
#include <iomanip>
using namespace std;

#include "Manager.h"
#include "Tester.h"

Manager man;
Tester tester;

void printResult(string filename, double* x, double* y, int n)
{
	ofstream out(filename);

	for (int i = 0; i < n; i++) {
		pair<double, double> res = man.getResult(x[i], y[i]);
		out.width(6);
		out.precision(4);
		out << x[i];
		out.width(6);
		out.precision(4);
		out << y[i];
		out.width(18);
		out.precision(10);
		out << res.first;
		out.width(18);
		out.precision(10);
		out << res.second;
		out << endl;
	}

	out.close();
}

int main()
{
	tester.generate(1);
	man.doIt();
	
	printResult("../resources/result.txt", tester.x, tester.y, tester.n);

	return 0;
}