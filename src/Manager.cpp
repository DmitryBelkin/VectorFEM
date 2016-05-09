#include "Manager.h"

Manager::Manager() {}

void Manager::doIt()
{
	grid = new Grid;
	RectangleGridGenerator generator;
	generator.generateFromFiles(*grid, "../resources/grid.txt", "../resources/boundary.txt");

	SlaeDirector director;
	SpraseSlaeBuilder builder;
	slae = director.build(builder, grid);
	slae->solve();
}

pair<double, double> Manager::getResult(double x, double y)
{
	pair<double, double> res = make_pair(0, 0);

	for (int i = 0; i < grid->element.size(); i++) {
		if (grid->element[i].x0 <= x && x <= grid->element[i].x1 &&
			grid->element[i].y0 <= y && y <= grid->element[i].y1) 
		{
			for (int j = 0; j < NUMBER_OF_LOCAL_BASIS_FUNCTIONS; j++) {
				Element base = grid->element[i];
				pair<double, double> psi = base.psi(j, x, y);

				res.first += slae->solution[grid->element[i].bf[j]] * psi.first;
				res.second += slae->solution[grid->element[i].bf[j]] * psi.second;
			}
			break;
		}
	}

	return res;
}