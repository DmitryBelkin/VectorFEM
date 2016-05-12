#ifndef TASK_H
#define TASK_H

#include "Grid.h"
#include "Slae.h"
#include "DiscreteSLAE.h"
#include <cstdint>

class Task
{
public:
	Grid* m_grid;
	LOC* m_slae;

	Task() {  m_grid = new Grid(); };
	virtual ~Task() {};

	void solveProblem();
	pair<double, double> getResult(const double, const double) const;
	void output(string filename, double* x, double* y, uint32_t n);
};

#endif TASK_H
