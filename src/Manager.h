#ifndef FILE_MANAGER_H
#define FILE_MANAGER_H

#include "Grid.h"
#include "Slae.h"
#include "SpraseSlaeBuilder.h"

class Manager
{
public:
	Grid* grid;
	Slae* slae;

	Manager();
	void doIt();
	pair<double, double> getResult(double, double);
};


#endif