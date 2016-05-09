#ifndef FILE_SPRASESLAEBUILDER_H
#define FILE_SPRASESLAEBUILDER_H

#include <functional>
#include <map>
#include <vector>
#include "Element.h"
#include "LA_Matrix.h"
#include "LA_Vector.h"
#include "Slae.h"

using namespace std;

class SpraseSlaeBuilder : public SlaeBuilder 
{
public:
	void createSlae();
	void fillUsingGrid(Grid*);

private:
	void addElementToGlobal(LA_Vector*, Border*);
	void addElementToGlobal(LA_Matrix*, LA_Vector*, Element*);
	void addSolutionToGlobal(LA_Vector*, Border*);
	void allocateMemoryForDiagonalAndRightSide();
	void allocateMemoryForSubmatrixOfElement(Element*);
	void buildSlaePortraitFor(Grid*);
	void clearColumn(int, double);
	void clearElement(int, double);
	void clearRow(int, double);
	void fillSubSlaeOfBorder(Border*);
	void fillSubSlaeOfDirichletBorder(Border*);
	void fillSubSlaeOfElement(Element*);
	void fillSubSlaeOfNeumannBorder(Border*);
	void setSizeOfSlae(unsigned int);
	LA_Matrix massMatrix(Element*);
	LA_Matrix stiffnessMatrix(Element*);
	LA_Vector rightSideVector(Element*);
};

#endif