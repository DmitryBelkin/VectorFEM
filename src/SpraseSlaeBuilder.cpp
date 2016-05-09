#include "SpraseSlaeBuilder.h"

void SpraseSlaeBuilder::createSlae() 
{
	slae = new Slae;
}

void SpraseSlaeBuilder::fillUsingGrid(Grid* grid) 
{
	buildSlaePortraitFor(grid);
	for (unsigned int k = 0; k < grid->element.size(); k++) {
		fillSubSlaeOfElement(&grid->element[k]);
	}
	for (unsigned int k = 0; k < grid->border.size(); k++) {
		fillSubSlaeOfBorder(&grid->border[k]);
	}
}

void SpraseSlaeBuilder::buildSlaePortraitFor(Grid* grid) 
{
	setSizeOfSlae(grid->number_of_bf);
	allocateMemoryForDiagonalAndRightSide();
	for (unsigned int k = 0; k < grid->element.size(); k++) {
		allocateMemoryForSubmatrixOfElement(&grid->element[k]);
	}
}

void SpraseSlaeBuilder::setSizeOfSlae(unsigned int size) 
{
	slae->size = size;
}

void SpraseSlaeBuilder::allocateMemoryForDiagonalAndRightSide() 
{
	slae->matrix.resize(slae->size);
	slae->diagonal.resize(slae->size);
	slae->rightSide.resize(slae->size);
	slae->solution.resize(slae->size);
}

void SpraseSlaeBuilder::allocateMemoryForSubmatrixOfElement(Element* current) 
{
	unsigned int local_size = NUMBER_OF_LOCAL_BASIS_FUNCTIONS;

	for (unsigned int i = 0; i < local_size * local_size; i++) {
		unsigned int row = current->bf[i / local_size];
		unsigned int col = current->bf[i % local_size];

		if (row > col) {
			slae->matrix[row][col] = 0;
		}
	}
}

void SpraseSlaeBuilder::fillSubSlaeOfElement(Element* current) 
{
	LA_Matrix M = massMatrix(current);
	LA_Matrix G = stiffnessMatrix(current);
	LA_Matrix local_matrix = M + G;
	LA_Vector local_right_side = rightSideVector(current);

	addElementToGlobal(&local_matrix, &local_right_side, current);
}

void SpraseSlaeBuilder::fillSubSlaeOfBorder(Border* current) 
{
	if (current->type == 1)
		fillSubSlaeOfDirichletBorder(current);

	else if (current->type == 2)
		fillSubSlaeOfNeumannBorder(current);
}

void SpraseSlaeBuilder::fillSubSlaeOfDirichletBorder(Border* current) 
{
	LA_Vector local_solution(NUMBER_OF_BORDER_BASIS_FUNCTIONS);
	BorderMethod::calculateDirichletVector(*current, local_solution);
	addSolutionToGlobal(&local_solution, current);
}

void SpraseSlaeBuilder::fillSubSlaeOfNeumannBorder(Border* current) 
{
	LA_Vector local_right_side(NUMBER_OF_BORDER_BASIS_FUNCTIONS);
	BorderMethod::calculateNeumannVector(*current, local_right_side);
	addElementToGlobal(&local_right_side, current);
}

void SpraseSlaeBuilder::addElementToGlobal(LA_Matrix* A, LA_Vector* b, Element* current) 
{
	for (unsigned int i = 0; i < A->size; i++) {
		int row = current->bf[i];
		slae->rightSide[row] += b->element[i];

		for (unsigned int j = 0; j < A->size; j++) {
			int col = current->bf[j];
			if (row > col) 
				slae->matrix[row][col] += A->element[i][j];
			else if (row == col) 
				slae->diagonal[row] += A->element[i][j];
		}
	}
}

void SpraseSlaeBuilder::addElementToGlobal(LA_Vector* b, Border* current) 
{
	for (unsigned int i = 0; i < b->size; i++) {
		int row = current->bf[i];
		slae->rightSide[row] += b->element[i];
	}
}

void SpraseSlaeBuilder::addSolutionToGlobal(LA_Vector* b, Border* current) 
{
	for (int i = 0; i < NUMBER_OF_BORDER_BASIS_FUNCTIONS; i++) {
		clearElement(current->bf[i], b->element[i]);
	}
}

void SpraseSlaeBuilder::clearElement(int id, double value) 
{
	clearRow(id, value);
	clearColumn(id, value);

	slae->diagonal[id] = 1;
	slae->rightSide[id] = value;
}

void SpraseSlaeBuilder::clearColumn(int id, double value) 
{
	for (unsigned int i = 0; i < slae->matrix.size(); i++) {
		for (auto &j : slae->matrix[i]) {
			if (j.first == id) {
				slae->rightSide[i] -= j.second * value;
				j.second = 0;
			}
		}
	}
}

void SpraseSlaeBuilder::clearRow(int row, double value) 
{
	for (auto &j : slae->matrix[row]) {
		slae->rightSide[j.first] -= j.second * value;
		j.second = 0;
	}
}

LA_Matrix SpraseSlaeBuilder::massMatrix(Element* current)
{
	LA_Matrix matrix(NUMBER_OF_LOCAL_BASIS_FUNCTIONS);
	ElementMethod::calculateMassMatrix(*current, matrix);
	return matrix;
}

LA_Matrix SpraseSlaeBuilder::stiffnessMatrix(Element* current) 
{
	LA_Matrix matrix(NUMBER_OF_LOCAL_BASIS_FUNCTIONS);
	ElementMethod::calculateStiffnessMatrix(*current, matrix);
	return matrix;
}

LA_Vector SpraseSlaeBuilder::rightSideVector(Element* current) 
{
	LA_Vector vect(NUMBER_OF_LOCAL_BASIS_FUNCTIONS);
	ElementMethod::calculateRightSideVector(*current, vect);
	return vect;
}


