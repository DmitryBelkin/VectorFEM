#include "DiscreteSLAE.h"
#include <cassert>

DiscreteSLAE::DiscreteSLAE() { m_slae = new LOC(); };

LOC* DiscreteSLAE::getSlae() { return m_slae; };

void DiscreteSLAE::generateDomainElements(Grid* grid) 
{
	buildSlaePortraitFor(grid);
	for (uint32_t k = 0; k < grid->element.size(); ++k) 
		fillSubSlaeOfElement(&grid->element[k]);
	for (uint32_t k = 0; k < grid->border.size(); ++k) 
		fillSubSlaeOfBorder(&grid->border[k]);
}

void DiscreteSLAE::buildSlaePortraitFor(Grid* grid) 
{
	setSizeOfSlae(grid->number_of_bf);
	allocateMemoryForDiagonalAndRightSide();
	for (uint32_t k = 0; k < grid->element.size(); ++k)
		allocateMemoryForSubmatrixOfElement(&grid->element[k]);
}

void DiscreteSLAE::setSizeOfSlae(uint32_t size) { m_slae->size = size; }

void DiscreteSLAE::allocateMemoryForDiagonalAndRightSide() 
{
	m_slae->matrix   .resize(m_slae->size);
	m_slae->diagonal .resize(m_slae->size);
	m_slae->rightSide.resize(m_slae->size);
	m_slae->solution .resize(m_slae->size);
}

void DiscreteSLAE::allocateMemoryForSubmatrixOfElement(FE* current) 
{
	const uint32_t local_size = amountLocalBf;

	for (uint32_t i = 0; i < local_size * local_size; ++i)
	{
		const uint32_t row = current->m_bf[i / local_size];
		const uint32_t col = current->m_bf[i % local_size];

		if (row > col)
			m_slae->matrix[row][col] = 0;
	}
}

void DiscreteSLAE::fillSubSlaeOfElement(FE* current) 
{
	Mtrx M = massMatrix     (current);
	Mtrx G = stiffnessMatrix(current);

	Mtrx local_matrix     = M + G;
	Vctr local_right_side = rightSideVector(current);

	addElementToGlobal(&local_matrix, &local_right_side, current);
}

void DiscreteSLAE::fillSubSlaeOfBorder(Border* current) 
{
	switch(current->m_type)
	{
	case 1 : fillSubSlaeOfDirichletBorder(current); break;
	case 2 : fillSubSlaeOfNeumannBorder  (current); break;
	default: assert(!"Error in 'fillSubSlaeOfBorder'.");
	}
}

void DiscreteSLAE::fillSubSlaeOfDirichletBorder(Border* current) 
{
	Vctr local_solution(amountBorderBf);
	BorderMethod::calculateDirichletVector(*current, local_solution);
	addSolutionToGlobal(&local_solution, current);
}

void DiscreteSLAE::fillSubSlaeOfNeumannBorder(Border* current) 
{
	Vctr local_right_side(amountBorderBf);
	BorderMethod::calculateNeumannVector(*current, local_right_side);
	addElementToGlobal(&local_right_side, current);
}

void DiscreteSLAE::addElementToGlobal(Mtrx* A, Vctr* b, FE* current) 
{
	for (uint32_t i = 0; i < A->size; ++i)
	{
		const uint32_t row = current->m_bf[i];
		m_slae->rightSide[row] += b->element[i];

		for (uint32_t j = 0; j < A->size; j++) 
		{
			const uint32_t col = current->m_bf[j];
			if (row > col) 
				m_slae->matrix[row][col] += A->element[i][j];
			else 
			if (row == col) 
				m_slae->diagonal[row] += A->element[i][j];
		}
	}
}

void DiscreteSLAE::addElementToGlobal(Vctr* b, Border* current) 
{
	for (uint32_t i = 0; i < b->size; ++i) 
	{
		const uint32_t row = current->m_bf[i];
		m_slae->rightSide[row] += b->element[i];
	}
}

void DiscreteSLAE::addSolutionToGlobal(Vctr* b, Border* current) 
{
	for (uint32_t i = 0; i < amountBorderBf; ++i)
		clearElement(current->m_bf[i], b->element[i]);
}

void DiscreteSLAE::clearElement(const uint32_t id, const double value) 
{
	clearRow   (id, value);
	clearColumn(id, value);

	m_slae->diagonal [id] = 1;
	m_slae->rightSide[id] = value;
}

void DiscreteSLAE::clearColumn(const uint32_t id, const double value) 
{
	for (uint32_t i = 0; i < m_slae->matrix.size(); ++i) 
	{
		for (auto &j : m_slae->matrix[i]) 
		{
			if (j.first == id) 
			{
				m_slae->rightSide[i] -= j.second * value;
				j.second = 0;
			}
		}
	}
}

void DiscreteSLAE::clearRow(const uint32_t row, const double value) 
{
	for (auto &j : m_slae->matrix[row]) 
	{
		m_slae->rightSide[j.first] -= j.second * value;
		j.second = 0;
	}
}

Mtrx DiscreteSLAE::massMatrix(FE* current)
{
	Mtrx matrix(amountLocalBf);
	ElementMethod::calculateMassMatrix(*current, matrix);
	return matrix;
}

Mtrx DiscreteSLAE::stiffnessMatrix(FE* current) 
{
	Mtrx matrix(amountLocalBf);
	ElementMethod::calculateStiffnessMatrix(*current, matrix);
	return matrix;
}

Vctr DiscreteSLAE::rightSideVector(FE* current) 
{
	Vctr vect(amountLocalBf);
	ElementMethod::calculateRightSideVector(*current, vect);
	return vect;
}
