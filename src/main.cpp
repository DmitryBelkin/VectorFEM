#include "Task.h"
#include "DomainModel.h"

int main()
{
	DomainModel * domain = new DomainModel();
	Task        * task   = new Task();

	domain ->inputSplitNumbers();
	domain ->createGrid       ();
	task   ->solveProblem     ();
	domain ->createCheckPoints();
	task   ->output("../resources/result.txt", domain->m_x, domain->m_y, domain->m_n);

	return 0;
}
