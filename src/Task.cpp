#include "Task.h"
#include <fstream>

void Task::solveProblem()
{
	RectangleGridGenerator * generator = new RectangleGridGenerator();
	DiscreteSLAE           * builder   = new DiscreteSLAE          ();

	generator->generateFromFiles(*m_grid, "../resources/grid.txt", "../resources/boundary.txt");
	builder  ->generateDomainElements(m_grid);
	m_slae   = builder->getSlae();
	m_slae   ->solveLOC();
}

pair<double, double> Task::getResult(const double x, const double y) const
{
	pair<double, double> res = make_pair(0, 0);

	for (uint32_t i = 0; i < m_grid->element.size(); ++i) 
	{
		if (   m_grid->element[i].m_x0 <= x && x <= m_grid->element[i].m_x1 
			&& m_grid->element[i].m_y0 <= y && y <= m_grid->element[i].m_y1) 
		{
			for (uint32_t j = 0; j < amountLocalBf; ++j) 
			{
				FE base = m_grid->element[i];
				pair<double, double> psi = base.psi(j, x, y);

				res.first  += m_slae->solution[m_grid->element[i].m_bf[j]] * psi.first;
				res.second += m_slae->solution[m_grid->element[i].m_bf[j]] * psi.second;
			}
			break;
		}
	}

	return res;
}

void Task::output(string filename, double* x, double* y, uint32_t n)
{
	ofstream out(filename);

	for (uint32_t i = 0; i < n; ++i) 
	{
		pair<double, double> res = getResult(x[i], y[i]);

		out.precision(4);
		out << x[i] << '\t';

		out.precision(4);
		out << y[i] << '\t';

		out.precision(10);
		out << res.first << '\t';

		out.precision(10);
		out << res.second;

		out << endl;
	}

	out.close();
}
