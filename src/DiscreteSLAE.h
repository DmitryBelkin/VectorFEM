#ifndef DISCRETE_SLAE_H
#define DISCRETE_SLAE_H

#include "FE.h"
#include "MatrixVectorOP.h"
#include "Slae.h"
#include <cstdint>

using namespace std;

class DiscreteSLAE
{
public:
	DiscreteSLAE();
	~DiscreteSLAE();

	LOC* m_slae;

	void generateDomainElements(Grid*);
	LOC* getSlae();

private:
	void addElementToGlobal          (Mtrx*         , Vctr*       , FE*);
	void addElementToGlobal          (Vctr*         , Border*          );
	void addSolutionToGlobal         (Vctr*         , Border*          );
	void clearColumn                 (const uint32_t, const double     );
	void clearElement                (const uint32_t, const double     );
	void clearRow                    (const uint32_t, const double     );
	void fillSubSlaeOfBorder         (Border*                          );
	void fillSubSlaeOfDirichletBorder(Border*                          );
	void fillSubSlaeOfElement        (FE*                              );
	void fillSubSlaeOfNeumannBorder  (Border*                          );
	Mtrx massMatrix                  (FE*                              );
	Mtrx stiffnessMatrix             (FE*                              );
	Vctr rightSideVector             (FE*                              );
	void setSizeOfSlae               (uint32_t                         );
	void buildSlaePortraitFor        (Grid*                            );
	void allocateMemoryForDiagonalAndRightSide(   );
	void allocateMemoryForSubmatrixOfElement  (FE*);
};

#endif DISCRETE_SLAE_H
