//---------------------------------------------------------------------------

#ifndef matrixalgorithmH
#define matrixalgorithmH
//---------------------------------------------------------------------------

#include "SATCVector.h"
#include "SATCMatrix.h"

void LU_Decomp ( Matrix& A, Vector& Indx );
void LU_BackSub ( Matrix& A, Vector& Indx, Vector& b );

#endif
