#pragma hdrstop
#pragma argsused

#include <tchar.h>
#include <stdio.h>
#include <iomanip>

/*
  Testprogramm für die Matrix/Vektor-Bilbiothek aus

  O.Montenbruck, E. Gill: Satellite Orbits - Models, Methods, and Applications. Springer Verlag, Heidelberg, (2000).

  (c) 02.05.2012 Fredie Kern
*/

#include "SATCVector.h"
#include "SATCMatrix.h"

int main(int argc,char* argv[])
{
	Vector a(1,2,3);
	Vector b(-2,1,0);
	Vector c;

	Matrix A,B,C;

	c = Cross(a,b);
	c = a+b;
	std::cout << std::setw(8) << c << "\n";

	A = Diag(c);
	std::cout << std::setw(8) << A << "\n";

	a = A * c;
	std::cout << std::setw(8) << a << "\n";

	A = Dyadic(a,b);
	std::cout << std::setw(8) << A << "\n";

	B = Inv(A);
	std::cout << std::setw(8) << B << "\n";

	C = B * A;
	std::cout << std::setw(8) << C << "\n";

	return 0;
}
