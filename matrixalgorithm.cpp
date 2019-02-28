//---------------------------------------------------------------------------

#pragma hdrstop

#include "matrixalgorithm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
//------------------------------------------------------------------------------
//
// SAT_VecMat.cpp
//
// Purpose:
//
//   Vector/matrix operations
//
// Notes:
//
//   This software is protected by national and international copyright.
//   Any unauthorized use, reproduction or modificaton is unlawful and
//   will be prosecuted. Commercial and non-private application of the
//   software in any form is strictly prohibited unless otherwise granted
//   by the authors.
//
//   The code is provided without any warranty; without even the implied
//   warranty of merchantibility or fitness for a particular purpose.
//
// Last modified:
//
//   2000/03/04  OMO  Final version (1st edition)
//   2000/11/25  OMO  Initialized imax in LU_Decomp
//   2001/08/17  OMO  Minor upgrades
//
// (c) 1999-2001  O. Montenbruck, E. Gill
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>

#include "GNU_iomanip.h"

using namespace std;


//------------------------------------------------------------------------------
//
// LU_Decomp
//
// Purpose:
//
//   LU-Decomposition.
//
//   Given an nxn matrix A, this routine replaces it by the LU decomposition
//   of a rowwise permutation of itself. A is output, arranged as in
//   equation (2.3.14) of Press et al. (1986); Indx is an ouput vector which
//   records the row permutation effected by partial pivoting. This routine is
//   used in combination with LU_BackSub to solve linear equations or invert
//   a matrix.
//
// Input/output:
//
//   A       Square matrix; replaced by LU decomposition of permutation of A
//           on output
//   Indx    Permutation index vector
//
// Note:
//
//   Adapted from LUDCMP of Press et al. (1986).
//
//------------------------------------------------------------------------------

void LU_Decomp ( Matrix& A, Vector& Indx )
{

  // Constants

  const int    n    = A.size1();
  const double tiny = 1.0e-20;       // A small number

  // Variables

  int     imax=0;
  int     i,j,k;
  double  aAmax, Sum, Dum;
  Vector  V(n);

  // Loop over rows to get scaling information

  for (i=0; i<n; i++) {
	aAmax = 0.0;
	for (j=0;j<n;j++) if (fabs(A(i,j)) > aAmax ) aAmax=fabs(A(i,j));
    if (aAmax==0.0) {
      // No nonzero largest element
      cerr << "ERROR: Singular matrix A in LU_Decomp";
      exit(1);
    };
    V(i) = 1.0/aAmax;           // V stores the implicit scaling of each row
  };

  // Loop over columns of Crout's method

  for ( j=0; j<n; j++ ) {

    if (j > 0) {
	  for ( i=0; i<j; i++ ) {   // This is equation 2.3.12 except for i=j
        Sum = A(i,j);
        if (i>0) {
          for ( k=0; k<i; k++ )  Sum -= A(i,k)*A(k,j);
          A(i,j) = Sum;
        };
      };
    };

    aAmax=0.0;                  // Initialize for the search of the largest
                                // pivot element

    for ( i=j; i<n; i++ ) {     // This is i=j of equation 2.3.12 and
      Sum = A(i,j);             // i=j+1..N of equation 2.3.13
      if (j > 0) {
        for ( k=0; k<j; k++ ) Sum -= A(i,k)*A(k,j);
        A(i,j) = Sum;
      };
      Dum = V(i)*fabs(Sum);     // Figure of merit for the pivot
      if (Dum >= aAmax) {       // Is it better than the best so far ?
        imax  = i;
        aAmax = Dum;
      };
    };

	if (j != imax) {            // Do we need to interchange rows?
      for ( k=0; k<n; k++) {    // Yes, do so ...
		Dum = A(imax,k);
        A(imax,k) = A(j,k);
		A(j,k) = Dum;
      }
      V(imax) = V(j);           // Also interchange the scale factor
    };

	Indx(j) = imax;

    if (j != n-1) {             // Now finally devide by the pivot element
      if (A(j,j) == 0.0) {      // If the pivot element is zero the matrix
        A(j,j) = tiny;          // is singular (at least to the precision of
      };                        // the algorithm). For some applications on
      Dum=1.0/A(j,j);           // singular matrices, it is desirable to
      for (i=j+1;i<n;i++) {     // substitude tiny for zero.
        A(i,j)=A(i,j)*Dum;
      };
    };

  };   // Go back for the next column in the reduction

  if (A(n-1,n-1)==0.0) A(n-1,n-1)=tiny;

};


//------------------------------------------------------------------------------
//
// LU_BackSub
//
// Purpose:
//
//   LU Backsubstitution
//
//   Solves the set of n linear equations Ax=b. Here A is input, not as the
//   matrix A but rather as its LU decomposition, determined by the function
//   LU_Decomp. b is input as the right-hand side vector b, and returns with
//   the solution vector x. A and Indx are not modified by this function and
//   can be left in place for successive calls with different right-hand
//   sides b. This routine takes into account the posssibility that B will
//   begin with many zero elements, so it is efficient for use in matrix
//   inversions.
//
// Input/output:
//
//   A       LU decomposition of permutation of A
//   Indx    Permutation index vector
//   b       Right-hand side vector b; replaced by solution x of Ax=b on output
//
//------------------------------------------------------------------------------

void LU_BackSub ( Matrix& A, Vector& Indx, Vector& b )
{

  // Constants

  const int  n = A.size1();

  // Local variables

  int     ii,i,ll,j;
  double  Sum;

  //
  // Start
  //

  ii = -1;                      // When ii is set to a nonegative value, it will
                                // become the first nonvanishing element of B.
  for (i=0; i<n; i++) {         // We now do the forward substitution.
    ll = (int) Indx(i);         // The only wrinkle is to unscramble the
	Sum = b(ll);                // permutation as we go.
    b(ll) = b(i);
    if (ii != -1) {
      for (j=ii; j<i; j++) Sum -= A(i,j)*b(j);
    }
    else {
      if (Sum != 0.0) ii = i;   // A nonzero element was encountered, so from
    };                          // now on we will have to do the sums in the
    b(i) = Sum;                 // loop above.
   };

   for (i=n-1; i>=0; i--) {     // Now we do the backsubstitution, eqn 2.3.7.
     Sum=b(i);
     if (i<n-1) {
       for (j=i+1;j<n;j++) {
          Sum = Sum-A(i,j)*b(j);
       };
     };
     b(i) = Sum/A(i,i);         // Store a component of the solution vector X.
   };

};





