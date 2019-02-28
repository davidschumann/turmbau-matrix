/*
   ursprünglich

	 SAT_VecMat.cpp

   aus: O.Montenbruck, E. Gill: Satellite Orbits - Models, Methods, and Applications. Springer Verlag, Heidelberg, (2000).

   02.05.2012
	  Abhängigkeiten zwischen class Matrix und class Vector entflochten, so dass jede Klasse für sich
	  ohne die andere genutzt werden kann.
*/
//---------------------------------------------------------------------------

#pragma hdrstop

#include "SATCMatrix.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

#include <iomanip>


#include "matrixalgorithm.h"
#include "GNU_iomanip.h"

using namespace std;


//------------------------------------------------------------------------------
//
// Matrix (class implementation)
//
// Purpose:
//
//   Matrix data type and associated operations
//
//------------------------------------------------------------------------------


// Constructors

Matrix::Matrix ()                          // Matrix without elements
  : n(0), m(0)
{
  M = 0;
}

Matrix::Matrix (int dim1, int dim2)        // Nullmatrix of specified shape
  : n(dim1), m(dim2)
{
  int i,j;
  // Memory allocation
  M = new double*[dim1];
  for (i=0; i<dim1; i++) M[i] = new double[dim2];
  // Initialization
  for (i=0; i<dim1; i++) {
	for (j=0; j<dim2; j++) M[i][j]=0.0;
  }
}

Matrix::Matrix (const Matrix& M_)          // Copy
{
  int i,j;
  n = M_.n;
  m = M_.m;
  // Memory allocation
  M = new double*[n];
  for (i=0; i<n; i++) M[i] = new double[m];
  // Initialization
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) M[i][j]=M_.M[i][j];
  }
}

Matrix::Matrix (const double* p, int dim1, int dim2)   // Array copy
{
  int i,j;
  n = dim1;
  m = dim2;
  // Memory allocation
  M = new double*[n];
  for (i=0; i<n; i++) M[i] = new double[m];
  // Initialization
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) M[i][j]=p[i*dim2+j];
  }

}

Matrix::~Matrix()
{
  for (int i=0; i<n; i++) delete[] M[i];
  delete [] M;
};


// Assignment

Matrix& Matrix::operator=(const double value)
{
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      M[i][j]=value;
  return (*this);
}

Matrix& Matrix::operator=(const Matrix& M_)
{
  int i,j;
  if (this == &M_) return (*this);
  // Allocate matrix if still empty
  if (M==0) {
    n = M_.n;
    m = M_.m;
    M = new double*[n];
    for (i=0; i<n; i++) M[i] = new double[m];
  };
  if ( (n!=M_.n) || (m!=M_.m) ) {
	cerr << "ERROR: Incompatible shapes in Matrix operator=(Matrix)" << endl;
    exit(1);
  };
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      M[i][j]=M_.M[i][j];
  return (*this);
}


// Component access

Vector Matrix::Col(int j) const
{
  Vector Res(n);
  for (int i=0; i<n; i++)
	Res.v[i]=M[i][j];
//	Res[i] = M[i][j];
  return Res;
}

Vector Matrix::Row(int i) const
{
  Vector Res(m);
  for (int j=0; j<m; j++)
	Res.v[j]=M[i][j];
//	Res[j] = M[i][j];
  return Res;
}

Vector Matrix::Diag() const
{
  if (n!=m) {
	cerr << "ERROR: Invalid shape in Matrix.Diag()" << endl;
	exit(1);
  };
  Vector Vec(n);
  for (int i=0; i<n; i++)
	Vec.v[i] = M[i][i];
	// Vec[i] = M[i][i];
  return Vec;
}


Matrix Matrix::slice(int first_row, int last_row, int first_col, int last_col)
{
  if (first_row<0 || last_row<first_row || n-1<last_row ||
	  first_col<0 || last_col<first_col || m-1<last_col) {
	cerr << "ERROR: Invalid arguments in Matrix.slice()" << endl;
    exit(1);
  };
  Matrix Aux(last_row-first_row+1,last_col-first_col+1);
  for (int i=0;i<=last_row-first_row;i++)
    for (int j=0;j<=last_col-first_col;j++)
       Aux(i,j) = M[i+first_row][j+first_col];
  return Aux;
}


void Matrix::SetCol(int j, const Vector& Col)
{
  if (Col.size()!=n) {
    cerr << "ERROR: Incompatible shapes in Matrix.SetCol()" << endl;
    exit(1);
  };
  if (j<0 || m<=j) {
    cerr << "ERROR: Column index out of range in Matrix.SetCol()" << endl;
    exit(1);
  };
  for (int i=0; i<n; i++) M[i][j]=Col(i);
}

void Matrix::SetRow(int i, const Vector& Row)
{
  if (Row.size()!=m) {
    cerr << "ERROR: Incompatible shapes in Matrix.SetRow()" << endl;
    exit(1);
  };
  if (i<0 || n<=i) {
    cerr << "ERROR: Row index out of range in Matrix.SetRow()" << endl;
    exit(1);
  };
  for (int j=0; j<n; j++) M[i][j]=Row(j);
}


// Unit matrix

Matrix Id(int Size)
{
  Matrix Aux(Size,Size);
  for (int i=0; i<Size; i++) Aux.M[i][i] = 1.0;
  return Aux;
}


// Diagonal matrix

Matrix Diag(const Vector& Vec)
{
  Matrix Mat(Vec.n,Vec.n);
  for (int i=0; i<Vec.n; i++) Mat.M[i][i] = Vec.v[i];
/*
  Matrix Mat(getn(Vec),getn(Vec));
  for (int i=0; i<getn(Vec); i++)
  {
	Mat.M[i][i] = Vec[i];;
  }
*/
  return Mat;
}


// Elementary rotations

Matrix R_x(double Angle)
{
  const double C = cos(Angle);
  const double S = sin(Angle);
  Matrix U(3,3);
  U.M[0][0] = 1.0;  U.M[0][1] = 0.0;  U.M[0][2] = 0.0;
  U.M[1][0] = 0.0;  U.M[1][1] =  +C;  U.M[1][2] =  +S;
  U.M[2][0] = 0.0;  U.M[2][1] =  -S;  U.M[2][2] =  +C;
  return U;
}

Matrix R_y(double Angle)
{
  const double C = cos(Angle);
  const double S = sin(Angle);
  Matrix U(3,3);
  U.M[0][0] =  +C;  U.M[0][1] = 0.0;  U.M[0][2] =  -S;
  U.M[1][0] = 0.0;  U.M[1][1] = 1.0;  U.M[1][2] = 0.0;
  U.M[2][0] =  +S;  U.M[2][1] = 0.0;  U.M[2][2] =  +C;
  return U;
}

Matrix R_z(double Angle)
{
  const double C = cos(Angle);
  const double S = sin(Angle);
  Matrix U(3,3);
  U.M[0][0] =  +C;  U.M[0][1] =  +S;  U.M[0][2] = 0.0;
  U.M[1][0] =  -S;  U.M[1][1] =  +C;  U.M[1][2] = 0.0;
  U.M[2][0] = 0.0;  U.M[2][1] = 0.0;  U.M[2][2] = 1.0;
  return U;
}


// Transposition

Matrix Transp(const Matrix& Mat)
{
  Matrix T(Mat.m,Mat.n);
  for ( int i=0; i<T.n; i++ )
    for ( int j=0; j<T.m; j++ )
      T.M[i][j] = Mat.M[j][i];
  return T;
}


// Inverse

Matrix Inv(const Matrix& Mat)
{
  const int n = Mat.n;

  Matrix LU(n,n), Inverse(n,n);
  Vector b(n), Indx(n);

  if (Mat.m!=Mat.n) {
    cerr << "ERROR: Invalid shape in Inv(Matrix)" << endl;
    exit(1);
  };

  // LU decomposition

  LU = Mat;
  LU_Decomp ( LU, Indx );

  // Solve Ax=b for  unit vectors b_1..b_n

  for (int j=0; j<n; j++ ) {
    b=0.0; b(j)= 1.0;                     // Set b to j-th unit vector
    LU_BackSub ( LU, Indx, b );           // Solve Ax=b
    Inverse.SetCol(j,b);                  // Copy result
  };

  return Inverse;

}


// Scalar multiplication and division of a vector

Matrix operator * (double value, const Matrix& Mat)
{
  Matrix Aux(Mat.n,Mat.m);
  for (int i=0; i<Mat.n; i++)
    for (int j=0; j<Mat.m; j++)
      Aux.M[i][j]=value*Mat.M[i][j];
  return Aux;
}

Matrix operator * (const Matrix& Mat, double value)
{
  return value*Mat;
}

Matrix operator / (const Matrix& Mat, double value)
{
  Matrix Aux(Mat.n,Mat.m);
  for (int i=0; i<Mat.n; i++)
    for (int j=0; j<Mat.m; j++)
      Aux.M[i][j]=Mat.M[i][j]/value;
  return Aux;
}


// Unary minus

Matrix operator - (const Matrix& Mat)
{
  Matrix Aux(Mat.n,Mat.m);
  for (int i=0; i<Mat.n; i++)
    for (int j=0; j<Mat.m; j++)
      Aux.M[i][j]=-Mat.M[i][j];
  return Aux;
}


// Matrix addition and subtraction

Matrix operator + (const Matrix& left, const Matrix& right)
{
  if ( (left.n!=right.n) || (left.m!=right.m) ) {
    cerr << "ERROR: Incompatible shape in +(Matrix,Matrix)" << endl;
    exit(1);
  };
  Matrix Aux(left.n,left.m);
  for (int i=0; i<left.n; i++)
    for (int j=0; j<left.m; j++)
      Aux.M[i][j] = left.M[i][j] + right.M[i][j];
  return Aux;
}

Matrix operator - (const Matrix& left, const Matrix& right)
{
  if ( (left.n!=right.n) || (left.m!=right.m) ) {
    cerr << "ERROR: Incompatible shape in -(Matrix,Matrix)" << endl;
    exit(1);
  };
  Matrix Aux(left.n,left.m);
  for (int i=0; i<left.n; i++)
    for (int j=0; j<left.m; j++)
      Aux.M[i][j] = left.M[i][j] - right.M[i][j];
  return Aux;
}


// Matrix product

Matrix operator * (const Matrix& left, const Matrix& right)
{
  if (left.m!=right.n) {
    cerr << "ERROR: Incompatible shape in *(Matrix,Matrix)" << endl;
    exit(1);
  };
  Matrix Aux(left.n,right.m);
  double Sum;
  for (int i=0; i<left.n; i++)
    for (int j=0; j<right.m; j++) {
      Sum = 0.0;
      for (int k=0; k<left.m; k++)
        Sum += left.M[i][k] * right.M[k][j];
      Aux.M[i][j] = Sum;
    }
  return Aux;
}


// Vector/matrix product

Vector operator * (const Matrix& Mat, const Vector& Vec)
{
  if (Mat.m!=Vec.n) {
    cerr << "ERROR: Incompatible shape in *(Matrix,Vector)" << endl;
    exit(1);
  };
  Vector Aux(Mat.n);
  double Sum;
  for (int i=0; i<Mat.n; i++) {
    Sum = 0.0;
    for (int j=0; j<Mat.m; j++)
      Sum += Mat.M[i][j] * Vec.v[j];
	Aux.v[i] = Sum;
  }
  return Aux;
}

Vector operator * (const Vector& Vec, const Matrix& Mat)
{
  if (Mat.n!=Vec.n) {
	cerr << "ERROR: Incompatible shape in *(Vector,Matrix)" << endl;
	exit(1);
  };
  Vector Aux(Mat.m);
  double Sum;
  for (int j=0; j<Mat.m; j++) {
	Sum = 0.0;
	for (int i=0; i<Mat.n; i++)
	  Sum += Vec.v[i] * Mat.M[i][j];
	Aux.v[j] = Sum;
  }
  return Aux;
}


// Dyadic product

Matrix Dyadic (const Vector& left, const Vector& right)
{
  Matrix Mat(left.n,right.n);
  for (int i=0;i<left.n;i++)
	for (int j=0;j<right.n;j++)
	  Mat.M[i][j] = left.v[i]*right.v[j];
  return Mat;
}


// Matrix output

ostream& operator << (ostream& os, const Matrix& Mat)
{
  int w = os.width();
  for (int i=0; i<Mat.size1(); i++) {
	for (int j=0; j<Mat.size2(); j++)
	  os << setw(w) << Mat(i,j);
	os << endl;
  }
  return os;
}


