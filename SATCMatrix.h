//---------------------------------------------------------------------------

#ifndef SATCMatrixH
#define SATCMatrixH
//---------------------------------------------------------------------------

//------------------------------------------------------------------------------
//
// Matrix (class definition)
//
// Purpose:
//
//   Matrix data type and associated operations
//
//------------------------------------------------------------------------------

#include "SATCVector.h"

class Matrix
{
  private:

    // Elements
	int      n;                       // First dimension (number of rows)
	int      m;                       // Second dimension (number of columns)
	double **M;                       // Matrix M(n,m)

  public:

	// Constructors
	Matrix ();                                      // Matrix without elements
	Matrix (int dim1, int dim2);                    // Nullmatrix
	Matrix (const Matrix& M_);                      // Matrix copy
	Matrix (const double* p, int dim1, int dim2);   // Array copy

	// Destructor
	~Matrix();

	// Assignment
	Matrix& operator=(const double value);
	Matrix& operator=(const Matrix& M_);

	// Size
	int size1() const { return n; };
	int size2() const { return m; };

	// Component access (Fortran notation)
	double  operator () (int i, int j) const { return M[i][j]; };
	double& operator () (int i, int j)       { return M[i][j]; };
    Vector Col(int j) const;
    Vector Row(int i) const;
    Vector Diag() const;
    Matrix slice(int first_row, int last_row, int first_col, int last_col);
    void SetCol(int j, const Vector& Col);
    void SetRow(int i, const Vector& Row);

    // Unit matrix
    friend Matrix Id(int Size);

    // Diagonal matrix
    friend Matrix Diag(const Vector& Vec);

    // Elementary rotations
    friend Matrix R_x(double Angle);
    friend Matrix R_y(double Angle);
    friend Matrix R_z(double Angle);

    // Transposition and inverse
    friend Matrix Transp(const Matrix& Mat);
    friend Matrix Inv(const Matrix& Mat);

    // Scalar multiplication and division of a matrix
    friend Matrix operator * (double value, const Matrix& Mat);
    friend Matrix operator * (const Matrix& Mat, double value);
    friend Matrix operator / (const Matrix& Mat, double value);

    // Unary minus
    friend Matrix operator - (const Matrix& Mat);

    // Matrix addition and subtraction
    friend Matrix operator + (const Matrix& left, const Matrix& right);
    friend Matrix operator - (const Matrix& left, const Matrix& right);

    // Matrix product
    friend Matrix operator * (const Matrix& left, const Matrix& right);

    // Vector/matrix product
    friend Vector operator * (const Matrix& Mat, const Vector& Vec);
    friend Vector operator * (const Vector& Vec, const Matrix& Mat);

    // Dyadic product
    friend Matrix Dyadic (const Vector& left, const Vector& right);

    // Output
    friend std::ostream& operator << (std::ostream& os, const Matrix& Mat);

};

#endif
