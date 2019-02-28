//---------------------------------------------------------------------------

#ifndef SATCVectorH
#define SATCVectorH
//---------------------------------------------------------------------------

//------------------------------------------------------------------------------
//
// Vector (class definition)
//
// Purpose:
//
//   Vector data type and associated operations
//
//------------------------------------------------------------------------------

#include <iostream>

class Vector
{

//  private:
  public:
	// Elements
	int    n;      // Dimension
	double *v;     // Vector v(n)


  public:

	//friend class Matrix;

	// Constructors
	Vector ();                              // Vector without elements
	Vector (int Size);                      // Nullvector of specified size
	Vector (const Vector& V);               // Vector copy
	Vector (const double* p, int N);        // Array copy
	Vector (double x, double y, double z);  // 3dim-Vector
	Vector (double x, double y, double z,   // 6dim-Vector
			double X, double Y, double Z);

	// Destructor
	~Vector();

	// Assignment
	Vector& operator=(const double value);
	Vector& operator=(const Vector& V);

	// Component access (Fortran notation)
	double  operator () (int i) const { return v[i]; };
	double& operator () (int i)       { return v[i]; };

/*
		// Index-Operator
	double& operator[](int i) {
	  return v[i];
	};
	friend int getn(const Vector& V);
*/


	Vector slice (int first, int last) const;

	// Square root of vector elements
	Vector Sqrt();

	// Size
	int size() const { return n; };

	// Concatenation
	friend Vector Stack (const Vector& a, const Vector& b);

	// Vector from polar angles
	friend Vector VecPolar (double azim, double elev, double r=1.0);

	// Vector addition/subtraction with assignment
	void operator += (const Vector& V);
	void operator -= (const Vector& V);

	// Dot product, norm, cross product
	friend double Dot (const Vector& left, const Vector& right);
	friend double Norm (const Vector& V);
	friend Vector Cross (const Vector& left, const Vector& right);

	// Scalar multiplication and division of a vector
	friend Vector operator * (double value, const Vector& V);
	friend Vector operator * (const Vector& V, double value);
	friend Vector operator / (const Vector& V, double value);

	// Negation of a vector (unary minus)
	friend Vector operator - (const Vector& V);

	// Vector addition and subtraction
	friend Vector operator + (const Vector& left, const Vector& right);
	friend Vector operator - (const Vector& left, const Vector& right);

/*
	// Diagonal matrix
	friend Matrix Diag(const Vector& Vec);

	// Vector/matrix product
	friend Vector operator * (const Matrix& Mat, const Vector& Vec);
	friend Vector operator * (const Vector& Vec, const Matrix& Mat);

	// Dyadic product
	friend Matrix Dyadic (const Vector& left, const Vector& right);
*/
	// Output
	friend std::ostream& operator << (std::ostream& os, const Vector& Vec);
};


#endif
