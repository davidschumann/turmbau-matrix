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

#include "SATCVector.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

#include <cmath>
#include <iomanip>

#include "GNU_iomanip.h"

using namespace std;

//------------------------------------------------------------------------------
//
// Vector (class implementation)
//
// Purpose:
//
//   Vector data type and associated operations
//
//------------------------------------------------------------------------------


// Constructors, destructor

Vector::Vector ()                              // Vector without elements
  : n(0)
{
  v = 0;
}

Vector::Vector (int Size)                      // Creates null vector
 : n(Size)
{
  v = new double [Size];
  for (int i=0; i<Size; i++) v[i]=0.0;
}

Vector::Vector (const Vector& V)               // Vector copy
 : n(V.n)
{
  v = new double [V.n];
  for (int i=0; i<V.n; i++) v[i]=V.v[i];
}

Vector::Vector (const double* p, int N)        // Array copy
  : n(N)
{
  v = new double [N];
  for (int i=0; i<N; i++) v[i]=p[i];
}

Vector::Vector (double x, double y, double z)  // 3dim-Vector
  : n(3)
{
  v = new double [3];
  v[0]=x; v[1]=y; v[2]=z;
}

Vector::Vector (double x, double y, double z,   // 6dim-Vector
                double X, double Y, double Z)
  : n(6)
{
  v = new double [6];
  v[0]=x; v[1]=y; v[2]=z;
  v[3]=X; v[4]=Y; v[5]=Z;
}

Vector::~Vector()
{
  delete [] v;
}

// Component access

Vector Vector::slice (int first, int last) const
{
  Vector Aux(last-first+1);
  for (int i=first; i<=last; i++) Aux.v[i-first]=v[i];
  return Aux;
}


// Square root of vector elements

Vector Vector::Sqrt()
{
  Vector Aux(n);
  for (int i=0; i<n; i++) Aux.v[i]=sqrt(v[i]);
  return Aux;
}


// Assignment

Vector& Vector::operator=(const double value)
{
  for (int i=0; i<n; i++) v[i]=value;
  return (*this);
}

Vector& Vector::operator=(const Vector& V)
{
  if (this == &V) return (*this);
  // Allocate vector if still empty
  if (v==0) {
    n = V.n;
    v = new double [V.n];
  };
  // Check dimension
  if (n!=V.n) {
    cerr << "ERROR: Incompatible sizes in Vector operator=(Vector)" << endl;
    exit(1);
  };
  // Copy elements
  for (int i=0; i<n; i++) v[i]=V.v[i];
  return (*this);
}

// Concatenation

Vector Stack (const Vector& a, const Vector& b)
{
  int    i;
  Vector c(a.size()+b.size());
  for (i=0;i<a.size();i++) c(i)=a(i);
  for (i=0;i<b.size();i++) c(i+a.size())=b(i);
  return c;
}


// Vector from polar angles

Vector VecPolar (double azim, double elev, double r)
{
  return Vector(r*cos(azim)*cos(elev),r*sin(azim)*cos(elev),r*sin(elev));
}


// Vector addition/subtraction with assignment

void Vector::operator += (const Vector& V)
{
  if (n!=V.n) {
	cerr << "ERROR: Incompatible shape in Vector operator+=(Vector)" << endl;
	exit(1);
  };
  for (int i=0; i<n; i++) v[i]+=V.v[i];
}

void Vector::operator -= (const Vector& V)
{
  if (n!=V.n) {
	cerr << "ERROR: Incompatible shape in Vector operator-=(Vector)" << endl;
	exit(1);
  };
  for (int i=0; i<n; i++) v[i]-=V.v[i];
}


// Dot product, norm, cross product

double Dot (const Vector& left, const Vector& right)
{
  if (left.n!=right.n) {
	cerr << "ERROR: Incompatible shape in Dot(Vector,Vector)" << endl;
	exit(1);
  };
  double Sum = 0.0;
  for (int i=0; i<left.n; i++) Sum+=left.v[i]*right.v[i];
  return Sum;
}

double Norm (const Vector& V)
{
  return sqrt(Dot(V,V));
}

Vector Cross (const Vector& left, const Vector& right)
{
  if ( (left.n!=3) || (right.n!=3) ) {
	cerr << "ERROR: Invalid dimension in Cross(Vector,Vector)" << endl;
	exit(1);
  };
  Vector Result(3);
  Result.v[0] = left.v[1]*right.v[2] - left.v[2]*right.v[1];
  Result.v[1] = left.v[2]*right.v[0] - left.v[0]*right.v[2];
  Result.v[2] = left.v[0]*right.v[1] - left.v[1]*right.v[0];
  return Result;
}


// Scalar multiplication and division of a vector

Vector operator * (double value, const Vector& V)
{
  Vector Aux(V.n);
  for (int i=0; i<V.n; i++) Aux.v[i]=value*V.v[i];
  return Aux;
}

Vector operator * (const Vector& V, double value)
{
  return value*V;
}

Vector operator / (const Vector& V, double value)
{
  Vector Aux(V.n);
  for (int i=0; i<V.n; i++) Aux.v[i]=V.v[i]/value;
  return Aux;
}


// Negation of a vector (unary minus)

Vector operator - (const Vector& V)
{
  Vector Aux(V.n);
  for (int i=0; i<V.n; i++) Aux.v[i]=-V.v[i];
  return Aux;
}


// Vector addition and subtraction

Vector operator + (const Vector& left, const Vector& right)
{
  if (left.n!=right.n) {
	cerr << "ERROR: Incompatible shape in +(Vector,Vector)" << endl;
	exit(1);
  };
  Vector Aux(left.n);
  for (int i=0; i<left.n; i++) Aux.v[i]=left.v[i]+right.v[i];
  return Aux;
}

Vector operator - (const Vector& left, const Vector& right)
{
  if (left.n!=right.n) {
	cerr << "ERROR: Incompatible shape in -(Vector,Vector)" << endl;
	exit(1);
  };
  Vector Aux(left.n);
  for (int i=0; i<left.n; i++) Aux.v[i]=left.v[i]-right.v[i];
  return Aux;
}

// Vector output

ostream& operator << (ostream& os, const Vector& Vec)
{
  int w = os.width();
  for (int i=0; i<Vec.size(); i++)
	os << setw(w) << Vec(i);
  os << endl;
  return os;
}


/*
int getn(const Vector& V) {
  return(V.n);
};

*/
