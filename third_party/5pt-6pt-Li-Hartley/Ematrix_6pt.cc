// Copyright Richard Hartley, 2010
static const char *copyright = "Copyright Richard Hartley, 2010";

//--------------------------------------------------------------------------
// LICENSE INFORMATION
//
// 1.  For academic/research users:
//
// This program is free for academic/research purpose:   you can redistribute
// it and/or modify  it under the terms of the GNU General Public License as 
// published by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// Under this academic/research condition,  this program is distributed in 
// the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along 
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
// 2.  For commercial OEMs, ISVs and VARs:
// 
// For OEMs, ISVs, and VARs who distribute/modify/use this software 
// (binaries or source code) with their products, and do not license and 
// distribute their source code under the GPL, please contact NICTA 
// (www.nicta.com.au), and NICTA will provide a flexible OEM Commercial 
// License. 
//
//---------------------------------------------------------------------------

#include <vector>
using namespace std;
#include "math.h"
#include "hidden6.h"

#ifdef NO_TARGETJR
void test_E_matrix (const double E[3][3]) {}
#else
void test_E_matrix (const double E[3][3]);
#endif

// #define RH_DEBUG

// Some forward declarations
class poly3_1;
class poly3_2;
class poly3_3;
class EmatrixSet;

const int NMatches = 6;		// Number of point matches
const int Nvar = 9-NMatches;	// Number of free dof for E

typedef double Matches_6[6][3];
typedef double Polynomial[Maxdegree+1];
typedef double Ematrix[3][3];

// Define some variables that can be used in equations, etc
// These are in alphabetical order, assuming the order z x y
const int z_ = 0;
const int x_ = 1;
const int y_ = 2;

const int zz_ = 0;
const int zx_ = 1;
const int zy_ = 2;
const int xx_ = 3;
const int xy_ = 4;
const int yy_ = 5;

const int zzz_ = 0;
const int zzx_ = 1;
const int zzy_ = 2;
const int zxx_ = 3;
const int zxy_ = 4;
const int zyy_ = 5;
const int xxx_ = 6;
const int xxy_ = 7;
const int xyy_ = 8;
const int yyy_ = 9;

class poly3_3
   {
   protected :
      double A_[10];

   public :
      static const int size_;
      poly3_3 operator + (poly3_3);
      void operator += (poly3_3);
      poly3_3 operator - (poly3_3);
      poly3_3 operator * (double);
      double &operator [] (int i) {return A_[i]; }

      void print ()
         {
	 printf ("%12.3e\n%12.3e %12.3e\n", A_[zzz_], A_[zzx_], A_[zzy_]);
	 printf ("%12.3e %12.3e %12.3e\n",  A_[zxx_], A_[zxy_], A_[zyy_]);
	 printf ("%12.3e %12.3e %12.3e %12.3e\n", 
				  A_[xxx_], A_[xxy_], A_[xyy_], A_[yyy_]);
         }

      void clear() 
         { 
         for (int i=0; i<size_; i++) 
            A_[i] = 0.0;
         }

   };

class poly3_2 
   {
   protected :
      double A_[6];

   public :
      static const int size_;
      poly3_3 operator * (poly3_1);
      poly3_2 operator + (poly3_2);
      void operator += (poly3_2);
      poly3_2 operator - (poly3_2);
      double &operator [] (int i) { return A_[i]; }

      void clear() 
         { 
         for (int i=0; i<size_; i++) 
            A_[i] = 0.0;
         }

      void print ()
         {
	 printf ("%12.3e\n%12.3e %12.3e\n", A_[zz_], A_[zx_], A_[zy_]);
	 printf ("%12.3e %12.3e %12.3e\n",  A_[xx_], A_[xy_], A_[yy_]);
         }
   };

class poly3_1
   {
   protected :
      double A_[Nvar];

   public :
      static const int size_;

      // Constructors
      poly3_1(){};
      poly3_1 (double x, double y, double z)
	 { A_[0] = x; A_[1] = y; A_[2] = z; }
      ~poly3_1 () {};

      // Operators
      poly3_2 operator * (poly3_1);
      poly3_1 operator + (poly3_1);
      poly3_1 operator - (poly3_1);
      double &operator [] (int i) { return A_[i]; }

      void clear() 
         { 
         for (int i=0; i<size_; i++) 
            A_[i] = 0.0;
         }

      void print ()
         {
	 printf ("%12.3e\n%12.3e %12.3e\n", A_[z_], A_[x_], A_[y_]);
         }
   };

class EmatrixSet
   {
   protected :
      poly3_1 E_[3][3];

   public :

      EmatrixSet () {};
      ~EmatrixSet() {};

      poly3_1 &operator () (int i, int j) { return E_[i][j]; }

      void print ()
         {
         for (int i=0; i<Nvar; i++)
            {
            for (int j=0; j<3; j++)
               {
               for (int k=0; k<3; k++)
	          printf ("%12.3e ", E_[j][k][i]);
               printf ("\n");
               }
            printf ("\n");
            }
         }
   };

//=============================================================================
//           Various operators on the polynomial classes
//=============================================================================

const int poly3_1::size_ = 3;
const int poly3_2::size_ = 6;
const int poly3_3::size_ = 10;

poly3_2 poly3_1::operator * (poly3_1 p2)
   {
   poly3_1 &p1 = *this;
   poly3_2 prod;

   prod[zz_] = p1[z_]*p2[z_];
   prod[zx_] = p1[z_]*p2[x_] + p1[x_]*p2[z_];
   prod[zy_] = p1[z_]*p2[y_] + p1[y_]*p2[z_];
   prod[xx_] = p1[x_]*p2[x_];
   prod[xy_] = p1[x_]*p2[y_] + p1[y_]*p2[x_];
   prod[yy_] = p1[y_]*p2[y_];

   return prod;
   }

poly3_3 poly3_2::operator * (poly3_1 p2)
   {
   poly3_2 &p1 = *this;
   poly3_3 prod;

   prod[zzz_] = p1[zz_]*p2[z_];
   prod[zzx_] = p1[zz_]*p2[x_] + p1[zx_]*p2[z_];
   prod[zzy_] = p1[zz_]*p2[y_] + p1[zy_]*p2[z_];
   prod[zxx_] = p1[zx_]*p2[x_] + p1[xx_]*p2[z_];
   prod[zxy_] = p1[zx_]*p2[y_] + p1[zy_]*p2[x_] + p1[xy_]*p2[z_];
   prod[zyy_] = p1[zy_]*p2[y_] + p1[yy_]*p2[z_];
   prod[xxx_] = p1[xx_]*p2[x_];
   prod[xxy_] = p1[xx_]*p2[y_] + p1[xy_]*p2[x_];
   prod[xyy_] = p1[xy_]*p2[y_] + p1[yy_]*p2[x_];
   prod[yyy_] = p1[yy_]*p2[y_];

#ifdef RH_DEBUG
   printf ("In poly3_2 * poly3_1\n");
   printf ("poly3_2 = \n");
   p1.print();
   printf ("poly3_1 = \n");
   p2.print();
   printf ("poly3_2 * poly3_2 = \n");
   prod.print();
#endif

   return prod;
   }

poly3_3 poly3_3::operator * (double k)
   {
   poly3_3 &p1 = *this;
   poly3_3 prod;

   for (int i=0; i<size_; i++)
      prod[i] = p1[i]*k;

   return prod;
   }

poly3_3 poly3_3::operator + (poly3_3 p2)
   {
   poly3_3 &p1 = *this;
   poly3_3 sum;

   for (int i=0; i<size_; i++)
      sum[i] = p1[i] + p2[i];

   return sum;
   }

void poly3_3::operator += (poly3_3 p2)
   {
   poly3_3 &p1 = *this;

   for (int i=0; i<size_; i++)
      p1[i] += p2[i];
   }

poly3_3 poly3_3::operator - (poly3_3 p2)
   {
   poly3_3 &p1 = *this;
   poly3_3 dif;

   for (int i=0; i<size_; i++)
      dif[i] = p1[i] - p2[i];

   return dif;
   }

poly3_2 poly3_2::operator + (poly3_2 p2)
   {
   poly3_2 &p1 = *this;
   poly3_2 sum;

   for (int i=0; i<size_; i++)
      sum[i] = p1[i] + p2[i];

   return sum;
   }

void poly3_2::operator += (poly3_2 p2)
   {
   poly3_2 &p1 = *this;
   for (int i=0; i<size_; i++)
      p1[i] += p2[i];
   }

poly3_2 poly3_2::operator - (poly3_2 p2)
   {
   poly3_2 &p1 = *this;
   poly3_2 dif;

   for (int i=0; i<size_; i++)
      dif[i] = p1[i] - p2[i];

   return dif;
   }

poly3_1 poly3_1::operator + (poly3_1 p2)
   {
   poly3_1 &p1 = *this;
   poly3_1 sum;

   for (int i=0; i<size_; i++)
      sum[i] = p1[i] + p2[i];

   return sum;
   }

poly3_1 poly3_1::operator - (poly3_1 p2)
   {
   poly3_1 &p1 = *this;
   poly3_1 dif;

   for (int i=0; i<size_; i++)
      dif[i] = p1[i] - p2[i];

   return dif;
   }

//=============================================================================

poly3_3 polydet3 (EmatrixSet E)
   {
   // Takes the determinant of a polynomial
   poly3_3 det = 
      (E(1,1)*E(2,2) - E(2,1)*E(1,2)) * E(0,0) +
      (E(2,1)*E(0,2) - E(0,1)*E(2,2)) * E(1,0) +
      (E(0,1)*E(1,2) - E(1,1)*E(0,2)) * E(2,0);

#ifdef RH_DEBUG
   printf ("Det =\n");
   det.print();
#endif

   return det;
   }

static poly3_2 traceEEt (EmatrixSet E, int deg)
   {
   // Takes the trace of E E' -- returns a quadratic polynomial,
   // giving the values of the trace in degree deg.
   // Trace of product is the elementwise product of the elements

   poly3_2 tr;

   switch (deg)
      {
      case 0:
   	 tr = E(2,2) * E(2, 2);
         break;

      case 1:

         tr = E(0,2) * E(0, 2) + E(1,2) * E(1, 2) 
            + E(2,0) * E(2, 0) + E(2,1) * E(2, 1);
         break;

      case 2:

         tr = E(0,0) * E(0, 0) + E(0,1) * E(0, 1) 
            + E(1,0) * E(1, 0) + E(1,1) * E(1, 1);
         break;
      }

#ifdef RH_DEBUG
   printf ("Trace is:\n");
   tr.print();
#endif

   return tr;
   }

static void mono_coeff (poly3_3 B, PolyMatrix A, int n, int deg)
   {
   // Extracts the monomial coefficients in x and y (with z = 1) from
   // a cubic homogeneous polynomial. Returns 4 vectors (degrees 0 to 3 in w)

   for (int i=0; i<B.size_; i++)
      A[n][i][deg] = B[i];
   }

static void EEeqns (EmatrixSet E, PolyMatrix A, PolyDegree degrees)
   {
   //
   // Computes the equations that will be used to input to polyeig.
   //    void EEeqns(E, A)
   // where E has dimensions E(3, 3, 4).  The output is a matrix
   // of dimension A(4, 10, 10, where A(i, :, :) is the coeffient of w^{i-1}
   //

   // Makes all the equations from the essential matrix E

   // First of all, set the equations to zero
   memset (&(A[0][0][0]), 0, sizeof(PolyMatrix));

   // First equation is from the determinant
   mono_coeff (polydet3(E), A, 0, 0);

   // Other equations from the equation 2 E*E'*E - tr(E*E') E = 0
   // In the following loop, we compute EE'E(i,j) = sum_pq E(i,p)*E(q,p)*E(q,j)
   // The way this is done is optimized for speed.  We compute first the matrix
   // EE'(i, q) and then use this to accumulate EE'E(i, j)

   // Find the trace - this is a quadratic polynomial
   poly3_2 tr0 = traceEEt(E, 0);	// Degree 0 in w
   poly3_2 tr1 = traceEEt(E, 1);	// Degree 1 in w
   poly3_2 tr2 = traceEEt(E, 2);	// Degree 2 in w

   // Constant term
   for (int i=0; i<3; i++)
      {
      // An array of cubic polynomials, one for each j = 0 ... 2
      poly3_3 EEE_i[3];  // Will hold (EE'E)(i,j)
      for (int j=0; j<3; j++) EEE_i[j].clear();

      // Compute each EE'(i,q) = sum_p E(i,p) E(q,p)
      // Only term that counts is p=2, q=2

      poly3_2 EE_i2  = E(i,2)*E(2,2);
      for (int j=0; j<3; j++) EEE_i[j] += EE_i2 * E(2,j);

      // Now, EE'E(i,j) is computed for this i and all j
      // We can complete the computation of the coefficients from EE'E(i, j)
      for (int j=0; j<3; j++)
         mono_coeff(EEE_i[j]*2.0 - tr0*E(i,j), A, 3*i+j+1, 0);
      }

   // Term in w
   for (int i=0; i<3; i++)
      {
      // An array of cubic polynomials, one for each j = 0 ... 2
      poly3_3 EEE_i[3];  // Will hold (EE'E)(i,j)
      for (int j=0; j<3; j++) EEE_i[j].clear();

      // Compute each EE'(i,q) = sum_p E(i,p) E(q,p)
      // Only term that counts is p=2, q=2

      poly3_2 EE_i0  = E(i,2)*E(0,2);
      for (int j=0; j<3; j++) EEE_i[j] += EE_i0 * E(0,j);

      poly3_2 EE_i1  = E(i,2)*E(1,2);
      for (int j=0; j<3; j++) EEE_i[j] += EE_i1 * E(1,j);

      poly3_2 EE_i2  = E(i,0)*E(2,0) + E(i,1)*E(2,1);
      for (int j=0; j<3; j++) EEE_i[j] += EE_i2 * E(2,j);

      // Now, EE'E(i,j) is computed for this i and all j
      // We can complete the computation of the coefficients from EE'E(i, j)
      for (int j=0; j<3; j++)
         mono_coeff(EEE_i[j]*2.0 - tr1*E(i,j), A, 3*i+j+1, 1);
      }

   // Term in w^2
   for (int i=0; i<3; i++)
      {
      // An array of cubic polynomials, one for each j = 0 ... 2
      poly3_3 EEE_i[3];  // Will hold (EE'E)(i,j)
      for (int j=0; j<3; j++) EEE_i[j].clear();

      // Compute each EE'(i,q) = sum_p E(i,p) E(q,p)
      // Only term that counts is p=2, q=2

      poly3_2 EE_i0  = E(i,0)*E(0,0) + E(i,1)*E(0,1);
      for (int j=0; j<3; j++)
         EEE_i[j] += EE_i0 * E(0,j);

      poly3_2 EE_i1  = E(i,0)*E(1,0) + E(i,1)*E(1,1);
      for (int j=0; j<3; j++)
         EEE_i[j] += EE_i1 * E(1,j);

      // Now, EE'E(i,j) is computed for this i and all j
      // We can complete the computation of the coefficients from EE'E(i, j)
      for (int j=0; j<3; j++)
         mono_coeff(EEE_i[j]*2.0 - tr2*E(i,j), A, 3*i+j+1, 2);
      }

   // Return also the degrees of the equations
   for (int j=0; j<Ncols; j++)
      {
      degrees[0][j] = 0;	// Equations from determinant have no w
      for (int i=1; i<Nrows; i++)
         degrees[i][j] = 2;	// Other equations have degree 2
      }
   }

static void null_space_solve_6x9 (
	double A[NMatches][9], EmatrixSet &E)
   {
   // This will compute the set of solutions for the equations
   // Sweep out one column at a time, starting with highest column number

   // We do Gaussian elimination to convert M to the form M = [X | I]
   // Then the null space will be [-I | X].

   // For present, this is done without pivoting.  
   // Mostly, do not need to actually change right hand part (that becomes I)

   const int lastrow  = NMatches-1; // Last real row in the matrix
   const int firstcol = Nvar; // First column to do elimination to make I
   const int lastcol  = 8; 

   // First sweep is to get rid of the above diagonal parts
   for (int col=lastcol; col>firstcol; col--)  // No need to do first col
      {
      // Remove column col
      const int row = col-firstcol;	// Row to pivot around
      const double pivot = A[row][col];

      // Sweep out all rows up to the current one 
      for (int i=0; i<row; i++)
         {
         // This factor of the pivot row is to subtract from row i
         const double fac = A[i][col] / pivot;

         // Constant terms
         for (int j=0; j<col; j++)
            A[i][j] -= fac * A[row][j];
         }
      }

   // Now, do backward sweep to clear below the diagonal
   for (int col=firstcol; col<lastcol; col++) // No need to do lastcol
      {
      // Remove column col
      const int row = col-firstcol;	// Row to pivot around
      const double pivot = A[row][col];

      // Sweep out all rows up to the current one 
      for (int i=row+1; i<=lastrow; i++)
         {
         // This factor of the pivot row is to subtract from row i
         const double fac = A[i][col] / pivot;

         // Constant terms
         for (int j=0; j<firstcol; j++)
            A[i][j] -= fac * A[row][j];
         }
      }

   // Make this into a matrix of solutions
   // This code is specific to 6x9.  Change for 5x9
   double fac;
   E(0, 0) = poly3_1(1.0, 0.0, 0.0);		// z
   E(0, 1) = poly3_1(0.0, 1.0, 0.0);		// x
   E(0, 2) = poly3_1(0.0, 0.0, 1.0);		// y

   fac = -1.0/A[0][3];
   E(1, 0) = poly3_1(fac*A[0][0], fac*A[0][1], fac*A[0][2]);
   fac = -1.0/A[1][4];
   E(1, 1) = poly3_1(fac*A[1][0], fac*A[1][1], fac*A[1][2]);
   fac = -1.0/A[2][5];
   E(1, 2) = poly3_1(fac*A[2][0], fac*A[2][1], fac*A[2][2]);
   fac = -1.0/A[3][6];
   E(2, 0) = poly3_1(fac*A[3][0], fac*A[3][1], fac*A[3][2]);
   fac = -1.0/A[4][7];
   E(2, 1) = poly3_1(fac*A[4][0], fac*A[4][1], fac*A[4][2]);
   fac = -1.0/A[5][8];
   E(2, 2) = poly3_1(fac*A[5][0], fac*A[5][1], fac*A[5][2]);

   // Now make them orthogonal
   for (int i=0; i<Nvar; i++)
      {
      // Take the inner product of the others
      for (int j=0; j<i; j++)
         {
         // Form the inner product of i-th and j-th layers
         double val = 0.0;
         for (int k=0; k<3; k++)
            for (int l=0; l<3; l++)
               val += E(k,l)[i] * E(k,l)[j];

         // Subtract it from E(.,.)[i]
         for (int k=0; k<3; k++)
            for (int l=0; l<3; l++)
                E(k,l)[i] -= val * E(k,l)[j];
         }

      // Now, normalize this one to length 1
      double lensq = 0.0;
      for (int k=0; k<3; k++)
         for (int l=0; l<3; l++)
             lensq += E(k,l)[i] * E(k,l)[i];
      double fac = 1.0 / sqrt(lensq);

      for (int k=0; k<3; k++)
         for (int l=0; l<3; l++)
             E(k,l)[i] *= fac;
      }

//#define USE_TEST_VALUES
#ifdef  USE_TEST_VALUES

   // Put an artificial value in 
   E(0,0)[0] =  2; E(0,1)[0] =   4; E(0,2)[0] = -1;
   E(1,0)[0] =  4; E(1,1)[0] =   5; E(1,2)[0] = -8;
   E(2,0)[0] =  2; E(2,1)[0] = -11; E(2,2)[0] =  8;

   E(0,0)[1] =  0; E(0,1)[1] =  -1; E(0,2)[1] =  2;
   E(1,0)[1] =  1; E(1,1)[1] =   7; E(1,2)[1] =  1;
   E(2,0)[1] = -2; E(2,1)[1] =   6; E(2,2)[1] =  7;

   E(0,0)[2] =  2; E(0,1)[2] =  -3; E(0,2)[2] =  7;
   E(1,0)[2] =  1; E(1,1)[2] =  -3; E(1,2)[2] = -9;
   E(2,0)[2] =  4; E(2,1)[2] =   1; E(2,2)[2] = -9;

#endif
   }


// Forward declaration
void print_polymatrix (PolyMatrix A, PolyDegree degrees, int maxdegree);

static void compute_E_A(
	Matches q, Matches qp, 			// Input points
	EmatrixSet &E, 				// Basis for E-matrix (returned)
	PolyMatrix &A, PolyDegree &degrees	// Matrix of polynomials
	)
   {
   // Computes the E-matrix from match inputs

   // A matrix to solve linearly for the ematrix
   double M[NMatches][9];
   memset (&(M[0][0]), 0, sizeof (M));

   for (int i=0; i<NMatches; i++)
      {
      M[i][0] = qp[i][0]*q[i][0];
      M[i][1] = qp[i][0]*q[i][1]; 
      M[i][2] = qp[i][0]*q[i][2];
      M[i][3] = qp[i][1]*q[i][0];
      M[i][4] = qp[i][1]*q[i][1]; 
      M[i][5] = qp[i][1]*q[i][2]; 
      M[i][6] = qp[i][2]*q[i][0];
      M[i][7] = qp[i][2]*q[i][1];
      M[i][8] = qp[i][2]*q[i][2]; 
      }

   // Solve using null_space_solve to get a basis for E
   null_space_solve_6x9 (M, E);

// #  define RH_DEBUG
#  ifdef RH_DEBUG
      printf ("E = \n");
      E.print();

   // Check that this is right
   for (int m=0; m<Nvar; m++)
      {
      printf ("Matrix %d\n", m);
      
      for (int pt=0; pt<NMatches; pt++)
         {
         double val = 0.0;
         for (int i=0; i<3; i++)  for (int j=0; j<3; j++)
            val += qp[pt][i] * E(i,j)[m] * q[pt][j];

         printf ("Point %d : %12.4e\n", pt, val);
         }
      }
#  endif

// #define USE_TEST_DATA
#ifdef USE_TEST_DATA
   // Override by getting the input data
   E(1,1) = poly3_1 ( 0,  2,  2);
   E(1,1) = poly3_1 (-1, -3,  4);
   E(1,1) = poly3_1 ( 2,  7, -1);

   E(1,1) = poly3_1 ( 1,  1,  4);
   E(1,1) = poly3_1 ( 7, -3,  5);
   E(1,1) = poly3_1 ( 1, -9, -8);

   E(1,1) = poly3_1 (-2,  4,  2);
   E(1,1) = poly3_1 ( 6,  1, -11);
   E(1,1) = poly3_1 ( 7, -9,  8);
#endif

   // Now, get the polynomial equations in A
   EEeqns(E, A, degrees);

// #define RH_DEBUG
#  ifdef RH_DEBUG
#  undef RH_DEBUG
   print_polymatrix (A, degrees, 2); 
#  endif
   }

void compute_E_A_6pt (
        Matches q, Matches qp,
        double EE[3][3][3],
        double AA[3][10][10])
   {
   // This is used by the Matlab interface.
   // It takes the matches and returns the basis for the E-matrices (EE)
   // along with a 3x3 matrix of polynomials, which allows us to solve
   // for w.  

   // Get the matrix set
   PolyMatrix A;
   PolyDegree degrees;
   EmatrixSet E;
   compute_E_A (q, qp, E, A, degrees);
   // print_polymatrix (A, degrees, 2);

   // Finally, get the 10-th degree polynomial out of this
   // if (poly) compute_determinant (A, poly);

   // Now, copy to the simple arrays
   if (EE)
      for (int d=0; d<3; d++) for (int i=0; i<3; i++) for (int j=0; j<3; j++)
         EE[d][i][j] = E(i,j)[d];  // Do not transpose - we want Ematrices thus

   if (AA)
      for (int d=0; d<3; d++) for (int i=0; i<10; i++) for (int j=0; j<10; j++)
         AA[d][i][j] = A[j][i][d]; // Transpose
   }

void print_polymatrix (PolyMatrix A, PolyDegree degrees, int maxdegree)
   {
   // Print out the matrix
   printf ("Equation matrix\n");
   for (int degree=0; degree<=maxdegree; degree++)
      {
      printf ("A%1d = {", degree);
      for (int i=0; i<10; i++)
         {
         if (i != 0) printf (",");
         printf ("{");
         for (int j=0; j<10; j++)
            {
            if (j != 0) printf (",");
	    if (degree <= degrees[i][j])
               printf ("%.12f ", A[i][j][degree]);
            else 
               printf ("%.12f ", 0.0);
            }
         printf ("\n");
         printf ("}");
         }
      printf ("};\n");
      } 
   }

// Declaration of the function to find roots
int find_real_roots_sturm(
   double *p, int order, double *roots, int *nroots, bool non_neg = true);

static double urandom()
   {
   // Returns a real random between -1 and 1
   const int MAXRAND = 65000;
   return 2.0*((rand()%MAXRAND)/((double) MAXRAND) - 0.5);
   }

static inline double pval (double *p, int deg, double x)
   {
   // Evaluates a polynomial at a given point x.  Assumes deg >= 0
   double val = p[deg];
   for (int i=deg-1; i>=0; i--)
      val = x*val + p[i];
   return val;
   }

static void compute_E_matrix (
	EmatrixSet &Es, 
	PolyMatrix A, 
        PolyDegree deg,		// Degree of each entry in A
        int rows[Nrows], 
        double w, 
	Ematrix &E
	)
   {
   // Compute the essential matrix corresponding to this root from
   // the matrix of equations A, assumed to be in row-echelon form
   // as defined by the array rows.  

   double a10 = pval(A[rows[1]][0], deg[rows[1]][0], w);
   double a11 = pval(A[rows[1]][1], deg[rows[1]][1], w);
   double a20 = pval(A[rows[2]][0], deg[rows[2]][0], w);
   double a21 = pval(A[rows[2]][1], deg[rows[2]][1], w);
   double a22 = pval(A[rows[2]][2], deg[rows[2]][2], w);

   double x = -a10/a11;
   double y = -(a20 + x*a21) / a22;

// #define RH_DEBUG
#ifdef RH_DEBUG
#undef RH_DEBUG
   printf ("In c_E_m: %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
           a10, a11, a20, a21, a22, x, y);

#endif

   // Multiply out the solution to get the essential matrix
   for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
         {
         poly3_1 &p = Es(i, j);
         E[i][j] = p[z_] + x*p[x_] + y*p[y_];
         }
   }

void compute_E_matrices_6pt (
     Matches q, Matches qp, 
     Ematrix Ematrices[Maxdegree], 
     double *flengths,
     int &nroots)
   {
   // Declare and clear the matrix of equations

   // Get the matrix set
   EmatrixSet E;
   PolyMatrix A;
   PolyDegree degrees;
   compute_E_A(q, qp, E, A, degrees);

// #define PRINT_RESULTS
#ifdef PRINT_RESULTS
   printf ("Polymatrix\n");
   print_polymatrix (A, degrees, 2);
#endif

   // Go ahead and find the polynomial determinant
   int rows[Nrows];
   double scale_factor = 1.0;

# define PRE_PROCESS
# ifdef PRE_PROCESS
   det_preprocess_6pt (A, degrees, 3);
   do_scale (A, degrees, scale_factor, true);
   // printf ("Scaled:  scale_factor = %f\n", scale_factor);
   // printf ("\nPolymatrix after preprocessing\n");
   // print_polymatrix (A, degrees, 2);
#endif

   find_polynomial_determinant (A, degrees, rows);
   double *poly = A[rows[0]][0];	
   int poly_degree = degrees[rows[0]][0];

   // Find the positive real roots
#ifdef PRE_PROCESS
   double roots[Maxdegree];
   find_real_roots_sturm(poly, poly_degree, roots, &nroots);
#else
   double roots[Maxdegree];
   find_real_roots_sturm(poly+3, poly_degree-3, roots, &nroots);
#endif

   // Only accept roots that are beyond a threshold - also multiply by scale
   int goodroots = 0;
   for (int i=0; i<nroots; i++)
      if (roots[i] > 1.0e-8)
         roots[goodroots++] = roots[i];
   nroots = goodroots;

   // Put their square roots in the array flengths
   for (int i=0; i<nroots; i++)
      flengths[i] = sqrt(roots[i]*scale_factor);

   // Now, get the ematrices
   for (int i=0; i<nroots; i++) 
      compute_E_matrix (E, A, degrees, rows, roots[i], Ematrices[i]);

// #define PRINT_RESULTS
#ifdef PRINT_RESULTS
#undef PRINT_RESULTS
   // printf ("Polynomial[0] = %13.6e\n", poly[0]);
   printf ("Polynomial\n");
   for (int i=0; i<=poly_degree; i++)
      printf ("\t%14.6f\n", poly[i]/poly[0]);
      // printf ("\t%14.6e\n", poly[i]);
#endif

// #define PRINT_RESULTS
#ifdef PRINT_RESULTS
#undef PRINT_RESULTS
   // Print out the roots
   printf ("Roots\n");
   for (int i=0; i<nroots; i++)
      printf ("\t%14.6f\n", roots[i]);
#endif

// #define PRINT_RESULTS
#ifdef PRINT_RESULTS
#undef PRINT_RESULTS
   // Print out the essential matrices
   printf ("Ematrices\n");
   for (int m=0; m<nroots; m++)
      {
      const Ematrix &E = Ematrices[m];
      for (int i=0; i<3; i++)
         printf ("\t%12.5f  %12.5f  %12.5f\n", E[i][0], E[i][1], E[i][2]);
      printf ("\n");

      // Now, compute to see if it has worked
      printf ("Verify: ");
      for (int pt=0; pt<NMatches; pt++) 
         {
         double sum = 0.0;
         for (int i=0; i<3; i++) for (int j=0; j<3; j++)
            sum += qp[pt][i] * E[i][j] * q[pt][j];
         printf ("%11.3e ", sum);
         }
      printf ("\n\n");
      }
#endif
   }

#ifndef BUILD_MEX

static int accuracy_test_main (int argc, char *argv[])
   {
   // Declare the data structure for the point matches
   const int NRepetitions = 10000;

   // Set up a histogram
   int histogram[20];
   for (int i=0; i<20; i++) histogram[i] = 0;

   int nhistogram[Maxdegree+1];
   for (int i=0; i<=Maxdegree; i++) nhistogram[i] = 0;

   double maxerr = 0.0; // Holds the maximum error, for verification

   for (int rep=0; rep<NRepetitions; rep++)
      {
      Matches_6 q, qp;

      // Fill the matches
      for (int i=0; i<NMatches; i++) for (int j=0; j<3; j++)
         q[i][j] = urandom();

      for (int i=0; i<NMatches; i++) for (int j=0; j<3; j++)
         qp[i][j] = urandom();

// #define USE_TEST_INPUT
#ifdef USE_TEST_INPUT

q[0][0]=9.70144633e-001; q[0][1]=2.15169613e-001; q[0][2]=7.60343346e-001;
q[1][0]=5.84102850e-001; q[1][1]=4.02951578e-001; q[1][2]=5.10040338e-001;
q[2][0]=4.95644219e-001; q[2][1]=6.51368343e-001; q[2][2]=7.43706145e-001;
q[3][0]=3.01953469e-001; q[3][1]=8.96119352e-002; q[3][2]=8.25965482e-001;
q[4][0]=3.89587009e-001; q[4][1]=7.75304853e-001; q[4][2]=1.79369095e-001;
q[5][0]=1.09361210e-001; q[5][1]=9.05158338e-001; q[5][2]=8.76351004e-001;

qp[0][0]=9.99792564e-001; qp[0][1]=8.64254686e-001; qp[0][2]=3.68775223e-002;
qp[1][0]=5.44682442e-001; qp[1][1]=9.97616182e-001; qp[1][2]=5.11008154e-001;
qp[2][0]=8.73512176e-001; qp[2][1]=7.02214723e-002; qp[2][2]=9.87538737e-001;
qp[3][0]=9.22713430e-001; qp[3][1]=5.64260920e-001; qp[3][2]=4.31462054e-001;
qp[4][0]=3.37844565e-001; qp[4][1]=7.20723428e-001; qp[4][2]=1.36555667e-002;
qp[5][0]=3.74060193e-001; qp[5][1]=9.22685220e-001; qp[5][2]=5.46479648e-001;

        // Make sure we do not repeat
        rep = NRepetitions;

#endif

      // Now, compute the ematrices and focal lengths
      Ematrix Ematrices[Maxdegree];
      double flengths[Maxdegree];
      int nroots;
      compute_E_matrices_6pt (q, qp, Ematrices, flengths, nroots);

      for (int m=0; m<nroots; m++)
         {
         // Get the focal length
         double f = flengths[m];

         // Multiply it out
         Ematrix E;
         E[0][0] = Ematrices[m][0][0] * f * f;
         E[0][1] = Ematrices[m][0][1] * f * f;
         E[0][2] = Ematrices[m][0][2] * f;

         E[1][0] = Ematrices[m][1][0] * f * f;
         E[1][1] = Ematrices[m][1][1] * f * f;
         E[1][2] = Ematrices[m][1][2] * f;

         E[2][0] = Ematrices[m][2][0] * f;
         E[2][1] = Ematrices[m][2][1] * f;
         E[2][2] = Ematrices[m][2][2];

         // Now, test it
         test_E_matrix (E);
         }

      // Keep histogram
      nhistogram[nroots] += 1;

      // Now, compute to see if it has worked
      for (int m=0; m<nroots; m++)
         {
         Ematrix &E = Ematrices[m];

         // Normalize the E matrix
         double sumsq = 0.0;
         for (int i=0; i<3; i++) for (int j=0; j<3; j++)
            sumsq += E[i][j]*E[i][j];
         double fac = 1.0 / sqrt(sumsq);
         for (int i=0; i<3; i++) for (int j=0; j<3; j++)
            E[i][j] *= fac;

         // Normalize the matrix
         for (int pt=0; pt<NMatches; pt++) 
            {
            double sum = 0.0;
            for (int i=0; i<3; i++) for (int j=0; j<3; j++)
               sum += qp[pt][i] * E[i][j] * q[pt][j];

            if (fabs(sum) > maxerr) maxerr = fabs(sum);

            // Get the logarithm
            int llog = (int) (-log10(fabs(sum) + 1.0e-100));
            if (llog >= 20) llog = 19;
            if (llog < 0) llog = 0;
            histogram[llog] += 1;

            if (llog == 0) printf ("Sum = %12.3e\n", sum);
            }
         }
      }

   printf ("Maximum error = %13.5e\n", maxerr);
   for (int i=0; i<20; i++)
      printf ("%2d: %d\n", i, histogram[i]);

   printf ("Number of solutions\n");
   for (int i=0; i<=Maxdegree; i++)
      printf ("%2d: %d\n", i, nhistogram[i]);

   return 0;
   }

static int speed_test_main (int argc, char *argv[])
   {
   // Declare the data structure for the point matches
   const int NRepetitions = 100;
   const int NData = 1000;

   printf ("\nTest for %.3e computations\n", (double) NRepetitions*NData);

   for (int dat=0; dat<NData; dat++)
      {
      Matches_6 q, qp;

      // Fill the matches
      for (int i=0; i<NMatches; i++) for (int j=0; j<3; j++)
         q[i][j] = urandom();

      for (int i=0; i<NMatches; i++) for (int j=0; j<3; j++)
         qp[i][j] = urandom();

      for (int rep=0; rep<NRepetitions; rep++)
         {
         // Now, compute the ematrices and focal lengths
         Ematrix Ematrices[Maxdegree];
         double flengths[Maxdegree];
         int nroots;
         compute_E_matrices_6pt (q, qp, Ematrices, flengths, nroots);
         }
      }

   return 0;
   }

static int once_off_main (int argc, char *argv[])
   {
   // Declare the data structure for the point matches

   Matches_6 q, qp;

   // Fill the matches
   for (int i=0; i<NMatches; i++) for (int j=0; j<3; j++)
      q[i][j] = urandom();

   for (int i=0; i<NMatches; i++) for (int j=0; j<3; j++)
      qp[i][j] = urandom();

#define USE_TEST_INPUT
#ifdef USE_TEST_INPUT
#undef USE_TEST_INPUT

q[0][0]=9.70144633e-001; q[0][1]=2.15169613e-001; q[0][2]=7.60343346e-001;
q[1][0]=5.84102850e-001; q[1][1]=4.02951578e-001; q[1][2]=5.10040338e-001;
q[2][0]=4.95644219e-001; q[2][1]=6.51368343e-001; q[2][2]=7.43706145e-001;
q[3][0]=3.01953469e-001; q[3][1]=8.96119352e-002; q[3][2]=8.25965482e-001;
q[4][0]=3.89587009e-001; q[4][1]=7.75304853e-001; q[4][2]=1.79369095e-001;
q[5][0]=1.09361210e-001; q[5][1]=9.05158338e-001; q[5][2]=8.76351004e-001;

qp[0][0]=9.99792564e-001; qp[0][1]=8.64254686e-001; qp[0][2]=3.68775223e-002;
qp[1][0]=5.44682442e-001; qp[1][1]=9.97616182e-001; qp[1][2]=5.11008154e-001;
qp[2][0]=8.73512176e-001; qp[2][1]=7.02214723e-002; qp[2][2]=9.87538737e-001;
qp[3][0]=9.22713430e-001; qp[3][1]=5.64260920e-001; qp[3][2]=4.31462054e-001;
qp[4][0]=3.37844565e-001; qp[4][1]=7.20723428e-001; qp[4][2]=1.36555667e-002;
qp[5][0]=3.74060193e-001; qp[5][1]=9.22685220e-001; qp[5][2]=5.46479648e-001;

#endif

   // Print out the values
   printf ("q1 = [");
   for (int j=0; j<NMatches; j++)
      {
      if (j != 0) printf (";");
      printf ("...\n");
      printf ("%14.11f, %14.11f, %14.11f", q[j][0], q[j][1], q[j][2]);
      }
   printf ("];\n");

   printf ("q2 = [");
   for (int j=0; j<NMatches; j++)
      {
      if (j != 0) printf (";");
      printf ("...\n");
      printf ("%14.11f, %14.11f, %14.11f", qp[j][0], qp[j][1], qp[j][2]);
      }
   printf ("];\n\n");

   // Now, compute the ematrices and focal lengths
   Ematrix Ematrices[Maxdegree];
   double flengths[Maxdegree];
   int nroots;
   compute_E_matrices_6pt (q, qp, Ematrices, flengths, nroots);

   for (int m=0; m<nroots; m++)
      {
      // Get the focal length
      double f = flengths[m];

      // Multiply it out
      Ematrix E;
      E[0][0] = Ematrices[m][0][0];
      E[0][1] = Ematrices[m][0][1];
      E[0][2] = Ematrices[m][0][2] / f;

      E[1][0] = Ematrices[m][1][0];
      E[1][1] = Ematrices[m][1][1];
      E[1][2] = Ematrices[m][1][2] / f;

      E[2][0] = Ematrices[m][2][0] / f;
      E[2][1] = Ematrices[m][2][1] / f;
      E[2][2] = Ematrices[m][2][2] / (f * f);

      // Now, test it
      test_E_matrix (E);
      }

   // Now, print out the focal lengths
   printf ("Focal lengths = \n");
   for (int i=0; i<nroots; i++)
      printf ("\t%13.5e\n", flengths[i]);
   
   return 0;
   }

#ifndef NO_TARGETJR
static int synth_main (int argc, char *argv[])
   {
   // Declare the data structure for the point matches
   void get_synthetic_data (double q[6][3], double qp[6][3], double *f);

   const int NRepetitions = 10000;
   for (int rep = 0; rep < NRepetitions; rep++)
      {
      Matches_6 q, qp;
      double flen;
      get_synthetic_data (q, qp, &flen);

      // Now, compute the ematrices and focal lengths
      Ematrix Ematrices[Maxdegree];
      double flengths[Maxdegree];
      int nroots;
      compute_E_matrices_6pt (q, qp, Ematrices, flengths, nroots);

      for (int m=0; m<nroots; m++)
         {
         // Get the focal length
         double f = flengths[m];

         // Multiply it out
         Ematrix E;
         E[0][0] = Ematrices[m][0][0];
         E[0][1] = Ematrices[m][0][1];
         E[0][2] = Ematrices[m][0][2] / f;

         E[1][0] = Ematrices[m][1][0];
         E[1][1] = Ematrices[m][1][1];
         E[1][2] = Ematrices[m][1][2] / f;

         E[2][0] = Ematrices[m][2][0] / f;
         E[2][1] = Ematrices[m][2][1] / f;
         E[2][2] = Ematrices[m][2][2] / (f * f);

         // Now, test it
         test_E_matrix (E);
         }

      // Now, print out the focal lengths
      printf ("Focal lengths = \n");
      for (int i=0; i<nroots; i++)
         printf ("\t%13.5e\n", flengths[i]);

      // Find the maximum error
      double maxerr = 1.0e6;
      for (int i=0; i<nroots; i++)
         {
         double err = fabs( (flen - flengths[i])/flen);
         if (err < maxerr) maxerr = err;
         }
      printf ("Error=\t%12.5e\n", maxerr);
      }

   return 0;
   }
#endif

#ifdef MAIN_6PT

int main (int argc, char *argv[])
   {
#ifdef NO_TARGETJR
   return once_off_main(argc, argv);
#else
   return synth_main(argc, argv);
#endif
   }

#endif // MAIN_6PT

#endif // BUILD_MEX
