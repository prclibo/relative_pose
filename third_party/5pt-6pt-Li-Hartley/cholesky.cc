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


/* Solves sets of linear equations, in the symmetric case.  */

/* Author : Richard Hartley */

#include <math.h>
#include "hidden6.h"

/*****************************************************************/
/*      Solving linear equations based on Choleski               */
/*****************************************************************/

void backsolve_choleski (
	KMatrix K, double *x, int n, int bandwidth)
   {
   /* Assumes that K is in Choleski form.  In other works, only the
    * bottom half L of K is used. It is also assumed that only "bandwidth"
    * rows below the diagonal are non-zero.
    * Solves L L^T b  = x.
    * The solution, b overwrites x.
    */
   int i, j;

   /* First step of back substitution */
   for (i=0; i<n; i++)
      {
      int startj = i-bandwidth+1;
      if (startj < 0) startj = 0;

      for (j=startj; j<i; j++)
         x[i] -= K[i][j] * x[j];
      x[i] /= K[i][i];
      }

   /* Second set of back substitution */
   for (i=n-1; i>=0; i--)
      {
      int endj = i+bandwidth;
      if (endj >= n) endj = n-1;

      for (j=i+1; j<=endj; j++)
	 x[i] -= K[j][i] * x[j];
      x[i] /= K[i][i];
      }
   }

int choleski_decomposition (KMatrix A, int n, int bandwidth)
   {
   /* Does the choleski decomposition of a matrix.
    * It is assumed that A has bandwidth as given (i.e number of rows above
    * or on the diagonal 
    */
   int i, j, k;

   /* Do it row by row */
   for (j=0; j<n; j++)
      {
      /* Reduce the j-th column */

      // Last non-zero row + 1 in this column
      int endi = j+bandwidth;
      if (endi > n) endi = n;

      if (j > 0)
         {
	 for (i=j; i<endi; i++)
            {
            int startk = i-bandwidth+1;
            if (startk < 0) startk = 0;

	    for (k=startk; k<j; k++)
	       A[i][j] -= A[i][k]*A[j][k];
            }
         }

      /* Divide by the square root */
      if (A[j][j] > 0.0)
         {
	 /* Work out the scaling factor */
	 double d = sqrt(A[j][j]);
	 double fac = 1.0 / d;

	 /* Update them */
	 A[j][j] = d;
	 for (i=j+1; i<endi; i++)
	    A[i][j] *= fac;
	 }
      else
	 return 0;
      }

   /* Return success */
   return 1;
   }

int solve_choleski (KMatrix A, int n, double *x, int bandwidth)
   {
   /* Solves a symmetric set of equations */
   /* On return, A is replaced by a matrix such that L L^T = A.
    * and x holds the solution */
   /*
    * This routine does not seem to be as good as solve_symmetric.
    *   + It is a bit slower.
    *   + It requires positive definite matrices.
    */
 
   /* First, get the decomposition */
   if (! choleski_decomposition (A, n, bandwidth)) return 0;

   /* Now, do back substitution */
   backsolve_choleski (A, x, n, bandwidth);

   /* Return success */
   return 1;
   }

