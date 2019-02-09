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

//
// Matlab usage:
//
//    [E, A] = partsolveE6 (q1, q2)
//
// where
//    q1, q2 are point matches in two views (dimension 3 x 6)
//    E has dimension 3 x 3 x 3 containing a basis of 3 E-matrices
//    A has dimension 3 x 10 x 10 representing a 10 x 10 matrix of polynomials
//      of degree 2.
//
// method used is the one described by Li and Hartley
//

#include "string.h"
#include "mex.h"
#include "hidden6.h"

typedef double Matches6x3[6][3];
void compute_E_A_6pt ( 
	Matches q, Matches qp, 
        double E[3][3][3], double A[3][10][10]);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray *prhs[])
   {
   // Reconstruction using the trifocal tensor

   // -------------------------------------------------------------------------
   // Checking the arguments.

   // Check right number of arguments
   if (nrhs  != 2)
      {
      mexPrintf ("solveE: Wrong number of arguments.\n");
      mexPrintf ("Usage: E = solveE(q1, q2)\n");
      return; 
      }

   // Check the input
   mxArray *q1 = prhs[0];
   mxArray *q2 = prhs[1];

   // Check the dimensions of the arguments
   int ndimensions1 = mxGetNumberOfDimensions (q1);
   const mwSize *q1dim = mxGetDimensions (q1);

   int ndimensions2 = mxGetNumberOfDimensions (q2);
   const mwSize *q2dim = mxGetDimensions (q2);

   // Now check them
   if (ndimensions1 != 2 || q1dim[0] != 3 || q1dim[1] != 6 ||
       ndimensions2 != 2 || q2dim[0] != 3 || q2dim[1] != 6)
      {
      mexPrintf ("Bad input to mex function solveE\n");
      mexPrintf ("Inputs q1 and q2 must have dimensions [3, 6]\n");
      return; 
      }

   // -------------------------------------------------------------------------
   // Read and reformat the input

   Matches6x3 q, qp;

   double *p1 = (double *) mxGetData(q1);
   memcpy (&(q[0][0]),  p1, sizeof(q));

   double *p2 = (double *) mxGetData(q2);
   memcpy (&(qp[0][0]), p2, sizeof(qp));

   // -------------------------------------------------------------------------
   // Do the computation
   double E[3][3][3];
   double A[3][10][10];
   compute_E_A_6pt (q, qp, E, A);

   // -------------------------------------------------------------------------
   // Return the results

   // First is the set of matrices
   if (nlhs > 0) 
      {
      // Return E
      int dims[3];
      dims[0] = 3;
      dims[1] = 3; 
      dims[2] = 3;
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy (mxGetData(plhs[0]), E, sizeof(E));
      }

   if (nlhs > 1) 
      {
      // Return A
      int dims[3];
      dims[0] = 10;
      dims[1] = 10; 
      dims[2] = 3;
      plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy (mxGetData(plhs[1]), A, sizeof(A));
      }
   }
