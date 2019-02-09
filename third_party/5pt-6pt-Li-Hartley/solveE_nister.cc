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
//    E = solveE_nister (q1, q2)
//
// where
//    q1, q2 are matched points (each of dimension 3 x 5)
//    E is a list of essential matrices returned
//
// method used is the one described by Li and Hartley
//

#include "string.h"
#include "mex.h"
#include "hidden6.h"

typedef double Ematrix[3][3];

void compute_E_matrices (
     Matches q, Matches qp, Ematrix Ematrices[10], int &nroots,
     bool optimized = true);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray *prhs[])
   {
   // Fast way to compute the Essential matrix from 5 points

   // -------------------------------------------------------------------------
   // Checking the arguments.

   // Check right number of arguments
   if (nrhs  != 2)
      {
      mexPrintf ("solveE: Wrong number of arguments.\n");
      mexPrintf ("Usage: [E,f] = solveE_nister(q1, q2)\n");
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
   if (ndimensions1 != 2 || q1dim[0] != 3 || q1dim[1] != 5 ||
       ndimensions2 != 2 || q2dim[0] != 3 || q2dim[1] != 5)
      {
      mexPrintf ("Bad input to mex function solveE\n");
      mexPrintf (
	"Inputs q1 and q2 must have dimensions [3, 5]\n");
      return; 
      }

   if (q1dim[1] != q2dim[1])
      {
      mexPrintf ("Bad input to mex function solveE\n");
      mexPrintf ( "Inputs q1 and q2 must have same dimensions.\n");
      return; 
      }

   // -------------------------------------------------------------------------
   // Read and reformat the input
   int npoints = q1dim[1];

   double q[5][3], qp[5][3];

   double *p1 = (double *) mxGetData(q1);
   memcpy (&(q[0][0]),  p1, 3*npoints*sizeof(double));

   double *p2 = (double *) mxGetData(q2);
   memcpy (&(qp[0][0]),  p2, 3*npoints*sizeof(double));

   // -------------------------------------------------------------------------
   // Do the computation
   Ematrix Ematrices[Maxdegree];
   double flengths[Maxdegree];
   int nroots;
   compute_E_matrices (q, qp, Ematrices, nroots);

   // Return the results
   // Return the cameras MM
   if (nlhs > 0) 
      {
      // Create an array of the right size and fill with the two matrices
      int dims[3];
      dims[0] = 3;
      dims[1] = 3; 
      dims[2] = nroots;
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy (mxGetData(plhs[0]), Ematrices, nroots*sizeof(Ematrix));
      }

   // Return the cameras MM
   if (nlhs > 1) 
      {
      // Create an array of the right size and fill with the two matrices
      int dims[2];
      dims[0] = nroots;
      dims[1] = 1; 
      plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy (mxGetData(plhs[1]), flengths, nroots*sizeof(double));
      }
   }
