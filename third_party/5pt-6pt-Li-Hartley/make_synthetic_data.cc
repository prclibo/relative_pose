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

#include "stdio.h"
#include "Utilities_CC/utilities_cc.h"

void get_synthetic_data (double q[6][3], double qp[6][3], double *f)
   {

   int npoints = 6;
   double depth0 = 5.0;
   double depth1 = 6.0;
   double flen = urandom(500.0, 100.0);
   printf ("Flength = %f\n", flen);
   *f = flen;

   // Get two camera matrices
   rhMatrix M0 = rhMatrix::I0();
   rhMatrix M1(3, 4);

   // Get a random rotation
   rhVector r(grandom(0.0, 0.1), grandom(0.0, 1.0), grandom(0.0, 1.0));
   rhMatrix R = Pmap::PvectorToMatrix(r);
   for2Dindex (i, j, R)
      M1[i][j] = R[i][j];

   // Place them looking at the origin
   M1[0][3] = M0[0][3] = 0.0;
   M1[1][3] = M0[1][3] = 0.0;
   M0[2][3] = depth0;
   M1[2][3] = depth1;

   // Now, multiply by focal length
   for (int i=0; i<4; i++)
      {
      M0[0][i] *= flen;
      M0[1][i] *= flen;
      M1[0][i] *= flen;
      M1[1][i] *= flen;
      }

   // Make points in a unit ball
   for (int i=0; i<npoints; i++)
      {
      // Get points in a unit square
      double x = 2* rh_irandom () - 1.0;
      double y = 2* rh_irandom () - 1.0;
      double z = 2* rh_irandom () - 1.0;
      rhVector xyz(x, y, z, 1.0);

      // Work out the projections
      rhVector u0 = (M0*xyz);
      rhVector u1 = (M1*xyz);

      for (int j=0; j<3; j++)
         {
         q[i][j] = u0[j];
         qp[i][j] = u1[j];
         }
      }
   }
