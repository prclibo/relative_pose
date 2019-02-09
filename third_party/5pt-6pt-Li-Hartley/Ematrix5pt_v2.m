
function EMatrices = Ematrix5pt_v2 (q, qp)

   % Computes possible E-matrices from matches
   % Call as follows:  E = Ematrix5pt_v2(q1, q2)
   %   where both q1 and q2 are of dimension (3, 5).
   % This version uses polyeig to solve.
   %   Copyright Richard Hartley: 2010

   %--------------------------------------------------------------------------
   % LICENSE INFORMATION
   %
   % 1.  For academic/research users:
   %
   % This program is free for academic/research purpose:   you can redistribute
   % it and/or modify  it under the terms of the GNU General Public License as 
   % published by the Free Software Foundation, either version 3 of the License,
   % or (at your option) any later version.
   %
   % Under this academic/research condition,  this program is distributed in 
   % the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
   % the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   % PURPOSE. See the GNU General Public License for more details.
   %
   % You should have received a copy of the GNU General Public License along 
   % with this program. If not, see <http://www.gnu.org/licenses/>.
   %
   % 2.  For commercial OEMs, ISVs and VARs:
   % 
   % For OEMs, ISVs, and VARs who distribute/modify/use this software 
   % (binaries or source code) with their products, and do not license and 
   % distribute their source code under the GPL, please contact NICTA 
   % (www.nicta.com.au), and NICTA will provide a flexible OEM Commercial 
   % License. 
   %
   %---------------------------------------------------------------------------

   % Use a C routine to get the E-matrix basis and the polynomial matrix
   [E, A] = partsolveE(q, qp);

   % Now, solve the polynomial eigenvalue problem.  You need to
   % ignore infinite solutions.
   [X, e] = polyeig (A(:,:,1), A(:,:,2), A(:,:,3), A(:,:,4), A(:,:,5));
   I = find(real(e) < Inf & not(imag(e)));
   e = e(I);
   X = X(:,I);

   % Multiply out w, x, y, z by the basis to get the solution
   s = size(e, 1);
   wxyz = [e'; X(2,:)./X(1,:); X(3,:)./X(1,:); ones(1, s)];

   EMatrices = reshape(reshape(E, 9, 4) * wxyz, 3, 3, s);
