
function Ematrices = Ematrix5pt_v3 (q, qp)

   % Computes possible E-matrices from matches
   % Call as follows:  E = Ematrix5-t(q1, q2)
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
   [E, A, poly] = partsolveE(q, qp);

   % Now, solve the polynomial by companion matrices. Take the real eigenvalues
   e = eig(companion(poly));
   I = find(not(imag(e)));
   e = e(I);

   s = size(e, 1);	% Number of solutions

   % Initialize the array of matrices
   Ematrices=zeros(3, 3, 0);

   % Compute the solutions one by one
   for i=1:s

      % Evaluate the matrix A at the values w
      w = e(i);
      M = w^4*A(:,:,5) + w^3*A(:,:,4) + w^2*A(:,:,3) + w*A(:,:,2) + A(:,:,1);

      % Find the minimum eigenvalue
      [X, f] = eig(M);
      ff = diag(f);
      ff = ff.*conj(ff);
      I = find(ff == min(ff));
      xy = X(:,I);

      % Extract solutions for x and y
      x = xy(2) / xy(1);
      y = xy(3) / xy(1);
      
      % Reconsitute the matrix
      Ematrices(:,:,i) = ...
	 w * E(:,:,1) + x * E(:,:,2) + y * E(:,:,3) + E(:,:,4);
   end

%-----------------------------------------------------------------------------

function M = companion(pol)

   s = size(pol);
   degree = s(2)-1;

   M = [[zeros(1, degree-1) ; eye(degree-1)] , -pol(1:degree)'/pol(degree+1)];

