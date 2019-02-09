
function EMatrices = Ematrix5pt (q, qp)

   % Computes possible E-matrices from matches
   % Call as follows:  E = Ematrix5-t(q1, q2)
   %   where both q1 and q2 are of dimension (3, 5).
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

   % Set up an array to keep the matrices
   EMatrices = zeros(3, 3, 0);

   % Now, form the matrix 
   M  = zeros(5,9); % get the five points
   for i = 1:5
       M (i,:) = [qp(1,i)*q(1,i), qp(1,i)*q(2,i), qp(1,i)*q(3,i), ...
                  qp(2,i)*q(1,i), qp(2,i)*q(2,i), qp(2,i)*q(3,i), ...
                  qp(3,i)*q(1,i) ,qp(3,i)*q(2,i), qp(3,i)*q(3,i) ...
                 ];
   end

   % Compute the null-space.  This will contain a basis for
   % the 4-dimensional space of E-matrices.
   N = reshape(null(M), 3, 3, 4);
   N = permute(N, [3, 1, 2]);

   %---------------------------------------------------
   % For testing only - compare with Mathematica file
   %
   % Fw = [2,  4, -1; 4,  5, -8;  2, -11,  8];
   % Fx = [0, -1,  2; 1,  7,  1; -2,   6,  7];
   % Fy = [2, -3,  7; 1, -3, -9;  4,   1, -9];
   % Fz = [5,  2,  7; 1, -2, -4;  5,  -1,  8];
   % N = [Fw, Fx, Fy, Fz]
   % N = reshape(N, 3, 3, 4)
   % N = permute(N, [3, 1, 2])
   %---------------------------------------------------

   % Now, get the equations from the null-space
   A = EEeqns(N);

   % Now, little trick - make smaller by reducing by inv(A1)
   Ainv = inv(A(:,:,1));

   A1 = Ainv * A(:,:,1);
   A2 = Ainv * A(:,:,2);
   A3 = Ainv * A(:,:,3);
   A4 = Ainv * A(:,:,4);

   A1 = A1(1:6, 1:6);
   A2 = A2(1:6, 1:6);
   A3 = A3(1:6, 1:6);
   A4 = A4(1:6, 1:6);

   %A1 = A(:,:,1);
   %A2 = A(:,:,2);
   %A3 = A(:,:,3);
   %A4 = A(:,:,4);

   % Now, solve the polynomial eigenvalue problem.  You need to
   % ignore infinite solutions.
   [X, e] = polyeig (A1, A2, A3, A4);
   I = find(real(e) < Inf & not(imag(e)));
   e = e(I);
   X = X(:,I);

   % Take the appropriate linear combination of the basis to get result
   s = size(e, 1);
   wxyz = [e'; X(2,:)./X(1,:); X(3,:)./X(1,:); ones(1, s)];

   EMatrices = reshape(reshape(N, 4, 9)' * wxyz, 3, 3, s);

%%-----------------------------------------------------------------------------

function A = EEeqns (E)

   %
   % Computes the equations that will be used to input to polyeig.
   %    A = EEeqns(E)
   % where E has dimensions E(3, 3, 4).  The output is a matrix
   % of dimension (10, 10, 4), where A(:, :, i) is the coeffient of w^{i-1}
   %
   
   % Makes all the equations from the essential matrix E
   
   % Set up equation set of appropriate size
   A = zeros(10, 10, 4);

   % Find the trace - this is a quadratic polynomial
   tr = traceEEt (E);

   % First equation is from the determinant
   A(1,:,:) = mono_coeff (polydet(E));
  
   % Other equations from the equation 2 E*E'*E - tr(E*E') E = 0
   % In the following loop, we compute EE'E(i,j) = sum_pq E(i,p)*E(q,p)*E(q,j)
   % The way this is done is optimized for speed.  We compute first the matrix
   % EE'(i, q) and then use this to accumulate EE'E(i, j)

   eqn = 2;  % Count on the next equation
   for i=1:3
      EEE_i = zeros(4, 4, 4, 3);  % Will hold (EE'E)(i,j)

      % Compute each EE'(i,q) = sum_p E(i,p) E(q,p)
      for q=1:3

         % Accumulate EE(i, q)
         EE_iq = zeros(4, 4);
         for p=1:3
	    EE_iq = EE_iq + E(:,i,p) * E(:,q,p)';
         end

         % Now, accumulate EEE(ij) = sum_q  EE'(i,q) * E(q, j)
         for j=1:3
            EEE_i(:,:,:,j) = EEE_i(:,:,:,j) + mvprod (EE_iq, E(:,q, j));
         end
      end

      % Now, EE'E(i,j) is computed for this i and all j
      % We can complete the computation of the coefficients from EE'E(i, j)
      for j=1:3
         A(eqn,:,:) = mono_coeff (2 * EEE_i(:,:,:,j) - mvprod(tr, E(:,i,j)));
         eqn = eqn+1;
      end
   end

%%-----------------------------------------------------------------------------

function det = polydet (E)

   % Takes the determinant of the essential matrix.  Returns a cubic
   % polynomial (that is, a 4 x 4 x 4 array).

   det =  ...
	mvprod (E(:,1,1), (E(:,2,2) * E(:,3,3)' - E(:,3,2) * E(:,2,3)')) + ...
	mvprod (E(:,2,1), (E(:,3,2) * E(:,1,3)' - E(:,1,2) * E(:,3,3)')) + ...
	mvprod (E(:,3,1), (E(:,1,2) * E(:,2,3)' - E(:,2,2) * E(:,1,3)'));

%%-----------------------------------------------------------------------------

function prod = mvprod(A, x)

   % Multiplies a quadratic by a linear polynomial
   prod = reshape(kron(A, x), 4, 4, 4);

%%-----------------------------------------------------------------------------

function C = mono_coeff (B)

   % Extracts the monomial coefficients in x and y (with z = 1) from
   % a cubic homogeneous polynomial. Returns 4 vectors (degrees 0 to 3 in w)

   % The size of the array to be returned
   C = zeros(10, 4);

   % Make some constants to make the code easier to read

   % Degrees of terms in w
   w0 = 1;
   w1 = 2;
   w2 = 3;
   w3 = 4;

   % Linear variables
   w = 1;
   x = 2;
   y = 3;
   z = 4;

   % Monomials
   xx  = 4;
   xy  = 5;
   yy  = 6;
   xxx = 7;
   xxy = 8;
   xyy = 9;
   yyy = 10;

   % Terms in w^0
   C(1,   w0) = B(z, z, z);
   C(x,   w0) = B(x, z, z) + B(z, x, z) + B(z, z, x);
   C(y,   w0) = B(y, z, z) + B(z, y, z) + B(z, z, y);
   C(xx,  w0) = B(x, x, z) + B(x, z, x) + B(z, x, x);   
   C(yy,  w0) = B(y, y, z) + B(y, z, y) + B(z, y, y);   
   C(xy,  w0) = B(x, y, z) + B(y, x, z) + ... 
                B(x, z, y) + B(y, z, x) + ...
                B(z, x, y) + B(z, y, x);
   C(xxx, w0) = B(x, x, x);
   C(xxy, w0) = B(x, x, y) + B(x, y, x) + B(y, x, x);
   C(xyy, w0) = B(x, y, y) + B(y, x, y) + B(y, y, x);
   C(yyy, w0) = B(y, y, y);

   % Terms in w^1
   C(1,   w1) = B(w, z, z) + B(z, w, z) + B(z, z, w);
   C(x,   w1) = B(x, w, z) + B(w, z, x) + B(z, x, w) + ...
                B(w, x, z) + B(x, z, w) + B(z, w, x);
   C(y,   w1) = B(y, w, z) + B(w, z, y) + B(z, y, w) + ...
                B(w, y, z) + B(y, z, w) + B(z, w, y);
   C(xx,  w1) = B(x, x, w) + B(x, w, x) + B(w, x, x);   
   C(yy,  w1) = B(y, y, w) + B(y, w, y) + B(w, y, y);   
   C(xy,  w1) = B(x, y, w) + B(y, x, w) + ... 
                B(x, w, y) + B(y, w, x) + ...
                B(w, x, y) + B(w, y, x);

   % Terms in w^2
   C(1,   w2) = B(w, w, z) + B(w, z, w) + B(z, w, w);
   C(x,   w2) = B(x, w, w) + B(w, x, w) + B(w, w, x);
   C(y,   w2) = B(y, w, w) + B(w, y, w) + B(w, w, y);

   % Terms in w^3
   C(1,   w3) = B(w, w, w);

%%-----------------------------------------------------------------------------

function tr = traceEEt (E)
   
   % Takes the trace of E E' -- returns a quadratic polynomial
   % Trace of product is the elementwise product of the elements

   tr = E(:,1,1) * E(:,1, 1)' ...
      + E(:,1,2) * E(:,1, 2)' ...
      + E(:,1,3) * E(:,1, 3)' ...
      + E(:,2,1) * E(:,2, 1)' ...
      + E(:,2,2) * E(:,2, 2)' ...
      + E(:,2,3) * E(:,2, 3)' ...
      + E(:,3,1) * E(:,3, 1)' ...
      + E(:,3,2) * E(:,3, 2)' ...
      + E(:,3,3) * E(:,3, 3)';
