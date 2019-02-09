
function [EMatrices, flengths, A] = Ematrix6pt (q, qp)

   % Computes possible E-matrices and focal lengths from six matches
   % Call as follows:  [E, f] = Ematrix6pt(q1, q2)
   %   where both q1 and q2 are of dimension (3, 6).
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
   flengths = zeros(0);

   % Now, form the matrix 
   M  = zeros(6,9); % get the six points
   for i = 1:6
       M (i,:) = [qp(1,i)*q(1,i), qp(1,i)*q(2,i), qp(1,i)*q(3,i), ...
                  qp(2,i)*q(1,i), qp(2,i)*q(2,i), qp(2,i)*q(3,i), ...
                  qp(3,i)*q(1,i) ,qp(3,i)*q(2,i), qp(3,i)*q(3,i) ...
                 ];
   end

   % Compute the null-space.  This will contain a basis for
   % the 3-dimensional space of E-matrices.
   N = reshape(null(M), 3, 3, 3);
   N = permute(N, [3, 1, 2]);

   %---------------------------------------------------
   % For testing only - compare with Mathematica file
   %
   % Fz = [2,  4, -1; 4,  5, -8;  2, -11,  8];
   % Fx = [0, -1,  2; 1,  7,  1; -2,   6,  7];
   % Fy = [2, -3,  7; 1, -3, -9;  4,   1, -9];
   % N = [Fz, Fx, Fy];
   % N = reshape(N, 3, 3, 3);
   % N = permute(N, [3, 1, 2]);
   %---------------------------------------------------

   % Now, get the equations from the null-space
   A = EEeqns6pt(N);

   % Now, solve the polynomial eigenvalue problem.  You need to
   % ignore infinite solutions.
   [X, e] = polyeig (A(:,:,1), A(:,:,2), A(:,:,3));
   I = find(real(e) < 10^6 & real(e) > 0.000001 & not(imag(e)));
   e = e(I);
   X = X(:,I);

   % Take the appropriate linear combination of the basis to get result
   s = size(e, 1);
   zxy = [ones(1, s); X(2,:)./X(1,:); X(3,:)./X(1,:)];

   EMatrices = reshape(reshape(N, 3, 9)' * zxy, 3, 3, s);
   flengths = sqrt(e);

   % Uncomment this if you want solution checked
   % check_solutions (EMatrices, flengths, q, qp);

%%-----------------------------------------------------------------------------

function result = check_solutions (EMatrices, flengths, q, qp)

   % 
   % Function for checking that the results are right
   %

   % First, check that the solutions fit the points
   nsolns = size(flengths, 1);

   for i=1:nsolns

      disp(['Testing solution ', num2str(i)]);
      disp('------------------');

      % Display the focal length
      f = (flengths(i))

      F = EMatrices(:,:,i);
      E = diag([f, f, 1]) * EMatrices(:,:,i) * diag([f, f, 1])

      % Display the residual -- should be near zero
      residual = diag (q' * F * qp)'

      % E-matrix should have two equal singular values, other one zero
      svals = svd(E)

   end

   result = 1;

%%-----------------------------------------------------------------------------

function A = EEeqns6pt (E)

   %
   % Computes the equations that will be used to input to polyeig.
   %    A = EEeqns(E)
   % where E has dimensions E(3, 3, 3).  The output is a matrix
   % of dimension (10, 10, 3), where A(:, :, i) is the coeffient of w^{i-1}
   %
   
   % Makes all the equations from the essential matrix E
   
   % Set up equation set of appropriate size
   A = zeros(10, 10, 3);

   % Find the trace - this is a quadratic polynomial
   tr = traceEwEtw (E);

   % First equation is from the determinant
   pdet = zeros(3, 3, 3, 3);
   pdet(:,:,:,1) = polydet(E);
   A(1,:,:) = mono_coeff_6pt (pdet);
  
   % Other equations from the equation 2 E*E'*E - tr(E*E') E = 0
   % In the following loop, we compute EE'E(i,j) = sum_pq E(i,p)*E(q,p)*E(q,j)
   % The way this is done is optimized for speed.  We compute first the matrix
   % EE'(i, q) and then use this to accumulate EE'E(i, j)

   % Basic way to hold a polynomial in w, x, y, z is as follows
   % A(a, b, c, i) holds the value of w^(i-1) * a * b * c,
   %               where a, b, c are one of x, y and z

   deg = [1, 1, 0];	% Degree in w of the diagonals of omega

   eqn = 2;  % Count on the next equation

   for i=1:3

      EEE_i = zeros(3, 3, 3, 3, 3);  % Will hold (EE'E)(i,j) - last index is j

      % Compute each EE'(i,q) = sum_p E(i,p) E(q,p)
      for q=1:3

         % Accumulate EE(i, q)
         EE_iq = zeros(3, 3, 3);
         for p=1:3
            d = deg(p)+deg(q)+1;
            EE_iq(:,:,d) = EE_iq(:,:,d) + E(:,i,p) * E(:,q,p)';
         end

         % Now, accumulate EEE(ij) = sum_q  EE'(i,q) * E(q, j)
         for j=1:3
            EEE_i(:,:,:,:,j) = EEE_i(:,:,:,:,j) + mvprod2 (EE_iq, E(:,q,j));
         end
      end

      % Now, EE'E(i,j) is computed for this i and all j
      % We can complete the computation of the coefficients from EE'E(i, j)
      for j=1:3
         A(eqn,:,:)=mono_coeff_6pt(2 * EEE_i(:,:,:,:,j) -mvprod2(tr, E(:,i,j)));
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
   prod = reshape(kron(A, x), 3, 3, 3);

%%-----------------------------------------------------------------------------

function prod = mvprod2(A, x)

   % Multiplies a quadratic by a linear polynomial
   prod (:,:,:,1) = reshape(kron(A(:,:,1), x), 3, 3, 3);
   prod (:,:,:,2) = reshape(kron(A(:,:,2), x), 3, 3, 3);
   prod (:,:,:,3) = reshape(kron(A(:,:,3), x), 3, 3, 3);

%%-----------------------------------------------------------------------------

function C = mono_coeff_6pt (B)

   % Extracts the monomial coefficients in x and y (with z = 1) from
   % a cubic homogeneous polynomial. Returns 3 vectors (degrees 0 to 2 in w)

   % The size of the array to be returned
   C = zeros(10, 3);

   % Make some constants to make the code easier to read

   % Linear variables
   z = 1;
   x = 2;
   y = 3;

   % Monomials
   xx  = 4;
   xy  = 5;
   yy  = 6;
   xxx = 7;
   xxy = 8;
   xyy = 9;
   yyy = 10;

   % Terms in w^0
   for w=1:3
      C(1,   w) = B(z, z, z, w);
      C(x,   w) = B(x, z, z, w) + B(z, x, z, w) + B(z, z, x, w);
      C(y,   w) = B(y, z, z, w) + B(z, y, z, w) + B(z, z, y, w);
      C(xx,  w) = B(x, x, z, w) + B(x, z, x, w) + B(z, x, x, w);   
      C(yy,  w) = B(y, y, z, w) + B(y, z, y, w) + B(z, y, y, w);   
      C(xy,  w) = B(x, y, z, w) + B(y, x, z, w) + ... 
                  B(x, z, y, w) + B(y, z, x, w) + ...
                  B(z, x, y, w) + B(z, y, x, w);
      C(xxx, w) = B(x, x, x, w);
      C(xxy, w) = B(x, x, y, w) + B(x, y, x, w) + B(y, x, x, w);
      C(xyy, w) = B(x, y, y, w) + B(y, x, y, w) + B(y, y, x, w);
      C(yyy, w) = B(y, y, y, w);
   end

%%-----------------------------------------------------------------------------

function C = mono_coeff (B)

   % Extracts the monomial coefficients in x and y (with z = 1) from
   % a cubic homogeneous polynomial. Returns 3 vectors (degrees 0 to 2 in w)

   % The size of the array to be returned
   C = zeros(10, 3);

   % Make some constants to make the code easier to read

   % Linear variables
   x = 2;
   y = 3;
   z = 1;

   % Monomials
   xx  = 4;
   xy  = 5;
   yy  = 6;
   xxx = 7;
   xxy = 8;
   xyy = 9;
   yyy = 10;

   % Terms in w^0
   w = 1;
   C(1,   w) = B(z, z, z);
   C(x,   w) = B(x, z, z) + B(z, x, z) + B(z, z, x);
   C(y,   w) = B(y, z, z) + B(z, y, z) + B(z, z, y);
   C(xx,  w) = B(x, x, z) + B(x, z, x) + B(z, x, x);   
   C(yy,  w) = B(y, y, z) + B(y, z, y) + B(z, y, y);   
   C(xy,  w) = B(x, y, z) + B(y, x, z) + ... 
               B(x, z, y) + B(y, z, x) + ...
               B(z, x, y) + B(z, y, x);
   C(xxx, w) = B(x, x, x);
   C(xxy, w) = B(x, x, y) + B(x, y, x) + B(y, x, x);
   C(xyy, w) = B(x, y, y) + B(y, x, y) + B(y, y, x);
   C(yyy, w) = B(y, y, y);

%%-----------------------------------------------------------------------------

function tr = traceEwEtw (E)
   
   % Takes the trace of E E' -- returns a quadratic polynomial
   % Trace of product is the elementwise product of the elements

   tr = zeros(3, 3, 3);	 % First subscript is degree in w.  Each is quadratic

   tr(:,:,3) = E(:,1,1) * E(:,1, 1)' ...
               + E(:,1,2) * E(:,1, 2)' ...
               + E(:,2,1) * E(:,2, 1)' ...
               + E(:,2,2) * E(:,2, 2)';

   tr(:,:,2) = E(:,2,3) * E(:,2, 3)' ...
               + E(:,3,1) * E(:,3, 1)' ...
               + E(:,3,2) * E(:,3, 2)' ...
               + E(:,1,3) * E(:,1, 3)';

   tr(:,:,1) = E(:,3,3) * E(:,3, 3)';
