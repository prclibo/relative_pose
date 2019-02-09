%
%  This compares my various methods for 5 point E-matrix estimation
%
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



% Set up state for the random number generator
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;

nsamples = 1000;

%-------------------------------------------------------

reps = 100;

disp (' ');
disp (['My Mex method with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 5);
   Q2 = rand(3, 5);

   for rept=1:reps
      E = solveE(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end
end

disp (['total solutions = ', num2str(totalsolutions)]);
toc

%-------------------------------------------------------

reps = 100;

disp (' ');
disp (['My version of Nister Mex method with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 5);
   Q2 = rand(3, 5);

   for rept=1:reps
      E = solveE_nister(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end
end

disp (['total solutions = ', num2str(totalsolutions)]);
toc


%-------------------------------------------------------

reps = 10;

disp (' ');
disp (['My 5pt Matlab method with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 5);
   Q2 = rand(3, 5);

   for rept=1:reps
      E = Ematrix5pt(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end
end

disp (['total solutions = ', num2str(totalsolutions)]);
toc

%-------------------------------------------------------

% Increase the number of repetitions
reps = 100;

disp (' ');
disp (['My hybrid method using polyeig with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 5);
   Q2 = rand(3, 5);

   for rept=1:reps
      E = Ematrix5pt_v2(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end
end

disp (['total solutions = ', num2str(totalsolutions)]);
toc

%-------------------------------------------------------
reps = 100;

disp (' ');
disp (['My hybrid method using companion-matrix solver with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 5);
   Q2 = rand(3, 5);

   for rept=1:reps
      E = Ematrix5pt_v3(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end
end

disp (['total solutions = ', num2str(totalsolutions)]);
toc

