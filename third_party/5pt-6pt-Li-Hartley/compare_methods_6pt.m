%
%  Compares various 6-point algorithms
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

reset = 1;

% Set up state for the random number generator
if (reset == 1)
   defaultStream = RandStream.getDefaultStream;
   savedState = defaultStream.State;

   for samples=1:78
      % Get some random values
      Q1 = rand(3, 6);
      Q2 = rand(3, 6);
   end

savedState = defaultStream.State;
end

nsamples = 100;
reps = 100;

%-------------------------------------------------------

reps = 100;

disp (' ');
disp (['My 6pt Matlab method with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 6);
   Q2 = rand(3, 6);

   for rept=1:reps
      [E, f] = Ematrix6pt(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end

   % check_solutions_6pt(E, f, Q1, Q2, '6pt_matlab.out');
end

disp (['total solutions = ', num2str(totalsolutions)]);

toc

%-------------------------------------------------------

reps = 100;

disp (' ');
disp (['My hybrid 6pt polyeig method with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 6);
   Q2 = rand(3, 6);

   for rept=1:reps
      [E, f] = Ematrix6pt_v2(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end

   % check_solutions_6pt(E, f, Q1, Q2, '6pt_matlab.out');
end

disp (['total solutions = ', num2str(totalsolutions)]);

toc


%-------------------------------------------------------
reps = 100;

disp (' ');
disp (['My 6point Mex method with ', num2str(reps*nsamples), ' iterations']);

% Go back to the same state
defaultStream.State = savedState;

tic
totalsolutions = 0;

for samples=1:nsamples

   % Get some random values
   Q1 = rand(3, 6);
   Q2 = rand(3, 6);

   for rept=1:reps
      [E, f] = solveE(Q1, Q2);
      nsolns = size(E, 3);
      totalsolutions = totalsolutions + nsolns;

      %for i=1:nsolns
      %   E(:,:,i)/E(3,3,i);
      %end
   end

   % check_solutions_6pt(E, f, Q1, Q2, '6pt_mex.out');
end

disp (['total solutions = ', num2str(totalsolutions)]);

toc
