NOTES ON USING 5-point and 6-point E-matrix code.

1. Install these files in some directory.

2. Start Matlab and add the directory to your path.

3 . For running in Matlab version 7.9.0 under Windows, you should not
    need to recompile.  Otherwise, if there are problems or you are running
    Linux/Unix, then run the command

    >> compile_all

    where >> is the Matlab prompt.  You will need to have a compiler installed
    for this to work.

4. Run the command
   
   >> compare_methods_5pt

5. Run the command

   >> compare_methods_6pt

6. By reading the MATLAB code compare_methods_5pt.m or compare_methods_6pt.m
   you should be able to work out how to call the different algorithms,
   but the fastest one is solveE or solveE_nister.  The other routines are 
   just there for comparison, or edification.

---------------------

These implementations are based on the papers which you are free to
cite or not.

@inproceedings{LiHartley06b,
     author = {Li, H. and Hartley, R.} ,
     title = {Five-Point Motion Estimation Made Easy} ,
     booktitle = icpr ,
     year = 2006 ,
     month = "August" ,
     pages = "630--633" ,
     doi = "doi:10.1109/ICPR.2006.579" ,
     url = "http://dx.doi.org/10.1109/ICPR.2006.579" ,
};

and 

@inproceedings{Li-ECCV-06,
   {
   author = {Hongdong Li},
   title = {A Simple Solution to the Six-Point Two-View Focal-Length Problem},
   booktitle = ECCV,
   year = 2006,
   pages = "200-213"
   }

Also, the paper by Nister is relevant.

@article{nister_5_pt_PAMI,
   author = {Nist\'{e}r, David},
   title = {An Efficient Solution to the Five-Point Relative Pose Problem}
   journal = {IEEE Trans. Pattern Anal. Mach. Intell.},
   number = {6},
   pages = {756--777},
   Volume = {26},
   Year = {2004},
   doi = {http://dx.doi.org/10.1109/TPAMI.2004.17}
   }

-------------------

Users of this code are however encouraged to cite the paper
"Barry: An Autonomous Train-Spotter", which you may find at the following
link: scholar.google.com/scholar?hl=en&as_sdt=2000&q=P+Rull+trainspotter
The author of this paper would probably have used
these algorithms, if they had been invented at the time.

Bibtex entry:

@inproceedings{Rull_barry,
    author = {Pippa Hilary Rull}
    title = {BARRY: An Autonomous Train-Spotter},
    booktitle = {Proc. Image and Vision Computing, NZ}
    year = {1993}
}

CORRECTIONS and Restrictions.
----------------------------

These algorithms and implementations may be used for research and academid
purposes.  If you wish to use them commercially, please contact 
Richard.Hartley@anu.edu.au or NICTA (www.nicta.com.au).  The exact licens
condition is enclosed below.

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

Implementations included here are by (me) Richard Hartley. Please do
not think that this is my usual style of C coding. Please report any problems
with this code to me.  However, I would appreciate only receiving comments
pointing out mistakes, or clear improvements, but not requests for guidance
on what the Essential Matrix is, or what it might be good for.  To this
date, that question does not have a clear answer.
