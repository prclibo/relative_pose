//#define DEBUG_ERR

int const MAXSOLS = 20;
#define MAXIT 30 // maximal number of iterations
#define MAXTREE (2*MAXIT + 1) // maximal depth of the isolating tree
#define ACC 1.0e-12 // root is polished to an approximate accuracy ACC
#define DEG MAXSOLS // degree of the polynomial



// all real roots of a polynomial p(x) lie in (bound[0], bound[1])
// Kioustelidis' method
void getBounds (const double p[DEG], double bound[2])
{
	double M = fabs(p[0]);
	bound[0] = bound[1]=0.0;
	for (int i = 0; i < DEG; ++i)
	{
		const double c = fabs(p[i]);
		if (c > M) M = c;
		if (p[i] < 0)
		{
			const double t = pow(-p[i], 1.0/(double)(DEG - i));
			if (t > bound[1]) bound[1] = t;
			if (t > bound[0] && !(i%2)) bound[0] = t;
		}
		else if (i%2)
		{
			const double t = pow(p[i], 1.0/(double)(DEG - i));
			if (t > bound[0]) bound[0] = t;
		}
	} // end for i
	M += 1.0; // Cauchy bound, since p[DEG] = 1.0;
	bound[0] *= -2.0;
	if (bound[0] < -M) bound[0] = -M;
	bound[1] *= 2.0;
	if (bound[1] > M) bound[1] = M;

} // end getBounds()



// get data to compute the Sturm sequence for a polynomial p(x) at any point
void getQuotients (const double p[DEG], const double dp[DEG - 1], double q[DEG + 1][2])
{
	double r_[DEG], r[DEG - 1];
	
	r_[DEG - 1] = 1.0;
	q[DEG][1] = p[DEG - 1] - dp[DEG - 2];
	r_[0] = dp[0];
	r[0] = dp[0]*q[DEG][1] - p[0];
	for (int j = 1; j < DEG - 1; ++j)
	{
		r_[j] = dp[j];
		r[j] = dp[j]*q[DEG][1] + r_[j - 1] - p[j];
	}

	for (int i = DEG - 1; i >= 2; --i)
	{
		const int i1 = i - 1;

		#ifdef DEBUG_ERR
		{
			if (!r[i1])
				std::cout << "division by zero in getQuotients()\n";
		}
		#endif

		const double ri1 = 1.0/r[i1];
		q[i][0] = r_[i]*ri1;
		q[i][1] = (r_[i1] - q[i][0]*r[i - 2])*ri1;
		const double t = r_[0];
		r_[0] = r[0];
		r[0] = r[0]*q[i][1] - t;
		for (int j = 1; j < i1; ++j)
		{
			const double t = r_[j];
			r_[j] = r[j];
			r[j] = r[j]*q[i][1] + q[i][0]*r_[j - 1] - t;
		}
		r_[i1] = r[i1];
	} // end for i

	q[1][0] = r_[1];
	q[1][1] = r_[0];
	q[0][1] = r[0];

} // end getQuotients()



// evaluate polynomial p(x) at point x0
double evalPoly (const double p[DEG], const double &x0)
{
	double s = x0 + p[DEG - 1];
	for (int i = DEG - 2; i >= 0; --i)
		s = s*x0 + p[i];
	return s;

} // end evalPoly()



// compute the number of sign changes in a sequence seq[]
int nchanges (const double seq[DEG + 1])
{
	int s = 0, s1, s2 = (seq[0] > 0.0)? 1 : ((seq[0] < 0.0)? -1 : 0);
	for (int i = 1; i < DEG + 1 && seq[i]; ++i)
	{
		s1 = s2;
		s2 = (seq[i] > 0)? 1 : -1;
		if (!(s1 + s2)) ++s;
	}
	return s;

} // end nchanges()



// evaluate Sturm sequence at point a
int evalSturmSeq (const double q[DEG + 1][2], const double &a)
{
	double sa[DEG + 1];
	// initialize sa
	sa[0] = q[0][1];
	sa[1] = q[1][0]*a + q[1][1];
	// compute sa recursively
	for (int i = 2; i < DEG; ++i)
		sa[i] = (q[i][0]*a + q[i][1])*sa[i - 1] - sa[i - 2];
    sa[DEG] = (a + q[DEG][1])*sa[DEG - 1] - sa[DEG - 2]; // since q[DEG][0] = 1
	return nchanges(sa);

} // end evalSturmSeq()



// isolate all real roots of a given polynomial p(x)
int isolateRoots (const double p[DEG], const double dp[DEG - 1], double Isol[DEG][2])
{
	int nIsol = 0, nTree = 1, nIters = 1, min = 0, sTree[MAXTREE][2];
	double Tree[MAXTREE][2], q[DEG + 1][2];

	// initialize the tree
	// all real roots of the polynomial p(x) lie in (Tree[0][0], Tree[0][1])
	getBounds(p, Tree[0]);

	getQuotients(p, dp, q);
	sTree[0][0] = evalSturmSeq(q, Tree[0][0]);
	sTree[0][1] = evalSturmSeq(q, Tree[0][1]);

	while (nTree > min)
	{
		const double a = Tree[min][0], b = Tree[min][1];
		const int sa = sTree[min][0], sb = sTree[min][1];
		const int s = sa - sb; // counts the number of real roots in (a, b)
		++min;

		if (s == 1) // an isolated root found
		{
			Isol[nIsol][0] = a;
			Isol[nIsol++][1] = b;
		}
		else if (s > 1) // proceed to make subdivision
		{
			const int nTree1 = nTree + 1;
			const double mid = 0.5*(a + b);
			// add intervals (a, mid] and (mid, b] to Isol
			Tree[nTree][1] = Tree[nTree1][0] = mid;
			sTree[nTree][1] = sTree[nTree1][0] = evalSturmSeq(q, mid);
			Tree[nTree][0] = a;
			Tree[nTree1][1] = b;
			sTree[nTree][0] = sa;
			sTree[nTree1][1] = sb;
			nTree += 2;
			++nIters;
			if (nIters > MAXIT)
			{ // perhaps some roots are too close

				#ifdef DEBUG_ERR
				{
					std::cout << "isolateRoots() exceeds maximum iterations\n";
				}
				#endif

				Isol[nIsol][0] = a;
				Isol[nIsol++][1] = mid;
				Isol[nIsol][0] = mid;
				Isol[nIsol++][1] = b;
				return nIsol;
			}
		} // end else if (s > 1)
	} // end while (nTree > min)
	
	return nIsol;

} // end isolateRoots()



// using Ridders' method, return the root of a polynomial p(x) known to lie between xl and x2
// this function is adopted from "Numerical recipes in C" by Press et al.
// the output is either 1, or 0 (no solution found)
bool polishRoots (const double p[DEG], const double &x1, const double &x2, double &ans)
{
	double fl = evalPoly(p, x1), fh = evalPoly(p, x2);
	if (!fh)
    {
        ans = x2;
        return 1;
    }

	if ((fl > 0)? (fh < 0) : (fh > 0))
	{
		double xl = x1, xh = x2;
		//ans = xl;
		ans = 0.5*(x1 + x2);
		for (int j = 1; j <= MAXIT; ++j)
		{
			const double xm = 0.5*(xl + xh), fm = evalPoly(p, xm), s = sqrt(fm*fm - fl*fh);
			if (!s) return 1;
			const double xnew = (fl < fh)? xm + (xl - xm)*fm/s : xm + (xm - xl)*fm/s;
			if (fabs(xnew - ans) <= ACC) return 1;
			ans = xnew;
			const double fnew = evalPoly(p, ans);
			if (!fnew) return 1;
			if (fnew >= 0? (fm<0) : (fm > 0))
			{
				xl = xm;
				fl = fm;
				xh = ans;
				fh = fnew;
			}
			else if (fnew >= 0? (fl < 0) : (fl > 0))
			{
				xh = ans;
				fh = fnew;
			}
			else
			{
				xl = ans;
				fl = fnew;
			}
			if (fabs(xh - xl) <= ACC) return 1;
		} // end for j

		#ifdef DEBUG_ERR
		{
			std::cout << "polishRoots() exceeds maximum iterations\n";
		}
		#endif

		return 0;
	} // end if
	else
	{
		#ifdef DEBUG_ERR
		{
			std::cout << "root must be bracketed in polishRoots()" << " " << x1 << " " << x2 << "\n";
		}
		#endif

		return (fabs(fl) < fabs(fh))? x1 : x2;
		//return 0;
	}

} // end polishRoots()



// find all real roots of the input square-free polynomial p(x) of degree DEG
int realRoots (const double p[DEG + 1], double roots[DEG])
{
	#ifdef DEBUG_ERR
	{
		if (!p[DEG])
			std::cout << "leading coefficient is zero in realRoots()\n";
	}
	#endif

	// copy and normalize the input polynomial p(x) and its derivative dp(x)
	double p1[DEG], dp1[DEG - 1];
	const double pdeg = 1.0/p[DEG], dpdeg = 1.0/(double)DEG;
	p1[0] = p[0]*pdeg;
	for (int i = 1; i < DEG; ++i)
	{
		p1[i] = p[i]*pdeg;
		dp1[i - 1] = (double)i*p1[i]*dpdeg;
	}	

	double Isol[DEG][2];

	// isolate all real roots of p(x)
	const int nIsol = isolateRoots(p1, dp1, Isol);

	int nr = 0;
	for (int i = 0; i < nIsol; ++i)
	{ // find an isolated real root of p(x)
		if (!polishRoots(p1, Isol[i][0], Isol[i][1], roots[nr])) continue;
		++nr;
	}

	return nr;

} // end realRoots()



#undef MAXIT
#undef MAXTREE
#undef ACC
#undef DEG

#undef DEBUG_ERR
