#include <cassert>
#include <cmath>
#include <cstring>

double const RREF_ZERO_THRESH = 1e-12;

// find the basis of the null space of an M x N matrix A (M < N) using specifically tailored QR factorization,
// result is N - M N-vectors Q
template <const int M, const int N>
static void nullQR (double A[M][N], double Q[N - M][N]) // matrix A changes during computation!!!
{
    const int NM = N - M, M1 = M - 1;

    // find Hauseholder vectors
    for (int j = 0; j < M; ++j)
    {
        double t = 0;
        for (int i = j + 1; i < N; ++i)
            t += pow(A[j][i], 2);
        double mu = sqrt(pow(A[j][j], 2) + t);
        mu = (A[j][j]<0)? 1.0/(A[j][j] - mu) : 1.0/(A[j][j] + mu);

        //A[j][j] = 1.0;
        for (int i = j + 1; i < N; ++i)
            A[j][i] *= mu;

        const double beta = -2.0/(1.0 + t*pow(mu, 2));
        for (int k = j + 1; k < M; ++k)
        {
            double w = A[k][j];
            for (int i = j + 1; i < N; ++i)
                w += A[k][i]*A[j][i];
            w *= beta;
            for (int i = j + 1; i < N; ++i)
                A[k][i] += A[j][i]*w;
        }
    }

    // multiply m Householder matrices, inverse order is more efficient
    // we only need last N - M columns of the resulting matrix
    for (int i = 0; i < NM; ++i)
        memset(Q[i], 0, sizeof(double)*M1);
    //    for (int k = 0; k < M1; ++k)
    //        Q[i][k] = 0.0;
    
    // start from Mth matrix
    double beta = 1.0;
    for (int i = M; i < N; ++i)
        beta += pow(A[M1][i], 2);

    beta = -2.0/beta;
    for (int k = 0; k < NM; ++k)
    {
        Q[k][M1] = A[M1][k + M]*beta;
        for (int i = M; i < N; ++i)
            Q[k][i] = (k == i - M)? 1.0 + A[M1][i]*Q[k][M1] : A[M1][i]*Q[k][M1];
    }
    
    // multiply by the remaining M-1 matrices
    for (int j = M1 - 1; j >= 0; --j)
    {
        double beta = 1.0;
        for (int i = j + 1; i < N; ++i)
            beta += pow(A[j][i], 2);

        beta = -2.0/beta;
        for (int k = 0; k < NM; ++k)
        {
            double w = Q[k][j];
            for (int i = j + 1; i < N; ++i)
                w += Q[k][i]*A[j][i];
            w *= beta;
            Q[k][j] += w;
            for (int i = j + 1; i < N; ++i)
                Q[k][i] += A[j][i]*w;
        }
    }

} // end nullQR()


// Specifically tailored Gauss-Jordan elimination with partial pivoting on M x N matrix B
template <const int N, const int M>
static int rref(double B[N][M])
{
    int piv_col = 0, piv_cols[N];
    for (int i = 0; i < N; ++i)
    {
        double piv = B[i][piv_col];
        int piv_row = i;
        while (true)
        {
            for (int j = i + 1; j < N; ++j)
            {
                if (std::abs(B[j][piv_col]) > std::abs(piv))
                {
                    piv = B[j][piv_col];
                    piv_row = j;
                }
            }
            if (std::abs(piv) > RREF_ZERO_THRESH || piv_col >= M) break;
            else ++piv_col;
        }
        piv_cols[i] = piv_col;
        if (piv_col >= M) break;

        double inv_piv = 1.0 / piv;
        for (int j = 0; j < M; ++j)
        {
            std::swap(B[i][j], B[piv_row][j]);
            B[i][j] *= inv_piv;
        }

        for (int j = i + 1; j < N; ++j)
        {
            double fac = B[j][piv_col];
            for (int k = 0; k < M; ++k)
                B[j][k] -= B[i][k] * fac;
        }
        ++piv_col;
    }
    int ri = 0;
    for (; ri < N; ++ri)
    {
        if (piv_cols[ri] >= M) break;
        for (int j = 0; j < ri; ++j)
        {
            double fac = B[j][piv_cols[ri]];
            for (int k = 0; k < M; ++k)
                B[j][k] -= B[ri][k] * fac;
        }
    }
    return ri;
}

// Specifically tailored Gauss-Jordan elimination with partial pivoting on M x N matrix B
template <const int M, const int N>
static void _rref(double B[M][N])
{
    for (int i = 0; i < M; ++i)
    { // index i numbers first M columns of B
        int ipiv = i, i1 = i + 1;
        double piv = B[i][i];
        for (int k = i1; k < M; ++k)
        {
            if (fabs(B[k][i]) > fabs(piv))
            {
                piv = B[k][i];
                ipiv = k;
            }
        }

        if (!piv) continue;
        // B[ipiv][i] = B[i][i];
        // B[i][i] = piv;
        piv = 1.0/piv;
        for (int j = 0; j < N; ++j)
        { // index j numbers columns of B from i1 to N - 1
            const double t = B[ipiv][j];
            B[ipiv][j] = B[i][j];
            B[i][j] = t*piv;
        }

        for (int k = i1; k < M; ++k)
        {
            double fac = B[k][i];
            for (int j = 0; j < N; ++j)
                B[k][j] -= B[i][j] * fac;
        }

    }
    for (int i = M - 1; i >= 0; --i)
    {
        assert(std::abs(B[i][i]) > RREF_ZERO_THRESH);
        for (int k = 0; k < i; ++k)
        {
            double fac = B[k][i] / B[i][i];
            for (int j = i; j < N; ++j)
                B[k][j] -= B[i][j] * fac;
        }
    }

} // end gj()


template <const int M, const int N>
static void gj(double B[M][N])
{
    for (int i = 0; i < M; ++i)
    { // index i numbers first M columns of B
        int ipiv = i, i1 = i + 1;
        double piv = B[i][i];
        for (int k = i1; k < M; ++k)
        {
            if (fabs(B[k][i]) > fabs(piv))
            {
                piv = B[k][i];
                ipiv = k;
            }
        }

        const int k1 = (i < 10)? i1 : 10; // optimization!!!

        if (!piv) continue;
        B[ipiv][i] = B[i][i];
        B[i][i] = piv;
        piv = 1.0/piv;
        for (int j = i1; j < N; ++j)
        { // index j numbers columns of B from i1 to N - 1
            const double t = B[ipiv][j];
            B[ipiv][j] = B[i][j];
            B[i][j] = t*piv;
            for (int k = k1; k < M; ++k)
            { // index k numbers rows of B from k1 to M - 1
                if (k == i) continue;
                B[k][j] -= B[i][j]*B[k][i];
            }
        }
    }

} // end gj()

        // compute rotation matrix
static void complementRt(double u, double v, double w, double s,
        double q[][3], double qq[][3], double rmat[9], double tvec[3])
{
        const double u2 = 2*u, v2 = 2*v, w2 = 2*w, ss = s * s;
        const double uu2 = u*u2, vv2 = v*v2, ww2 = w*w2, ss2 = 2*ss;
        const double uv2 = u2*v, vw2 = v2*w, uw2 = u2*w, us2 = u2*s, vs2 = v2*s, ws2 = w2*s;

        rmat[0] = uu2 + ss2 - 1;
        rmat[1] = ws2 + uv2;
        rmat[2] = uw2 - vs2;
        rmat[3] = uv2 - ws2;
        rmat[4] = vv2 + ss2 - 1;
        rmat[5] = us2 + vw2;
        rmat[6] = vs2 + uw2;
        rmat[7] = vw2 - us2;
        rmat[8] = ww2 + ss2 - 1;

        // compute translation vector t
        double S[2][3], t[1][3];
        const double t1 = rmat[6]*q[0][0] + rmat[7]*q[0][1] + rmat[8]*q[0][2];
        const double t2 = rmat[3]*q[0][0] + rmat[4]*q[0][1] + rmat[5]*q[0][2];
        const double t3 = rmat[0]*q[0][0] + rmat[1]*q[0][1] + rmat[2]*q[0][2];
        const double t4 = rmat[6]*q[1][0] + rmat[7]*q[1][1] + rmat[8]*q[1][2];
        const double t5 = rmat[3]*q[1][0] + rmat[4]*q[1][1] + rmat[5]*q[1][2];
        const double t6 = rmat[0]*q[1][0] + rmat[1]*q[1][1] + rmat[2]*q[1][2];
        S[0][0] = t2*qq[0][2] - t1*qq[0][1];
        S[0][1] = t1*qq[0][0] - t3*qq[0][2];
        S[0][2] = t3*qq[0][1] - t2*qq[0][0];
        S[1][0] = t5*qq[1][2] - t4*qq[1][1];
        S[1][1] = t4*qq[1][0] - t6*qq[1][2];
        S[1][2] = t6*qq[1][1] - t5*qq[1][0];
        nullQR<2, 3>(S, t);

        // normalize translation vectors so that ||t2|| = 1
        const double n = sqrt(pow(t[0][0], 2) + pow(t[0][1], 2) + pow(t[0][2], 2));

        tvec[0] = t[0][0] / n;
        tvec[1] = t[0][1] / n;
        tvec[2] = t[0][2] / n;
}
