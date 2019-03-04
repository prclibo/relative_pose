#include <cassert>
#include <cmath>

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
template <const int M, const int N>
static void rref(double B[M][N])
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
        assert(std::abs(B[i][i]) > 1e-12);
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

