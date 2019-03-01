#include <chrono>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include "relative_pose/relative_pose.hpp"
#include "precomp.hpp"
#include "relative_pose_estimator.hpp"
#include "nullqr.h"

int const DEG = 12;
int const MAXSOLS = 12;
#include "sturm.h"

int const NVIEWS = 2;
int const SAMPLE = 3;
double const IMAGE_THRESH = 1e-5;

namespace pc_3prast0_t2d
{

struct auxArrays
{
	double C[13][25], AM[12][12];
	double D[3][3][6];
	double p[13];

	void matrix13x25(const double [NVIEWS][SAMPLE], const double [NVIEWS][SAMPLE], 
		const double [NVIEWS][SAMPLE], const double &, const double &);
    void actmatrix();
	void det3x3();
	void poly12();
};



// construct 13 x 25 matrix C
void auxArrays::matrix13x25 (const double x[NVIEWS][SAMPLE], const double y[NVIEWS][SAMPLE], 
	const double z[NVIEWS][SAMPLE], const double &ss, const double &s)
{
	const int ii[3] = {0, 0, 1}, jj[3] = {1, 2, 2};
	const double t = ss - 0.5;

	double A[35], B[3][19];

	for (int m = 0; m < 3; ++m)
	{
		const int i = ii[m], j = jj[m], k = 2;

        const double u1 = -z[0][i]*y[0][j] + y[0][i]*z[0][j];
        const double u2 = z[0][i]*x[0][j] - x[0][i]*z[0][j];
        const double u3 = -y[0][i]*x[0][j] + x[0][i]*y[0][j];
        const double u4 = -z[1][i]*y[1][j] + y[1][i]*z[1][j];
        const double u5 = z[1][i]*x[1][j] - x[1][i]*z[1][j];
        const double u6 = -y[1][i]*x[1][j] + x[1][i]*y[1][j];
        const double u7 = -z[0][i]*y[0][k] + y[0][i]*z[0][k];
        const double u8 = z[0][i]*x[0][k] - x[0][i]*z[0][k];
        const double u9 = -y[0][i]*x[0][k] + x[0][i]*y[0][k];
        const double u10 = -z[1][i]*y[1][k] + y[1][i]*z[1][k];
        const double u11 = z[1][i]*x[1][k] - x[1][i]*z[1][k];
        const double u12 = -y[1][i]*x[1][k] + x[1][i]*y[1][k];

        const double S1 = x[1][j]*u1;
        const double S2 = u4*x[0][j];
        const double S3 = x[1][k]*u7;
        const double S4 = u10*x[0][k];
        const double S5 = y[1][j]*u1 + x[1][j]*u2;
        const double S6 = u5*x[0][j] + u4*y[0][j];
        const double S7 = y[1][k]*u7 + x[1][k]*u8;
        const double S8 = u11*x[0][k] + u10*y[0][k];
        const double S9 = y[1][j]*u2;
        const double S10 = u5*y[0][j];
        const double S11 = y[1][k]*u8;
        const double S12 = u11*y[0][k];
        const double S13 = z[1][j]*u1 + x[1][j]*u3;
        const double S14 = u6*x[0][j] + u4*z[0][j];
        const double S15 = z[1][k]*u7 + x[1][k]*u9;
        const double S16 = u12*x[0][k] + u10*z[0][k];
        const double S17 = z[1][j]*u2 + y[1][j]*u3;
        const double S18 = u6*y[0][j] + u5*z[0][j];
        const double S19 = z[1][k]*u8 + y[1][k]*u9;
        const double S20 = u12*y[0][k] + u11*z[0][k];
        const double S21 = z[1][j]*u3;
        const double S22 = u6*z[0][j];
        const double S23 = z[1][k]*u9;
        const double S24 = u12*z[0][k];
        const double S25 = s*(y[1][j]*u3 - z[1][j]*u2);
        const double S26 = -s*(y[0][j]*u6 - z[0][j]*u5);
        const double S27 = -s*(z[1][k]*u8 - y[1][k]*u9);
        const double S28 = s*(z[0][k]*u11 - y[0][k]*u12);
        const double S29 = -s*(x[1][j]*u3 - z[1][j]*u1);
        const double S30 = s*(x[0][j]*u6 - z[0][j]*u4);
        const double S31 = s*(z[1][k]*u7 - x[1][k]*u9);
        const double S32 = -s*(z[0][k]*u10 - x[0][k]*u12);
        const double S33 = -s*(y[1][j]*u1 - x[1][j]*u2);
        const double S34 = -s*(x[0][j]*u5 - y[0][j]*u4);
        const double S35 = -s*(y[1][k]*u7 - x[1][k]*u8);
        const double S36 = -s*(x[0][k]*u11 - y[0][k]*u10);
        const double S37 = t*(y[1][j]*u2 + x[1][j]*u1 + z[1][j]*u3);
        const double S38 = t*(x[0][j]*u4 + y[0][j]*u5 + z[0][j]*u6);
        const double S39 = t*(z[1][k]*u9 + y[1][k]*u8 + x[1][k]*u7);
        const double S40 = t*(z[0][k]*u12 + x[0][k]*u10 + y[0][k]*u11);

		if (!m)
		{
            A[0] = S1*S4 - S3*S2;
            A[1] = S1*S8 - S7*S2 - S6*S3 + S4*S5;
            A[2] = S1*S12 - S3*S10 - S2*S11 + S4*S9 + S5*S8 - S7*S6;
            A[3] = -S7*S10 - S6*S11 + S12*S5 + S8*S9;
            A[4] = -S10*S11 + S12*S9;
            A[5] = S1*S16 + S13*S4 - S14*S3 - S15*S2;
            A[6] = S1*S20 + S13*S8 - S14*S7 - S15*S6 + S16*S5 + S17*S4 - S18*S3 - S19*S2;
            A[7] = -S10*S15 - S11*S14 + S12*S13 + S16*S9 + S17*S8 - S18*S7 - S19*S6 + S20*S5;
            A[8] = -S10*S19 - S11*S18 + S12*S17 + S20*S9;
            A[9] = S1*S24 + S13*S16 - S14*S15 - S2*S23 + S21*S4 - S22*S3;
            A[10] = S13*S20 - S14*S19 - S15*S18 + S16*S17 + S21*S8 - S22*S7 - S23*S6 + S24*S5;
            A[11] = -S10*S23 - S11*S22 + S12*S21 + S17*S20 - S18*S19 + S24*S9;
            A[12] = S13*S24 - S14*S23 - S15*S22 + S16*S21;
            A[13] = S17*S24 - S18*S23 - S19*S22 + S20*S21;
            A[14] = S21*S24 - S22*S23;
            A[15] = S1*S28 - S2*S27 + S25*S4 - S26*S3;
            A[16] = S1*S32 - S2*S31 + S25*S8 - S26*S7 - S27*S6 + S28*S5 + S29*S4 - S3*S30;
            A[17] = -S10*S27 - S11*S26 + S12*S25 + S28*S9 + S29*S8 - S30*S7 - S31*S6 + S32*S5;
            A[18] = -S10*S31 - S11*S30 + S12*S29 + S32*S9;
            A[19] = S1*S36 + S13*S28 - S14*S27 - S15*S26 + S16*S25 - S2*S35 - S3*S34 + S33*S4;
            A[20] = S13*S32 - S14*S31 - S15*S30 + S16*S29 + S17*S28 - S18*S27 - S19*S26 + S20*S25 + S33*S8 - S34*S7 - S35*S6 + S36*S5;
            A[21] = -S10*S35 - S11*S34 + S12*S33 + S17*S32 - S18*S31 - S19*S30 + S20*S29 + S36*S9;
            A[22] = S13*S36 - S14*S35 - S15*S34 + S16*S33 + S21*S28 - S22*S27 - S23*S26 + S24*S25;
            A[23] = S17*S36 - S18*S35 - S19*S34 + S20*S33 + S21*S32 - S22*S31 - S23*S30 + S24*S29;
            A[24] = S21*S36 - S22*S35 - S23*S34 + S24*S33;
            A[25] = S1*S40 - S2*S39 + S25*S28 - S26*S27 - S3*S38 + S37*S4;
            A[26] = S25*S32 - S26*S31 - S27*S30 + S28*S29 + S37*S8 - S38*S7 - S39*S6 + S40*S5;
            A[27] = -S10*S39 - S11*S38 + S12*S37 + S29*S32 - S30*S31 + S40*S9;
            A[28] = S13*S40 - S14*S39 - S15*S38 + S16*S37 + S25*S36 - S26*S35 - S27*S34 + S28*S33;
            A[29] = S17*S40 - S18*S39 - S19*S38 + S20*S37 + S29*S36 - S30*S35 - S31*S34 + S32*S33;
            A[30] = S21*S40 - S22*S39 - S23*S38 + S24*S37 + S33*S36 - S34*S35;
            A[31] = S25*S40 - S26*S39 - S27*S38 + S28*S37;
            A[32] = S29*S40 - S30*S39 - S31*S38 + S32*S37;
            A[33] = S33*S40 - S34*S39 - S35*S38 + S36*S37;
            A[34] = S37*S40 - S38*S39;
		}
        B[m][0] = -x[0][i]*S2 - x[1][i]*S1;
        B[m][1] = -x[0][i]*S6 - y[0][i]*S2 - x[1][i]*S5 - y[1][i]*S1;
        B[m][2] = -y[0][i]*S6 - x[1][i]*S9 - y[1][i]*S5 - x[0][i]*S10;
        B[m][3] = -y[0][i]*S10 - y[1][i]*S9;
        B[m][4] = -z[0][i]*S2 - x[1][i]*S13 - z[1][i]*S1 - x[0][i]*S14;
        B[m][5] = -x[1][i]*S17 - y[1][i]*S13 - z[1][i]*S5 - x[0][i]*S18 - y[0][i]*S14 - z[0][i]*S6;
        B[m][6] = -y[0][i]*S18 - z[0][i]*S10 - y[1][i]*S17 - z[1][i]*S9;
        B[m][7] = -z[0][i]*S14 - x[1][i]*S21 - z[1][i]*S13 - x[0][i]*S22;
        B[m][8] = -y[0][i]*S22 - z[0][i]*S18 - y[1][i]*S21 - z[1][i]*S17;
        B[m][9] = -z[1][i]*S21 - z[0][i]*S22;
        B[m][10] = -x[1][i]*S25 - x[0][i]*S26;
        B[m][11] = -x[1][i]*S29 - y[1][i]*S25 - x[0][i]*S30 - y[0][i]*S26;
        B[m][12] = -y[1][i]*S29 - y[0][i]*S30;
        B[m][13] = -x[1][i]*S33 - z[1][i]*S25 - x[0][i]*S34 - z[0][i]*S26;
        B[m][14] = -y[1][i]*S33 - z[1][i]*S29 - y[0][i]*S34 - z[0][i]*S30;
        B[m][15] = -z[1][i]*S33 - z[0][i]*S34;
        B[m][16] = -x[1][i]*S37 - x[0][i]*S38;
        B[m][17] = -y[1][i]*S37 - y[0][i]*S38;
        B[m][18] = -z[1][i]*S37 - z[0][i]*S38;


	}

	const double r = ss - 1;
	const double rr = pow(r, 2);
    C[0][0] = A[3] - A[1];
    C[0][1] = A[0] - A[2] + A[4];
    C[0][2] = A[7] - A[5];
    C[0][3] = A[8] - A[6];
    C[0][4] = A[10] - A[1];
    C[0][5] = 2*A[0] - A[2] - A[9] + A[11];
    C[0][6] = A[17] - A[15];
    C[0][7] = A[18] - A[16];
    C[0][8] = A[12] - A[5];
    C[0][9] = A[13] - A[6];
    C[0][10] = A[0] - A[9] + A[14];
    C[0][11] = A[20];
    C[0][12] = A[21] - A[19];
    C[0][13] = A[22] - A[15];
    C[0][14] = A[23] - A[16];
    C[0][15] = A[24] - A[19];
    C[0][16] = A[26] - A[1]*r;
    C[0][17] = A[27] - A[25] + 2*A[0]*r - A[2]*r;
    C[0][18] = A[28] - A[5]*r;
    C[0][19] = A[29] - A[6]*r;
    C[0][20] = A[30] - A[25] + 2*A[0]*r - A[9]*r;
    C[0][21] = A[31] - A[15]*r;
    C[0][22] = A[32] - A[16]*r;
    C[0][23] = A[33] - A[19]*r;
    C[0][24] = A[34] - A[25]*r + A[0]*rr;
    C[1][0] = B[0][3] - B[0][1];
    C[1][1] = B[0][0] - B[0][2];
    C[1][2] = B[0][6] - B[0][4];
    C[1][3] = -B[0][5];
    C[1][4] = B[0][8] - B[0][1];
    C[1][5] = 2*B[0][0] - B[0][2] - B[0][7];
    C[1][6] = B[0][12] - B[0][10];
    C[1][7] = -B[0][11];
    C[1][8] = B[0][9] - B[0][4];
    C[1][9] = -B[0][5];
    C[1][10] = B[0][0] - B[0][7];
    C[1][11] = B[0][14];
    C[1][12] = -B[0][13];
    C[1][13] = B[0][15] - B[0][10];
    C[1][14] = -B[0][11];
    C[1][15] = -B[0][13];
    C[1][16] = B[0][17] - B[0][1]*r;
    C[1][17] = 2*B[0][0]*r - B[0][16] - B[0][2]*r;
    C[1][18] = B[0][18] - B[0][4]*r;
    C[1][19] = -B[0][5]*r;
    C[1][20] = 2*B[0][0]*r - B[0][16] - B[0][7]*r;
    C[1][21] = -B[0][10]*r;
    C[1][22] = -B[0][11]*r;
    C[1][23] = -B[0][13]*r;
    C[1][24] = -r*(B[0][16] - B[0][0]*r);
    C[2][0] = B[1][3] - B[1][1];
    C[2][1] = B[1][0] - B[1][2];
    C[2][2] = B[1][6] - B[1][4];
    C[2][3] = -B[1][5];
    C[2][4] = B[1][8] - B[1][1];
    C[2][5] = 2*B[1][0] - B[1][2] - B[1][7];
    C[2][6] = B[1][12] - B[1][10];
    C[2][7] = -B[1][11];
    C[2][8] = B[1][9] - B[1][4];
    C[2][9] = -B[1][5];
    C[2][10] = B[1][0] - B[1][7];
    C[2][11] = B[1][14];
    C[2][12] = -B[1][13];
    C[2][13] = B[1][15] - B[1][10];
    C[2][14] = -B[1][11];
    C[2][15] = -B[1][13];
    C[2][16] = B[1][17] - B[1][1]*r;
    C[2][17] = 2*B[1][0]*r - B[1][16] - B[1][2]*r;
    C[2][18] = B[1][18] - B[1][4]*r;
    C[2][19] = -B[1][5]*r;
    C[2][20] = 2*B[1][0]*r - B[1][16] - B[1][7]*r;
    C[2][21] = -B[1][10]*r;
    C[2][22] = -B[1][11]*r;
    C[2][23] = -B[1][13]*r;
    C[2][24] = -r*(B[1][16] - B[1][0]*r);
    C[3][0] = B[2][3] - B[2][1];
    C[3][1] = B[2][0] - B[2][2];
    C[3][2] = B[2][6] - B[2][4];
    C[3][3] = -B[2][5];
    C[3][4] = B[2][8] - B[2][1];
    C[3][5] = 2*B[2][0] - B[2][2] - B[2][7];
    C[3][6] = B[2][12] - B[2][10];
    C[3][7] = -B[2][11];
    C[3][8] = B[2][9] - B[2][4];
    C[3][9] = -B[2][5];
    C[3][10] = B[2][0] - B[2][7];
    C[3][11] = B[2][14];
    C[3][12] = -B[2][13];
    C[3][13] = B[2][15] - B[2][10];
    C[3][14] = -B[2][11];
    C[3][15] = -B[2][13];
    C[3][16] = B[2][17] - B[2][1]*r;
    C[3][17] = 2*B[2][0]*r - B[2][16] - B[2][2]*r;
    C[3][18] = B[2][18] - B[2][4]*r;
    C[3][19] = -B[2][5]*r;
    C[3][20] = 2*B[2][0]*r - B[2][16] - B[2][7]*r;
    C[3][21] = -B[2][10]*r;
    C[3][22] = -B[2][11]*r;
    C[3][23] = -B[2][13]*r;
    C[3][24] = -r*(B[2][16] - B[2][0]*r);
    C[4][0] = B[0][2] - B[0][0];
    C[4][1] = B[0][3] - B[0][1];
    C[4][2] = B[0][5];
    C[4][3] = B[0][6] - B[0][4];
    C[4][4] = B[0][7] - B[0][0];
    C[4][5] = B[0][8] - B[0][1];
    C[4][6] = B[0][11];
    C[4][7] = B[0][12] - B[0][10];
    C[4][8] = 0;
    C[4][9] = B[0][9] - B[0][4];
    C[4][10] = 0;
    C[4][11] = B[0][13];
    C[4][12] = B[0][14];
    C[4][13] = 0;
    C[4][14] = B[0][15] - B[0][10];
    C[4][15] = 0;
    C[4][16] = B[0][16] - B[0][0]*r;
    C[4][17] = B[0][17] - B[0][1]*r;
    C[4][18] = 0;
    C[4][19] = B[0][18] - B[0][4]*r;
    C[4][20] = 0;
    C[4][21] = 0;
    C[4][22] = -B[0][10]*r;
    C[4][23] = 0;
    C[4][24] = 0;
    C[5][0] = B[1][2] - B[1][0];
    C[5][1] = B[1][3] - B[1][1];
    C[5][2] = B[1][5];
    C[5][3] = B[1][6] - B[1][4];
    C[5][4] = B[1][7] - B[1][0];
    C[5][5] = B[1][8] - B[1][1];
    C[5][6] = B[1][11];
    C[5][7] = B[1][12] - B[1][10];
    C[5][8] = 0;
    C[5][9] = B[1][9] - B[1][4];
    C[5][10] = 0;
    C[5][11] = B[1][13];
    C[5][12] = B[1][14];
    C[5][13] = 0;
    C[5][14] = B[1][15] - B[1][10];
    C[5][15] = 0;
    C[5][16] = B[1][16] - B[1][0]*r;
    C[5][17] = B[1][17] - B[1][1]*r;
    C[5][18] = 0;
    C[5][19] = B[1][18] - B[1][4]*r;
    C[5][20] = 0;
    C[5][21] = 0;
    C[5][22] = -B[1][10]*r;
    C[5][23] = 0;
    C[5][24] = 0;
    C[6][0] = B[2][2] - B[2][0];
    C[6][1] = B[2][3] - B[2][1];
    C[6][2] = B[2][5];
    C[6][3] = B[2][6] - B[2][4];
    C[6][4] = B[2][7] - B[2][0];
    C[6][5] = B[2][8] - B[2][1];
    C[6][6] = B[2][11];
    C[6][7] = B[2][12] - B[2][10];
    C[6][8] = 0;
    C[6][9] = B[2][9] - B[2][4];
    C[6][10] = 0;
    C[6][11] = B[2][13];
    C[6][12] = B[2][14];
    C[6][13] = 0;
    C[6][14] = B[2][15] - B[2][10];
    C[6][15] = 0;
    C[6][16] = B[2][16] - B[2][0]*r;
    C[6][17] = B[2][17] - B[2][1]*r;
    C[6][18] = 0;
    C[6][19] = B[2][18] - B[2][4]*r;
    C[6][20] = 0;
    C[6][21] = 0;
    C[6][22] = -B[2][10]*r;
    C[6][23] = 0;
    C[6][24] = 0;
    C[7][0] = 0;
    C[7][1] = 0;
    C[7][2] = B[0][2] - B[0][0];
    C[7][3] = B[0][3] - B[0][1];
    C[7][4] = B[0][5];
    C[7][5] = B[0][6] - B[0][4];
    C[7][6] = 0;
    C[7][7] = 0;
    C[7][8] = B[0][7] - B[0][0];
    C[7][9] = B[0][8] - B[0][1];
    C[7][10] = B[0][9] - B[0][4];
    C[7][11] = B[0][11];
    C[7][12] = B[0][12] - B[0][10];
    C[7][13] = B[0][13];
    C[7][14] = B[0][14];
    C[7][15] = B[0][15] - B[0][10];
    C[7][16] = 0;
    C[7][17] = 0;
    C[7][18] = B[0][16] - B[0][0]*r;
    C[7][19] = B[0][17] - B[0][1]*r;
    C[7][20] = B[0][18] - B[0][4]*r;
    C[7][21] = 0;
    C[7][22] = 0;
    C[7][23] = -B[0][10]*r;
    C[7][24] = 0;
    C[8][0] = 0;
    C[8][1] = 0;
    C[8][2] = B[1][2] - B[1][0];
    C[8][3] = B[1][3] - B[1][1];
    C[8][4] = B[1][5];
    C[8][5] = B[1][6] - B[1][4];
    C[8][6] = 0;
    C[8][7] = 0;
    C[8][8] = B[1][7] - B[1][0];
    C[8][9] = B[1][8] - B[1][1];
    C[8][10] = B[1][9] - B[1][4];
    C[8][11] = B[1][11];
    C[8][12] = B[1][12] - B[1][10];
    C[8][13] = B[1][13];
    C[8][14] = B[1][14];
    C[8][15] = B[1][15] - B[1][10];
    C[8][16] = 0;
    C[8][17] = 0;
    C[8][18] = B[1][16] - B[1][0]*r;
    C[8][19] = B[1][17] - B[1][1]*r;
    C[8][20] = B[1][18] - B[1][4]*r;
    C[8][21] = 0;
    C[8][22] = 0;
    C[8][23] = -B[1][10]*r;
    C[8][24] = 0;
    C[9][0] = 0;
    C[9][1] = 0;
    C[9][2] = B[2][2] - B[2][0];
    C[9][3] = B[2][3] - B[2][1];
    C[9][4] = B[2][5];
    C[9][5] = B[2][6] - B[2][4];
    C[9][6] = 0;
    C[9][7] = 0;
    C[9][8] = B[2][7] - B[2][0];
    C[9][9] = B[2][8] - B[2][1];
    C[9][10] = B[2][9] - B[2][4];
    C[9][11] = B[2][11];
    C[9][12] = B[2][12] - B[2][10];
    C[9][13] = B[2][13];
    C[9][14] = B[2][14];
    C[9][15] = B[2][15] - B[2][10];
    C[9][16] = 0;
    C[9][17] = 0;
    C[9][18] = B[2][16] - B[2][0]*r;
    C[9][19] = B[2][17] - B[2][1]*r;
    C[9][20] = B[2][18] - B[2][4]*r;
    C[9][21] = 0;
    C[9][22] = 0;
    C[9][23] = -B[2][10]*r;
    C[9][24] = 0;
    C[10][0] = 0;
    C[10][1] = 0;
    C[10][2] = 0;
    C[10][3] = 0;
    C[10][4] = 0;
    C[10][5] = 0;
    C[10][6] = B[0][2] - B[0][0];
    C[10][7] = B[0][3] - B[0][1];
    C[10][8] = 0;
    C[10][9] = 0;
    C[10][10] = 0;
    C[10][11] = B[0][5];
    C[10][12] = B[0][6] - B[0][4];
    C[10][13] = B[0][7] - B[0][0];
    C[10][14] = B[0][8] - B[0][1];
    C[10][15] = B[0][9] - B[0][4];
    C[10][16] = B[0][11];
    C[10][17] = B[0][12] - B[0][10];
    C[10][18] = B[0][13];
    C[10][19] = B[0][14];
    C[10][20] = B[0][15] - B[0][10];
    C[10][21] = B[0][16] - B[0][0]*r;
    C[10][22] = B[0][17] - B[0][1]*r;
    C[10][23] = B[0][18] - B[0][4]*r;
    C[10][24] = -B[0][10]*r;
    C[11][0] = 0;
    C[11][1] = 0;
    C[11][2] = 0;
    C[11][3] = 0;
    C[11][4] = 0;
    C[11][5] = 0;
    C[11][6] = B[1][2] - B[1][0];
    C[11][7] = B[1][3] - B[1][1];
    C[11][8] = 0;
    C[11][9] = 0;
    C[11][10] = 0;
    C[11][11] = B[1][5];
    C[11][12] = B[1][6] - B[1][4];
    C[11][13] = B[1][7] - B[1][0];
    C[11][14] = B[1][8] - B[1][1];
    C[11][15] = B[1][9] - B[1][4];
    C[11][16] = B[1][11];
    C[11][17] = B[1][12] - B[1][10];
    C[11][18] = B[1][13];
    C[11][19] = B[1][14];
    C[11][20] = B[1][15] - B[1][10];
    C[11][21] = B[1][16] - B[1][0]*r;
    C[11][22] = B[1][17] - B[1][1]*r;
    C[11][23] = B[1][18] - B[1][4]*r;
    C[11][24] = -B[1][10]*r;
    C[12][0] = 0;
    C[12][1] = 0;
    C[12][2] = 0;
    C[12][3] = 0;
    C[12][4] = 0;
    C[12][5] = 0;
    C[12][6] = B[2][2] - B[2][0];
    C[12][7] = B[2][3] - B[2][1];
    C[12][8] = 0;
    C[12][9] = 0;
    C[12][10] = 0;
    C[12][11] = B[2][5];
    C[12][12] = B[2][6] - B[2][4];
    C[12][13] = B[2][7] - B[2][0];
    C[12][14] = B[2][8] - B[2][1];
    C[12][15] = B[2][9] - B[2][4];
    C[12][16] = B[2][11];
    C[12][17] = B[2][12] - B[2][10];
    C[12][18] = B[2][13];
    C[12][19] = B[2][14];
    C[12][20] = B[2][15] - B[2][10];
    C[12][21] = B[2][16] - B[2][0]*r;
    C[12][22] = B[2][17] - B[2][1]*r;
    C[12][23] = B[2][18] - B[2][4]*r;
    C[12][24] = -B[2][10]*r;

} // end auxArrays::matrix13x25()

void auxArrays::actmatrix()
{

    rref<13, 25>(C);

    for (int i = 0; i < 5; ++i)
    {
        std::copy(C[i + 8] + 13, C[i + 8] + 25, AM[i]);
        for(int j = 0; j < 12; ++j) AM[i][j] *= -1;
    }
    for (int i = 5; i < 12; ++i)
        std::fill(AM[i], AM[i] + 12, 0.0);

    AM[5][0] = AM[6][1] = AM[7][2] = AM[8][5] = AM[9][6] = AM[10][7] = AM[11][10] = 1.0;
}

// compute determinant of 3x3 matrix, the result is a polynomial p(x) of degree 12.
void auxArrays::det3x3 ()
{
	double q[3][10];
	const int i2[3][2] = {{0, 1}, {0, 2}, {1, 2}};

	for (int i = 0; i < 3; ++i)
	{ // 3 2-minors
		const int j = i2[i][0], k = i2[i][1];
		q[i][0] = D[j][1][0]*D[k][2][0]-D[j][2][0]*D[k][1][0];
		q[i][1] = D[j][1][0]*D[k][2][1]+D[j][1][1]*D[k][2][0]-D[j][2][0]*D[k][1][1]-D[j][2][1]*D[k][1][0];
		q[i][2] = D[j][1][0]*D[k][2][2]+D[j][1][1]*D[k][2][1]+D[j][1][2]*D[k][2][0]-D[j][2][0]*D[k][1][2]-D[j][2][1]*D[k][1][1]-D[j][2][2]*D[k][1][0];
		q[i][3] = D[j][1][0]*D[k][2][3]+D[j][1][1]*D[k][2][2]+D[j][1][2]*D[k][2][1]+D[j][1][3]*D[k][2][0]-D[j][2][0]*D[k][1][3]-D[j][2][1]*D[k][1][2]-D[j][2][2]*D[k][1][1]-D[j][2][3]*D[k][1][0];
		q[i][4] = D[j][1][0]*D[k][2][4]+D[j][1][1]*D[k][2][3]+D[j][1][2]*D[k][2][2]+D[j][1][3]*D[k][2][1]+D[j][1][4]*D[k][2][0]-D[j][2][0]*D[k][1][4]-D[j][2][1]*D[k][1][3]-D[j][2][2]*D[k][1][2]-D[j][2][3]*D[k][1][1]-D[j][2][4]*D[k][1][0];
		q[i][5] = D[j][1][0]*D[k][2][5]+D[j][1][1]*D[k][2][4]+D[j][1][2]*D[k][2][3]+D[j][1][3]*D[k][2][2]+D[j][1][4]*D[k][2][1]-D[j][2][1]*D[k][1][4]-D[j][2][2]*D[k][1][3]-D[j][2][3]*D[k][1][2]-D[j][2][4]*D[k][1][1]-D[j][2][5]*D[k][1][0];
		q[i][6] = D[j][1][1]*D[k][2][5]+D[j][1][2]*D[k][2][4]+D[j][1][3]*D[k][2][3]+D[j][1][4]*D[k][2][2]-D[j][2][2]*D[k][1][4]-D[j][2][3]*D[k][1][3]-D[j][2][4]*D[k][1][2]-D[j][2][5]*D[k][1][1];
		q[i][7] = D[j][1][2]*D[k][2][5]+D[j][1][3]*D[k][2][4]+D[j][1][4]*D[k][2][3]-D[j][2][3]*D[k][1][4]-D[j][2][4]*D[k][1][3]-D[j][2][5]*D[k][1][2];
		q[i][8] = D[j][1][3]*D[k][2][5]+D[j][1][4]*D[k][2][4]-D[j][2][4]*D[k][1][4]-D[j][2][5]*D[k][1][3];
		q[i][9] = D[j][1][4]*D[k][2][5]-D[j][2][5]*D[k][1][4];
	}

	p[0] = q[0][0]*D[2][0][0]-q[1][0]*D[1][0][0]+q[2][0]*D[0][0][0];
	p[1] = q[0][0]*D[2][0][1]+q[0][1]*D[2][0][0]-q[1][0]*D[1][0][1]-q[1][1]*D[1][0][0]+q[2][0]*D[0][0][1]+q[2][1]*D[0][0][0];
	p[2] = q[0][0]*D[2][0][2]+q[0][1]*D[2][0][1]+q[0][2]*D[2][0][0]-q[1][0]*D[1][0][2]-q[1][1]*D[1][0][1]-q[1][2]*D[1][0][0]+q[2][0]*D[0][0][2]+q[2][1]*D[0][0][1]+q[2][2]*D[0][0][0];
	p[3] = q[0][0]*D[2][0][3]+q[0][1]*D[2][0][2]+q[0][2]*D[2][0][1]+q[0][3]*D[2][0][0]-q[1][0]*D[1][0][3]-q[1][1]*D[1][0][2]-q[1][2]*D[1][0][1]-q[1][3]*D[1][0][0]+q[2][0]*D[0][0][3]+q[2][1]*D[0][0][2]+q[2][2]*D[0][0][1]+q[2][3]*D[0][0][0];
	p[4] = q[0][1]*D[2][0][3]+q[0][2]*D[2][0][2]+q[0][3]*D[2][0][1]+q[0][4]*D[2][0][0]-q[1][1]*D[1][0][3]-q[1][2]*D[1][0][2]-q[1][3]*D[1][0][1]-q[1][4]*D[1][0][0]+q[2][1]*D[0][0][3]+q[2][2]*D[0][0][2]+q[2][3]*D[0][0][1]+q[2][4]*D[0][0][0];
	p[5] = q[0][2]*D[2][0][3]+q[0][3]*D[2][0][2]+q[0][4]*D[2][0][1]+q[0][5]*D[2][0][0]-q[1][2]*D[1][0][3]-q[1][3]*D[1][0][2]-q[1][4]*D[1][0][1]-q[1][5]*D[1][0][0]+q[2][2]*D[0][0][3]+q[2][3]*D[0][0][2]+q[2][4]*D[0][0][1]+q[2][5]*D[0][0][0];
	p[6] = q[0][3]*D[2][0][3]+q[0][4]*D[2][0][2]+q[0][5]*D[2][0][1]+q[0][6]*D[2][0][0]-q[1][3]*D[1][0][3]-q[1][4]*D[1][0][2]-q[1][5]*D[1][0][1]-q[1][6]*D[1][0][0]+q[2][3]*D[0][0][3]+q[2][4]*D[0][0][2]+q[2][5]*D[0][0][1]+q[2][6]*D[0][0][0];
	p[7] = q[0][4]*D[2][0][3]+q[0][5]*D[2][0][2]+q[0][6]*D[2][0][1]+q[0][7]*D[2][0][0]-q[1][4]*D[1][0][3]-q[1][5]*D[1][0][2]-q[1][6]*D[1][0][1]-q[1][7]*D[1][0][0]+q[2][4]*D[0][0][3]+q[2][5]*D[0][0][2]+q[2][6]*D[0][0][1]+q[2][7]*D[0][0][0];
	p[8] = q[0][5]*D[2][0][3]+q[0][6]*D[2][0][2]+q[0][7]*D[2][0][1]+q[0][8]*D[2][0][0]-q[1][5]*D[1][0][3]-q[1][6]*D[1][0][2]-q[1][7]*D[1][0][1]-q[1][8]*D[1][0][0]+q[2][5]*D[0][0][3]+q[2][6]*D[0][0][2]+q[2][7]*D[0][0][1]+q[2][8]*D[0][0][0];
	p[9] = q[0][6]*D[2][0][3]+q[0][7]*D[2][0][2]+q[0][8]*D[2][0][1]+q[0][9]*D[2][0][0]-q[1][6]*D[1][0][3]-q[1][7]*D[1][0][2]-q[1][8]*D[1][0][1]-q[1][9]*D[1][0][0]+q[2][6]*D[0][0][3]+q[2][7]*D[0][0][2]+q[2][8]*D[0][0][1]+q[2][9]*D[0][0][0];
	p[10] = q[0][7]*D[2][0][3]+q[0][8]*D[2][0][2]+q[0][9]*D[2][0][1]-q[1][7]*D[1][0][3]-q[1][8]*D[1][0][2]-q[1][9]*D[1][0][1]+q[2][7]*D[0][0][3]+q[2][8]*D[0][0][2]+q[2][9]*D[0][0][1];
	p[11] = q[0][8]*D[2][0][3]+q[0][9]*D[2][0][2]-q[1][8]*D[1][0][3]-q[1][9]*D[1][0][2]+q[2][8]*D[0][0][3]+q[2][9]*D[0][0][2];
	p[12] = q[0][9]*D[2][0][3]-q[1][9]*D[1][0][3]+q[2][9]*D[0][0][3];

} // end auxArrays::det3x3()



// compute 12th degree univariate polynomial
void auxArrays::poly12 ()
{
	gj<13, 25>(C);

	// w^0
	D[0][0][0] = -C[8][15];
	D[0][1][0] = -C[8][19];
	D[0][2][0] = -C[8][24];
	D[1][0][0] = -C[10][15];
	D[1][1][0] = -C[10][19];
	D[1][2][0] = -C[10][24];
	D[2][0][0] = -C[11][15];
	D[2][1][0] = -C[11][19];
	D[2][2][0] = -C[11][24];
	// w^1
	D[0][0][1] = -C[8][14] + C[9][15];
	D[0][1][1] = -C[8][18] + C[9][19];
	D[0][2][1] = -C[8][23] + C[9][24];
	D[1][0][1] = -C[10][14] + C[11][15];
	D[1][1][1] = -C[10][18] + C[11][19];
	D[1][2][1] = -C[10][23] + C[11][24];
	D[2][0][1] = -C[11][14] + C[12][15];
	D[2][1][1] = -C[11][18] + C[12][19];
	D[2][2][1] = -C[11][23] + C[12][24];
	// w^2
	D[0][0][2] = -C[8][13] + C[9][14];
	D[0][1][2] = -C[8][17] + C[9][18];
	D[0][2][2] = -C[8][22] + C[9][23];
	D[1][0][2] = -C[10][13] + C[11][14];
	D[1][1][2] = -C[10][17] + C[11][18];
	D[1][2][2] = -C[10][22] + C[11][23];
	D[2][0][2] = -C[11][13] + C[12][14];
	D[2][1][2] = -C[11][17] + C[12][18];
	D[2][2][2] = -C[11][22] + C[12][23];
	// w^3
	D[0][0][3] = C[9][13];
	D[0][1][3] = -C[8][16] + C[9][17];
	D[0][2][3] = -C[8][21] + C[9][22];
	D[1][0][3] = C[11][13];
	D[1][1][3] = -C[10][16] + C[11][17];
	D[1][2][3] = -C[10][21] + C[11][22];
	D[2][0][3] = C[12][13];
	D[2][1][3] = -C[11][16] + C[12][17];
	D[2][2][3] = -C[11][21] + C[12][22];
	// w^4
	D[0][0][4] = 0;
	D[0][1][4] = C[9][16];
	D[0][2][4] = -C[8][20] + C[9][21];
	D[1][0][4] = 0;
	D[1][1][4] = C[11][16];
	D[1][2][4] = -C[10][20] + C[11][21];
	D[2][0][4] = 0;
	D[2][1][4] = C[12][16];
	D[2][2][4] = -C[11][20] + C[12][21];
	// w^5
	D[0][0][5] = 0;
	D[0][1][5] = 0;
	D[0][2][5] = C[9][20];
	D[1][0][5] = 0;
	D[1][1][5] = 0;
	D[1][2][5] = C[11][20];
	D[2][0][5] = 0;
	D[2][1][5] = 0;
	D[2][2][5] = C[12][20];

	det3x3();

} // end auxArrays::poly12()

// main function
// the output is the number of real solutions found
int pc3prast0 (const double x[NVIEWS][SAMPLE], const double y[NVIEWS][SAMPLE], 
	const double z[NVIEWS][SAMPLE], const double &ss,
    double rbuf[MAXSOLS * 9], double tbuf[MAXSOLS * 3])
{
	auxArrays S;
	const double s = sqrt(ss), b = sqrt(1 - ss);
	S.matrix13x25(x, y, z, ss, s);
    S.actmatrix(); exit(0);
	S.poly12();

	double w1[MAXSOLS];
	const int ns = realRoots(S.p, w1);

	for (int k = 0; k < ns; ++k)
	{ // cycle over all solutions

		const int dd[3] = {3, 4, 5};
		double w = w1[k], D1[2][3];
		for (int j = 0; j < 3; ++j)
		{
			const int d1 = dd[j];
			for (int i = 0; i < 2; ++i)
			{
				D1[i][j] = S.D[i][j][d1];
				for (int d = d1 - 1; d >= 0; --d)
					D1[i][j] = D1[i][j]*w + S.D[i][j][d];
			}
		}

		double q[1][3];
		nullQR<2, 3>(D1, q);
		double den = 1.0/q[0][2], u = q[0][0]*den, v = q[0][1]*den;
		const double fac = b/sqrt(pow(u, 2) + pow(v, 2) + pow(w, 2));
		u *= fac;
		v *= fac;
		w *= fac;

		// compute rotation matrix
		const double u2 = 2*u, v2 = 2*v, w2 = 2*w;
		const double uu2 = u*u2, vv2 = v*v2, ww2 = w*w2, ss2 = 2*ss;
		const double uv2 = u2*v, vw2 = v2*w, uw2 = u2*w, us2 = u2*s, vs2 = v2*s, ws2 = w2*s;

        double* const rmat = rbuf + 9 * k;
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
        const double t1 = rmat[6]*x[0][0] + rmat[7]*y[0][0] + rmat[8]*z[0][0];
        const double t2 = rmat[3]*x[0][0] + rmat[4]*y[0][0] + rmat[5]*z[0][0];
        const double t3 = rmat[0]*x[0][0] + rmat[1]*y[0][0] + rmat[2]*z[0][0];
        const double t4 = rmat[6]*x[0][1] + rmat[7]*y[0][1] + rmat[8]*z[0][1];
        const double t5 = rmat[3]*x[0][1] + rmat[4]*y[0][1] + rmat[5]*z[0][1];
        const double t6 = rmat[0]*x[0][1] + rmat[1]*y[0][1] + rmat[2]*z[0][1];
		S[0][0] = t2*z[1][0] - t1*y[1][0];
		S[0][1] = t1*x[1][0] - t3*z[1][0];
		S[0][2] = t3*y[1][0] - t2*x[1][0];
		S[1][0] = t5*z[1][1] - t4*y[1][1];
		S[1][1] = t4*x[1][1] - t6*z[1][1];
		S[1][2] = t6*y[1][1] - t5*x[1][1];
		nullQR<2, 3>(S, t);

		// normalize translation vectors so that ||t2|| = 1
        const double n = sqrt(pow(t[0][0], 2) + pow(t[0][1], 2) + pow(t[0][2], 2));

        tbuf[k * 3 + 0] = t[0][0] / n;
        tbuf[k * 3 + 1] = t[0][1] / n;
        tbuf[k * 3 + 2] = t[0][2] / n;


	} // end of cycle over all solutions

	return ns;

} // end pc3prast0()

} // namespace pc_3prast0_t2d

void recoverRt(double u, double v, double w, double s,
        double x[NVIEWS][SAMPLE], double y[NVIEWS][SAMPLE],
        double z[NVIEWS][SAMPLE], double rmat[9], double tvec[3])
{
    // compute rotation matrix
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
    const double t1 = rmat[6]*x[0][0] + rmat[7]*y[0][0] + rmat[8]*z[0][0];
    const double t2 = rmat[3]*x[0][0] + rmat[4]*y[0][0] + rmat[5]*z[0][0];
    const double t3 = rmat[0]*x[0][0] + rmat[1]*y[0][0] + rmat[2]*z[0][0];
    const double t4 = rmat[6]*x[0][1] + rmat[7]*y[0][1] + rmat[8]*z[0][1];
    const double t5 = rmat[3]*x[0][1] + rmat[4]*y[0][1] + rmat[5]*z[0][1];
    const double t6 = rmat[0]*x[0][1] + rmat[1]*y[0][1] + rmat[2]*z[0][1];
    S[0][0] = t2*z[1][0] - t1*y[1][0];
    S[0][1] = t1*x[1][0] - t3*z[1][0];
    S[0][2] = t3*y[1][0] - t2*x[1][0];
    S[1][0] = t5*z[1][1] - t4*y[1][1];
    S[1][1] = t4*x[1][1] - t6*z[1][1];
    S[1][2] = t6*y[1][1] - t5*x[1][1];
    nullQR<2, 3>(S, t);

    // normalize translation vectors so that ||t2|| = 1
    const double n = sqrt(pow(t[0][0], 2) + pow(t[0][1], 2) + pow(t[0][2], 2));

    tvec[0] = t[0][0] / n;
    tvec[1] = t[0][1] / n;
    tvec[2] = t[0][2] / n;
}

namespace cv
{

class PC3PRAST0EstimatorCallback CV_FINAL : public RelativePoseEstimatorCallback
{
protected:
    float angle_, dist_thresh_;

public:
    PC3PRAST0EstimatorCallback(float angle, float dist_thresh = 100)
        : angle_(angle)
        , dist_thresh_(dist_thresh) {}

    int runKernel( InputArray _m1, InputArray _m2, OutputArray _model ) const CV_OVERRIDE
    {
        Mat3d q1 = _m1.getMat(), q2 = _m2.getMat();
        CV_Assert(q1.cols == 1 && q2.cols == 1);

        double x[NVIEWS][SAMPLE], y[NVIEWS][SAMPLE], z[NVIEWS][SAMPLE];
        for (int si = 0; si < SAMPLE; ++si)
        {
            x[0][si] = q1(si, 0)[0];
            y[0][si] = q1(si, 0)[1];
            z[0][si] = q1(si, 0)[2];
            x[1][si] = q2(si, 0)[0];
            y[1][si] = q2(si, 0)[1];
            z[1][si] = q2(si, 0)[2];
        }

        const double s = std::cos(angle_ / 2);
        const double ss = s * s;

        pc_3prast0_t2d::auxArrays aux;
	    aux.matrix13x25(x, y, z, ss, s);
        aux.actmatrix();

        Eigen::Map<Eigen::MatrixXd> m(&aux.AM[0][0], 12, 12);
        m.transposeInPlace();

        Eigen::EigenSolver<Eigen::MatrixXd> eigen(m);
        Eigen::MatrixXcd V = eigen.eigenvectors().block<3, 12>(8, 0);
        Eigen::MatrixXd uvecs =
            V.real().colwise().normalized() * std::sqrt(1.0 - ss);
        Eigen::MatrixXd Via = V.imag().cwiseAbs();

        Mat1d model;
        for (int si = 0; si < V.cols(); ++si)
        {
            if (!Via.col(si).isApproxToConstant(0, IMAGE_THRESH)) continue;
            double u = uvecs(0, si), v = uvecs(1, si), w = uvecs(2, si);
            Mat1d rmat(3, 3), tvec(3, 1);
            // FIXME inconsistent rotation<->quaterion conversion.
            recoverRt(-u, -v, -w, s, x, y, z,
                    rmat.ptr<double>(), tvec.ptr<double>());
            Mat1d E = skew(tvec) * rmat;

            model.push_back(E);
        }
        model = model.reshape(1);
        _model.assign(model);

        return model.rows / 3;
    }
};

Mat estimateRelativePose_PC3PRAST0_T2D(double angle,
        InputArray _rays1, InputArray _rays2,
        int method, double prob, double threshold, OutputArray _mask)
{
    // CV_INSTRUMENT_REGION();
    Mat rays1, rays2;
    processInputArray(_rays1, _rays2, rays1, rays2);

    Mat models;
    if( method == RANSAC )
        createRANSACPointSetRegistrator(
                makePtr<PC3PRAST0EstimatorCallback>(angle), 3, threshold, prob)->run(
                rays1, rays2, models, _mask);
    else
        createLMeDSPointSetRegistrator(
                makePtr<PC3PRAST0EstimatorCallback>(angle), 3, prob)->run(
                rays1, rays2, models, _mask);

    return models;
}
}

