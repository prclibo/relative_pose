#include <chrono>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include "relative_pose/relative_pose.hpp"
#include "precomp.hpp"
#include "relative_pose_estimator.hpp"
#include "nullqr.h"

int const DEG = 10;
int const MAXSOLS = 10;
#include "sturm.h"

int const NVIEWS = 2;
int const SAMPLE = 4;
double const IMAGE_THRESH = 1e-5;

void helper(double E[5][9], double C[10][20], double Z[11]);
void finalize(double C[10][20], double w4, double Y[3][3]);

namespace pc_4pst0_nulle_poly
{

void nulle(const double x[NVIEWS][SAMPLE],
        const double y[NVIEWS][SAMPLE], const double z[NVIEWS][SAMPLE],
        double EE[5][9])
{
    double Q[5][9];
    for (int si = 0; si < SAMPLE; ++si)
    {
        Q[si][0] = x[0][si] * x[1][si];
        Q[si][1] = y[0][si] * x[1][si];
        Q[si][2] = z[0][si] * x[1][si];
        Q[si][3] = x[0][si] * y[1][si];
        Q[si][4] = y[0][si] * y[1][si];
        Q[si][5] = z[0][si] * y[1][si];
        Q[si][6] = x[0][si] * z[1][si];
        Q[si][7] = y[0][si] * z[1][si];
        Q[si][8] = z[0][si] * z[1][si];
    }
    // std::fill(Q[4], Q[4] + 9, 0.0);
    // Q[4][0] = Q[4][4] = Q[4][8] = 1.0;

    nullQR<4, 9>(Q, EE);
}

using namespace cv;

class PC4PST0NullEEstimatorCallback CV_FINAL : public RelativePoseEstimatorCallback
{
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

        // std::cerr << "Q" << std::endl;
        // for (int si = 0; si < SAMPLE; ++si)
        //     std::cerr << x[0][si] << " " << y[0][si] << " " << z[0][si] << " " << std::endl;

        // for (int si = 0; si < SAMPLE; ++si)
        //     std::cerr << x[1][si] << " " << y[1][si] << " " << z[1][si] << " " << std::endl;

        double EE[5][9], Z[11], C[10][20];
        nulle(x, y, z, EE);

// EE[0][0] = 0.05151177341997269232276579487007;
// EE[0][1] = -0.27131267891559346372432059979474;
// EE[0][2] = -0.32000010957229924057898529099475;
// EE[0][3] = 0.20410410857739594292503682027018;
// EE[0][4] = 0.29189462387842973756235664950509;
// EE[0][5] = -0.21326368279976704034339718418778;
// EE[0][6] = 0.79809029670486375529492306668544;
// EE[0][7] = 0.03312254090616562063331684839795;
// EE[0][8] = -0.10463420059105300874424671064844;
// EE[1][0] = -0.21418643998509331871638039501704;
// EE[1][1] = 0.41456806734945200743069904092408;
// EE[1][2] = -0.09263256702017599875098596839962;
// EE[1][3] = -0.36773680117281942747453626907372;
// EE[1][4] = 0.13812275300945475731140277275699;
// EE[1][5] = 0.09910654493209901383377058436963;
// EE[1][6] = 0.13772734102294767466467817484954;
// EE[1][7] = 0.75535174808030747239939728387981;
// EE[1][8] = -0.14149647737866341556944860258227;
// EE[2][0] = -0.45913492799681659972677039149858;
// EE[2][1] = 0.18979396022946168343104034192947;
// EE[2][2] = 0.60817209477164102526813849181053;
// EE[2][3] = 0.22525796955189714143585888450616;
// EE[2][4] = -0.33831836672617110473026968975319;
// EE[2][5] = 0.12003809251331505614235339862717;
// EE[2][6] = 0.42635838175445112119277268902806;
// EE[2][7] = -0.10185425428906370870496544966954;
// EE[2][8] = -0.10739972588722650204129394069241;
// EE[3][0] = 0.68874130862296567556768422946334;
// EE[3][1] = 0.44752536773460083185938174210605;
// EE[3][2] = -0.05414946695543717408716233308041;
// EE[3][3] = -0.08307607091898437656762865799465;
// EE[3][4] = -0.34338391613088670162312610045774;
// EE[3][5] = 0.26103700914389565967965722848021;
// EE[3][6] = 0.32261208034194083227319538309530;
// EE[3][7] = -0.10528370560869046435037432729587;
// EE[3][8] = 0.11961024057334911085970219346564;
// EE[4][0] = -0.21884911246984342647614596444328;
// EE[4][1] = -0.19788010204476447206900502351345;
// EE[4][2] = 0.02454966201074772202961327138837;
// EE[4][3] = -0.39110170539660710220530859260180;
// EE[4][4] = 0.10722659666333320127584727288195;
// EE[4][5] = 0.29365321891625240091627802030416;
// EE[4][6] = 0.20230367679133648417533208885288;
// EE[4][7] = -0.08923217683296037761042640568121;
// EE[4][8] = 0.78279396157572611603114864919917;

        helper(EE, C, Z);
        double w3[DEG], w[1][3];
        int ns = realRoots(Z, w3);
        Mat1d models, Evec(5, 9, &EE[0][0]);
        Mat1d Etr = Evec.col(0) + Evec.col(4) + Evec.col(8);
        for (int i = 0; i < ns; ++i)
        {
            double Y[3][3];
            finalize(C, w3[i], Y);
            nullQR<2, 3>(Y, w);
            double w2 = w[0][2] / w[0][0], w1 = w[0][1] / w[0][0];
            double w0 = -(w1 * Etr(1, 0) + w2 * Etr(2, 0) +
                w3[i] * Etr(3, 0) + 1.0 * Etr(4, 0)) / Etr(0, 0);
            Mat1d wxyz1(1, 5);
            wxyz1 << w0, w1, w2, w3[i], 1;
            Mat1d E = wxyz1 * Evec;
            models.push_back(E.reshape(1, 3));
        }

        _model.assign(models);

        return models.rows / 3;
    }
};

} // namespace pc_4pst0_nulle_poly

namespace cv
{

Mat estimateRelativePose_PC4PST0_NullE_Poly(
        InputArray _rays1, InputArray _rays2,
        int method, double prob, double threshold, OutputArray _mask)
{
    // CV_INSTRUMENT_REGION();
    Mat rays1, rays2;
    processInputArray(_rays1, _rays2, rays1, rays2);

    Mat models;
    if( method == RANSAC )
        createRANSACPointSetRegistrator(
                makePtr<pc_4pst0_nulle_poly::PC4PST0NullEEstimatorCallback>(), 4, threshold, prob)->run(
                rays1, rays2, models, _mask);
    else
        createLMeDSPointSetRegistrator(
                makePtr<pc_4pst0_nulle_poly::PC4PST0NullEEstimatorCallback>(), 4, prob)->run(
                rays1, rays2, models, _mask);

    return models;
}

}
