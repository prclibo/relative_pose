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

// void helper(double E[5][9], double C[10][20], double Z[11]);
void helper(double E[4][9], double C[9][19], double Z[11]);
// void finalize(double C[10][20], double w4, double Y[3][3]);
void finalize(double C[9][19], double w3, double Y[3][3]);

namespace pc_4pst0_nulle_poly
{

void nulle(const double x[NVIEWS][SAMPLE],
        const double y[NVIEWS][SAMPLE], const double z[NVIEWS][SAMPLE],
        double EE[4][9])
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
    std::fill(Q[4], Q[4] + 9, 0.0);
    Q[4][0] = Q[4][4] = Q[4][8] = 1.0;

    nullQR<5, 9>(Q, EE);
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

        double EE[4][9], Z[11], C[9][19];
        nulle(x, y, z, EE);

        helper(EE, C, Z);
        
        double w2[DEG], w[1][3];
        int ns = realRoots(Z, w2);
        Mat1d models, Evec(4, 9, &EE[0][0]);
        for (int i = 0; i < ns; ++i)
        {
            double Y[3][3];
            finalize(C, w2[i], Y);

            nullQR<2, 3>(Y, w);
            double w1 = w[0][2] / w[0][0], w0 = w[0][1] / w[0][0];
            Mat1d xyz1(1, 4);
            xyz1 << w0, w1, w2[i], 1;
            Mat1d E = xyz1 * Evec;
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
