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

int const SAMPLE = 3;
double const IMAGE_THRESH = 1e-5;

void helper(double s, double q[3][3], double qq[3][3], double C[13][25], double Z[13]);
void finalize(double C[13][25], double u3, double Y[3][3]);

namespace pc_3prast0_t2d_poly
{
using namespace cv;

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

        double s = std::cos(angle_ / 2);
        double q[3][3] = {
            {q1(0, 0)[0], q1(0, 0)[1], q1(0, 0)[2]},
            {q1(1, 0)[0], q1(1, 0)[1], q1(1, 0)[2]},
            {q1(2, 0)[0], q1(2, 0)[1], q1(2, 0)[2]}};
        double qq[3][3] = {
            {q2(0, 0)[0], q2(0, 0)[1], q2(0, 0)[2]},
            {q2(1, 0)[0], q2(1, 0)[1], q2(1, 0)[2]},
            {q2(2, 0)[0], q2(2, 0)[1], q2(2, 0)[2]}};

//         double q[3][3] = {
// {0.39244246646503067044164936305606, 0.61805234725280921992407456855290, 0.41193009485874598762933374018758},
// {0.00246488120088617090885918514687, 0.88403218237210479113485916968784, 0.88494753837647355254603098728694},
// {0.30040968936648093645658263994846, 0.58958186521471456220666595982038, 0.97842691601483089414159621810541}};
//         double qq[3][3] = {
// {0.02000510582935421943773235398112, 0.61217335215489687705314736376749, 0.67089634957624932898312408724451},
// {0.53363133850084076836850499603315, 0.96363332492597364442588059318950, 0.42750926647191278551218829306890},
// {0.13472254716347076275440031167818, 1.12194906174371134000011807074770, 0.43179828736144076906100508495001}};
// 
// 
//         s = 0.20875236849360884194837240102061;

        double C[17][25], Z[DEG + 1];
        std::fill_n(&C[0][0], 17 * 25, 0.0);
        helper(s, q, qq, C, Z);
        
        double u2[DEG], u[1][3];
        int ns = realRoots(Z, u2);
        Mat1d model;
        for (int i = 0; i < ns; ++i)
        {
            double Y[3][3];
            finalize(C, u2[i], Y);

            nullQR<2, 3>(Y, u);
            double u1 = u[0][2] / u[0][0], u0 = u[0][1] / u[0][0];

            double rbuf[9], tbuf[3];
            complementRt(u0, u1, u2[i], s, q, qq, rbuf, tbuf);
            Mat1d rmat(3, 3, rbuf), tvec(3, 1, tbuf);
            Mat1d E = skew(tvec) * rmat;

            model.push_back(E);
        }
        model = model.reshape(1);
        _model.assign(model);

        return model.rows / 3;
    }
};
} // namespace pc_3prast0_t2d_poly

namespace cv
{

Mat estimateRelativePose_PC3PRAST0_T2D_Poly(double angle,
        InputArray _rays1, InputArray _rays2,
        int method, double prob, double threshold, OutputArray _mask)
{
    // CV_INSTRUMENT_REGION();
    Mat rays1, rays2;
    processInputArray(_rays1, _rays2, rays1, rays2);

    Mat models;
    if( method == RANSAC )
        createRANSACPointSetRegistrator(
                makePtr<pc_3prast0_t2d_poly::PC3PRAST0EstimatorCallback>(angle), 3, threshold, prob)->run(
                rays1, rays2, models, _mask);
    else
        createLMeDSPointSetRegistrator(
                makePtr<pc_3prast0_t2d_poly::PC3PRAST0EstimatorCallback>(angle), 3, prob)->run(
                rays1, rays2, models, _mask);

    return models;
}
}

