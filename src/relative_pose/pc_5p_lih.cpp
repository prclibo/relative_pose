#include <chrono>
#include <opencv2/opencv.hpp>
#include "relative_pose/relative_pose.hpp"
#include "precomp.hpp"
#include "relative_pose_estimator.hpp"

void compute_E_matrices (double q[][3], double qp[][3],
        double Ematrices[10][3][3], int &nroots, bool optimized);

namespace cv
{
class PC5PLiHEstimatorCallback CV_FINAL : public RelativePoseEstimatorCallback// RelativePoseEstimatorCallback
{
protected:
    float dist_thresh_;

public:
    PC5PLiHEstimatorCallback(float dist_thresh = 100)
        : dist_thresh_(dist_thresh) {}

    int runKernel( InputArray _m1, InputArray _m2, OutputArray _model ) const CV_OVERRIDE
    {
        Mat3d q1 = _m1.getMat(), q2 = _m2.getMat();
        CV_Assert(q1.cols == 1 && q2.cols == 1);

        double q[5][3], qp[5][3];
        for (int si = 0; si < 5; ++si)
        {
            q[si][0] = q1(si, 0)[0];
            q[si][1] = q1(si, 0)[1];
            q[si][2] = q1(si, 0)[2];
            qp[si][0] = q2(si, 0)[0];
            qp[si][1] = q2(si, 0)[1];
            qp[si][2] = q2(si, 0)[2];
        }

        double ematrices[10][3][3];
        int num_roots;
        compute_E_matrices(q, qp, ematrices, num_roots, false);

        auto models = Mat1d(num_roots * 3, 3, &ematrices[0][0][0]).clone();
        _model.assign(models);

        return models.rows / 3;
    }
};

Mat estimateRelativePose_PC5P_LiH(InputArray _rays1, InputArray _rays2,
        int method, double prob, double threshold, OutputArray _mask)
{
    Mat rays1, rays2;
    processInputArray(_rays1, _rays2, rays1, rays2);

    Mat models;
    if( method == RANSAC )
        createRANSACPointSetRegistrator(
                makePtr<PC5PLiHEstimatorCallback>(), 5, threshold, prob)->run(
                rays1, rays2, models, _mask);
    else
        createLMeDSPointSetRegistrator(
                makePtr<PC5PLiHEstimatorCallback>(), 5, prob)->run(
                rays1, rays2, models, _mask);

    return models;

}

}
