#include <chrono>
#include <opencv2/opencv.hpp>
#include "relative_pose/relative_pose.hpp"
#include "precomp.hpp"
#include "relative_pose_estimator.hpp"

void compute_E_matrices (double q[][3], double qp[][3],
        double Ematrices[10][3][3], int &nroots, bool optimized);

namespace cv
{

class PC5PLiHEstimatorCallback CV_FINAL : public RelativePoseEstimatorCallback 
{
protected:
    float dist_thresh_;

public:
    PC5PLiHEstimatorCallback(float dist_thresh = 100)
        : dist_thresh_(dist_thresh) {}

    int runKernel( InputArray _m1, InputArray _m2, OutputArray _model ) const CV_OVERRIDE
    {
        Mat2d q1 = _m1.getMat(), q2 = _m2.getMat();
        CV_Assert(q1.type() == CV_64FC2 && q2.type() == CV_64FC2);
        CV_Assert(q1.cols == 1 && q2.cols == 1);


        double q[5][3], qp[5][3];
        for (int si = 0; si < 5; ++si)
        {
            q[si][0] = q1(si, 0)[0];
            q[si][1] = q1(si, 0)[1];
            q[si][2] = 1;
            qp[si][0] = q2(si, 0)[0];
            qp[si][1] = q2(si, 0)[1];
            qp[si][2] = 1;
        }

        double ematrices[10][3][3];
        int num_roots;
        compute_E_matrices(q, qp, ematrices, num_roots, false);

        auto models = Mat1d(num_roots * 3, 3, &ematrices[0][0][0]).clone();
        _model.assign(models);

        return models.rows / 3;
    }
};

Mat estimateRelativePose_PC5P_LiH(InputArray _points1, InputArray _points2,
        InputArray _cameraMatrix, int method, double prob, double threshold,
        OutputArray _mask)
{
    Mat points1, points2, cameraMatrix;
    processInputArray(_points1, _points2, _cameraMatrix, threshold,
            points1, points2, cameraMatrix, threshold);

    Mat models;
    if( method == RANSAC )
        createRANSACPointSetRegistrator(
                makePtr<PC5PLiHEstimatorCallback>(), 5, threshold, prob)->run(
                points1, points2, models, _mask);
    else
        createLMeDSPointSetRegistrator(
                makePtr<PC5PLiHEstimatorCallback>(), 5, prob)->run(
                points1, points2, models, _mask);

    return models;

}

}
