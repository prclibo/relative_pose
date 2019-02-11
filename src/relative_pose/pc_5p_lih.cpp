#include <opencv2/opencv.hpp>
#include "relative_pose/relative_pose.hpp"
#include "precomp.hpp"

void compute_E_matrices (double q[][3], double qp[][3],
        double Ematrices[10][3][3], int &nroots, bool optimized);

namespace cv
{

class PC5PLiHEstimatorCallback CV_FINAL : public PointSetRegistrator::Callback
{
protected:
    float dist_thresh_;

public:
    PC5PLiHEstimatorCallback(float dist_thresh = 100)
        : dist_thresh_(dist_thresh) {}

    int runKernel( InputArray _m1, InputArray _m2, OutputArray _model ) const CV_OVERRIDE
    {
        Mat2d q1 = _m1.getMat(), q2 = _m2.getMat();
        std::cerr << "q1 = " << q1 << std::endl;
        std::cerr << "q2 = " << q2 << std::endl;
        CV_Assert(q1.type() == CV_64FC2 && q2.type() == CV_64FC2);
        CV_Assert(q1.cols == 1 && q2.cols == 1);

        Mat _q1 = q1.reshape(1), _q2 = q2.reshape(1);
        std::cerr << "q1.size() = " << _q1.size() << std::endl;
        Mat1d mag1, mag2;
        magnitude(_q1.col(0), _q1.col(1), mag1);
        magnitude(_q2.col(0), _q2.col(1), mag2);
        // _q1 /= repeat(mag1, 1, 2);
        // _q2 /= repeat(mag2, 1, 2);
        // q1 = _q1.reshape(2); q2 = _q2.reshape(2);


        double q[5][3], qp[5][3];
        for (int si = 0; si < 5; ++si)
        {
            q[si][0] = q1(si, 0)[0] / mag1(si, 0);
            q[si][1] = q1(si, 0)[1] / mag1(si, 0);
            q[si][2] = 1 / mag1(si, 0);
            qp[si][0] = q2(si, 0)[0] / mag2(si, 0);
            qp[si][1] = q2(si, 0)[1] / mag2(si, 0);
            qp[si][2] = 1 / mag2(si, 0);
        }

        double ematrices[10][3][3];
        int num_roots;
        compute_E_matrices(q, qp, ematrices, num_roots, false);

        std::cerr << "num_roots = " << num_roots << std::endl;
        for (int i = 0; i < num_roots; ++i)
        {
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    std::cerr << ematrices[i][r][c] / ematrices[i][2][2] << " ";
            std::cerr << std::endl;
        }

        return 0;

        // Mat1d model;
        // for (int i = 0; i < ns; ++i)
        // {
        //     Vec3d rvec(rbuf + i * 3), tvec(tbuf + i * 3);
        //     // FIXME(li): Evgeniy's 4p2v code seems using transposed
        //     // quaternion<->rotation matrix representation instead of the
        //     // wikipedia convention.
        //     // https://math.stackexchange.com/questions/383754/are-there-different-conventions-for-representing-rotations-as-quaternions.
        //     rvec *= -1;

        //     model.push_back(rvec);
        //     model.push_back(tvec);
        //     model.push_back(rvec);
        //     model.push_back(-tvec);
        // }
        // model = model.reshape(1);
        // _model.assign(model);

        // return model.rows / 2;
    }

    void computeError( InputArray _m1, InputArray _m2, InputArray _model,
            OutputArray _err ) const CV_OVERRIDE
    {
        Mat x1, x2;
        _m1.getMat().convertTo(x1, CV_32F);
        _m2.getMat().convertTo(x2, CV_32F);
        Mat x1h, x2h;
        cv::convertPointsToHomogeneous(x1, x1h);
        cv::convertPointsToHomogeneous(x2, x2h);
        x1h = x1h.reshape(1);
        x2h = x2h.reshape(1);

        Mat model = _model.getMat(), rmat, rvec, tvec;
        model.row(0).convertTo(rvec, CV_32F);
        model.row(1).convertTo(tvec, CV_32F);
        Rodrigues(rvec, rmat);
        tvec = tvec.t();

        Mat1f E = skew(tvec) * rmat;
        Mat1f x2tE = x2h * E,
              x1tE = x1h * E;
        Mat1f x1tE_sqr = x1tE.mul(x1tE),
              x2tE_sqr = x2tE.mul(x2tE);
        Mat1f x2tEx1 = x2tE.mul(x1h);
        reduce(x2tEx1, x2tEx1, 1, REDUCE_SUM);
        Mat1f x2tEx1_sqr = x2tEx1.mul(x2tEx1);

        Mat1f errs;
        divide(x2tEx1_sqr, x1tE_sqr.col(0) + x1tE_sqr.col(1) +
                x2tE_sqr.col(0) + x2tE_sqr.col(1), errs);
        errs.convertTo(errs, CV_32F);
        errs = errs.t();

        // Cheirality check. c.f. cv::recovePose().
        Mat1f P1 = Mat1d::eye(3, 4), P2(3, 4);
        rmat.copyTo(P2.colRange(0, 3));
        tvec.copyTo(P2.col(3));
        Mat1f tlated(0, 0);
        triangulatePoints(P1, P2, x1, x2, tlated);

        tlated = tlated / repeat(tlated.row(3), 4, 1);
        auto mask1 = (tlated.row(2) > 0) & (tlated.row(2) < dist_thresh_);
        tlated = P2 * tlated;
        auto mask2 = (tlated.row(2) > 0) & (tlated.row(2) < dist_thresh_);

        errs.setTo(std::numeric_limits<float>::quiet_NaN(), mask1 | mask2);

        _err.assign(errs);
    }
};

void estimateRelativePose_PC5P_LiH(InputArray _points1, InputArray _points2,
        InputArray _cameraMatrix, int method, double prob, double threshold,
        OutputArray _rvecs, OutputArray _tvecs, OutputArray _mask)
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

    // Mat1d rvecs, tvecs;
    // for (int i = 0; i < models.rows; i += 2)
    // {
    //     rvecs.push_back(models.row(i) * 1.0);
    //     tvecs.push_back(models.row(i + 1) * 1.0);
    // }
    // _rvecs.assign(rvecs);
    // _tvecs.assign(tvecs);
}

}
