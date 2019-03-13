#ifndef RELATIVE_POSE_ESTIMATOR_HPP
#define RELATIVE_POSE_ESTIMATOR_HPP
#include <opencv2/opencv.hpp>
#include "precomp.hpp"

namespace cv
{
/*
class RelativePoseEstimatorCallback: public PointSetRegistrator::Callback
{
protected:
public:
    void computeError( InputArray _m1, InputArray _m2, InputArray _model,
            OutputArray _err ) const CV_OVERRIDE
    {
        Mat x1, x2;
        _m1.getMat().convertTo(x1, CV_32F);
        _m2.getMat().convertTo(x2, CV_32F);
        CV_Assert(x1.cols == 1 && x2.cols == 1);
        x1 = x1.reshape(1);
        x2 = x2.reshape(1);

        Mat model = _model.getMat(), E;
        model.convertTo(E, CV_32F);

        Mat1f x2tE = x2 * E,
              x1tE = x1 * E;
        Mat1f x1tE_sqr = x1tE.mul(x1tE),
              x2tE_sqr = x2tE.mul(x2tE);
        Mat1f x2tEx1 = x2tE.mul(x1);
        reduce(x2tEx1, x2tEx1, 1, REDUCE_SUM);
        Mat1f x2tEx1_sqr = x2tEx1.mul(x2tEx1);

        Mat1f errs;
        divide(x2tEx1_sqr, x1tE_sqr.col(0) + x1tE_sqr.col(1) +
                x2tE_sqr.col(0) + x2tE_sqr.col(1), errs);
        errs.convertTo(errs, CV_32F);
        errs = errs.t();
        _err.assign(errs);
    }
};
*/

static Mat1f simpleInverse(Mat1f const& m)
{
    // computes the inverse of a matrix m
    float det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
                 m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
                 m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    
    float invdet = 1.0 / det;
    
    Mat1f minv(3, 3); // inverse of matrix m
    minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
    minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
    minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
    minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
    minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
    minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
    minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
    minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
    minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
    return minv;
}
static void checkPositiveDepth( InputArray _m1, InputArray _m2, InputArray _model,
        Mat const& err_mask, OutputArray _mask)
{
    float dist_thresh_ = 50;
    bool zero_transl_ = true;

    Mat x1, x2;
    _m1.getMat().convertTo(x1, CV_32F);
    _m2.getMat().convertTo(x2, CV_32F);
    CV_Assert(x1.cols == 1 && x2.cols == 1);
    x1 = x1.reshape(1);
    x2 = x2.reshape(1);

    Mat model = _model.getMat(), E;
    model.convertTo(E, CV_32F);

    Mat R1, R2, t0;
    decomposeEssentialMat(E, R1, R2, t0);
    Mat best_mask;
    int best_count = 0;
    Mat best_X1s, best_X2s;
    for (int ri = 0; ri < 2; ++ri)
        for (int ti = 0; ti < 2; ++ti)
        {
            Mat1f R = (ri % 2 == 0 ? R1 : R2);
            Mat1f t = (ti % 2 == 0 ? t0 : -t0);
            Mat1f X1s(3, x1.rows), X2s(3, x1.rows);
            if (zero_transl_ && trace(R)[0] < -0.99) continue;
            for (int i = 0; i < x1.rows; ++i)
            {
                // std::cerr << "x1 = " << x1.row(i) << std::endl;
                // std::cerr << "x2 = " << x2.row(i) << std::endl;
                // std::cerr << "R = " << R << std::endl;
                // std::cerr << "t = " << t << std::endl;
                // std::cerr << x1.type() << " " << R.type() << std::endl;
                Mat1f A1 = skew(x1.row(i));
                // std::cerr << "A1 = " << A1 << std::endl;
                Mat1f A2 = skew(x2.row(i)) * R; // R can be integrated into x2.
                // std::cerr << "A2 = " << A2 << std::endl;
                Mat1f AtA = A1.t() * A1 + A2.t() * A2;
                // std::cerr << "AtA = " << AtA << std::endl;
                Mat1f AtAinv = simpleInverse(AtA);
                // std::cerr << "AtAinv = " << AtAinv << std::endl;
                Mat1f Atb = -A2.t() * skew(x2.row(i)) * t;
                // std::cerr << "Atb = " << Atb << std::endl;
                X1s.col(i) = AtAinv * Atb;                    
                // std::cerr << "X = " << X << std::endl;
                // exit(0);
                X2s.col(i) = R * X1s.col(i) + t;
            }
            Mat mask = ((X1s.row(2) > 0) & (X1s.row(2) < dist_thresh_)
                    & (X2s.row(2) > 0) & (X2s.row(2) < dist_thresh_));
            mask &= err_mask;
            int count = countNonZero(mask);
            if (count > best_count)
            {
                best_mask = mask;
                best_count = count;
                best_X1s = X1s;
                best_X2s = X2s;
            }
        }
    // std::cerr << "x1s" << x1.rowRange(0, 10) << std::endl;
    // std::cerr << "x2s" << x2.rowRange(0, 10) << std::endl;
    // std::cerr << "X1s = " << best_X1s.colRange(0, 10) << std::endl;
    // std::cerr << "X2s = " << best_X2s.colRange(0, 10) << std::endl;
    // std::cerr << "best_count = " << best_count << std::endl;
    _mask.assign(best_mask);
}

class RelativePoseEstimatorCallback: public PointSetRegistrator::Callback
{
protected:
    
public:
    RelativePoseEstimatorCallback() {}
    void computeError( InputArray _m1, InputArray _m2, InputArray _model,
            OutputArray _err ) const CV_OVERRIDE
    {
        Mat x1, x2;
        _m1.getMat().convertTo(x1, CV_32F);
        _m2.getMat().convertTo(x2, CV_32F);
        CV_Assert(x1.cols == 1 && x2.cols == 1);
        x1 = x1.reshape(1);
        x2 = x2.reshape(1);

        Mat model = _model.getMat(), E;
        model.convertTo(E, CV_32F);

        Mat1f x2tE = x2 * E,
              x1tEt = x1 * E.t();
        Mat1f x1tEt_sqr = x1tEt.mul(x1tEt),
              x2tE_sqr = x2tE.mul(x2tE);
        Mat1f x2tEx1 = x2tE.mul(x1);
        reduce(x2tEx1, x2tEx1, 1, REDUCE_SUM);
        Mat1f x2tEx1_sqr = x2tEx1.mul(x2tEx1);

        Mat1f errs;
        divide(x2tEx1_sqr, x1tEt_sqr.col(0) + x1tEt_sqr.col(1) +
                x2tE_sqr.col(0) + x2tE_sqr.col(1), errs);
        divide(x2tEx1_sqr,
                x1tEt_sqr.col(0) + x1tEt_sqr.col(1) + x1tEt_sqr.col(2) +
                x2tE_sqr.col(0) + x2tE_sqr.col(1) + x2tE_sqr.col(2) -
                2.0 * x2tEx1_sqr, errs);
        errs.convertTo(errs, CV_32F);
        errs = errs.t();
        _err.assign(errs);

        // std::cerr << "E = " << E << std::endl;
        // std::cerr << "x1 = " << x1.row(0) << std::endl;
        // std::cerr << "x2 = " << x2.row(0) << std::endl;
        // std::cerr << "x1tEt = " << x1tEt(0) << std::endl;
        // std::cerr << "x2tE = " << x2tE(0) << std::endl;
        // std::cerr << "x2tEx1 = " << x2tEx1(0) << std::endl;
        // std::cerr << "x1tEt_sqr = " << x1tEt_sqr(0) << std::endl;
        // std::cerr << "x2tE_sqr = " << x2tE_sqr(0) << std::endl;
        // std::cerr << "x2tEx1_sqr = " << x2tEx1_sqr(0) << std::endl;
        // std::cerr << "errs = " << errs(0) << std::endl;

    }
};



}
#endif
