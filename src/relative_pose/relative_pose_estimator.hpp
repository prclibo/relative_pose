#ifndef RELATIVE_POSE_ESTIMATOR_HPP
#define RELATIVE_POSE_ESTIMATOR_HPP
#include <opencv2/opencv.hpp>
#include "precomp.hpp"

namespace cv
{

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

}
#endif
