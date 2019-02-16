#include "mex.h"
#include "opencvmex.hpp"

#include "relative_pose/relative_pose.hpp"
#include "../src/relative_pose/precomp.hpp"
#include "../test/data_sampler.hpp"

#define _DO_NOT_EXPORT
#if defined(_DO_NOT_EXPORT)
#define DllExport  
#else
#define DllExport __declspec(dllexport)
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    using namespace cv;
    assert(nrhs == 4);
    Ptr<Mat> points1 = ocvMxArrayToImage_double(prhs[0], true);
    std::cout << "points1 = " << *points1 << std::endl;
    Ptr<Mat> points2 = ocvMxArrayToImage_double(prhs[1], true);
    std::cout << "points2 = " << *points2 << std::endl;
    double thresh = mxGetScalar(prhs[2]);
    double prob = mxGetScalar(prhs[3]);
    std::cout << "thresh = " << thresh << std::endl;

    estimateRelativePose_PC4PST0_NullE(*points1, *points2, Mat1d::eye(3, 3),
            RANSAC, prob, thresh, noArray());

}
