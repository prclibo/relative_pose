#include "mex.h"
#include "matrix.h"

LIBMMWMATRIX_PUBLISHED_API_EXTERN_C void *mxGetImagData(const mxArray *pa /* pointer to array */
                                                        );
#include <opencv2/matlab/bridge.hpp>
#include <opencv2/matlab/mxarray.hpp>

#include "relative_pose/relative_pose.hpp"
#include "../src/relative_pose/precomp.hpp"

#define _DO_NOT_EXPORT
#if defined(_DO_NOT_EXPORT)
#define DllExport  
#else
#define DllExport __declspec(dllexport)
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    using namespace cv;
    matlab::conditionalError(nrhs == 4, "Need 4 input arguments");

    matlab::MxArrayVector raw(prhs, prhs+nrhs);
    matlab::conditionalError(raw.at(2).size() == 1 && raw.at(3).size() == 1,
            "Input 2 or 3 is not scalar");

    matlab::ArgumentParser parser("what");
    parser.addVariant("estimateRelativePose_PC4PST0_NullE", 4, 0);

    matlab::MxArrayVector reordered = parser.parse(raw);
    bridge::BridgeVector inputs(reordered.begin(), reordered.end());

    Mat rays1 = inputs.at(0).toMat(),
        rays2 = inputs.at(1).toMat();
    double prob = inputs.at(2).toDouble(),
           thresh = inputs.at(3).toDouble();

    Mat E = estimateRelativePose_PC4PST0_NullE(rays1, rays2,
            RANSAC, prob, thresh, noArray());

    bridge::Bridge output;
    output = E;
    plhs[0] = output.toMxArray().releaseOwnership();
}

// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// {
//     using namespace cv;
//     assert(nrhs == 4);
//     Ptr<Mat> points1 = ocvMxArrayToImage_double(prhs[0], true);
//     Ptr<Mat> points2 = ocvMxArrayToImage_double(prhs[1], true);
//     double thresh = mxGetScalar(prhs[2]);
//     double prob = mxGetScalar(prhs[3]);
// 
//     Mat E = estimateRelativePose_PC4PST0_NullE(*points1, *points2,
//             RANSAC, prob, thresh, noArray());
// 
//     assert(nlhs == 1);
//     plhs[0] = ocvMxArrayFromMat_double(E);
// }
