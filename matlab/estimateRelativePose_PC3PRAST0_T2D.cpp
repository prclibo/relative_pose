#include "mex.h"
#include "matrix.h"

LIBMMWMATRIX_PUBLISHED_API_EXTERN_C void *mxGetImagData(const mxArray *pa /* pointer to array */
                                                        );
#include <opencv2/matlab/bridge.hpp>
#include <opencv2/matlab/mxarray.hpp>

#include "relative_pose/relative_pose.hpp"
#include "../src/relative_pose/precomp.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    using namespace cv;
    matlab::conditionalError(nrhs == 5, "Need 5 input arguments");

    matlab::MxArrayVector raw(prhs, prhs+nrhs);
    matlab::conditionalError(raw.at(0).size() == 1 &&
            raw.at(3).size() == 1 && raw.at(4).size() == 1,
            "Input 0, 3 or 4 is not scalar");

    matlab::ArgumentParser parser("what");
    parser.addVariant("estimateRelativePose_PC3PRAST0_T2D", 5, 0);

    matlab::MxArrayVector reordered = parser.parse(raw);
    bridge::BridgeVector inputs(reordered.begin(), reordered.end());

    double angle = inputs.at(0).toDouble();
    Mat rays1 = inputs.at(1).toMat(),
        rays2 = inputs.at(2).toMat();
    double prob = inputs.at(3).toDouble(),
           thresh = inputs.at(4).toDouble();

    Mat E = estimateRelativePose_PC3PRAST0_T2D(angle, rays1, rays2,
            RANSAC, prob, thresh, noArray());

    bridge::Bridge output;
    output = E;
    plhs[0] = output.toMxArray().releaseOwnership();
}

