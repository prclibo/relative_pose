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
    Ptr<Mat> points1 = ocvMxArrayToImage_uint8(prhs[0], true);
    std::cout << "points1 = " << points1->size() << std::endl;

    DataSampler sampler;
    Vec3d rvec, tvec;
    std::vector<Point2d> image_points1, image_points2;
    sampler.sample(4, rvec, tvec, image_points1, image_points2);

    Mat rmat;
    Rodrigues(rvec, rmat);
    Mat E0 = skew(tvec) * rmat;
    std::cerr << "E0 = " << E0 << std::endl;

    auto camera_matrix = sampler.cameraMatrix();
    std::cerr << "rvecs = " << rvec << std::endl;
    Mat rvecs, tvecs, mask;
    Mat E = estimateRelativePose_PC4PRA(cv::norm(rvec),
            image_points1, image_points2,
            camera_matrix, cv::RANSAC, 0.99, 1e-2, mask);

    std::cerr << "type = " << E0.type() << " " << E.type() << std::endl;
    std::cerr << "E0 = " << E0 << std::endl;
    E0 /= E0.at<double>(2, 2);
    double err = std::numeric_limits<double>::max();
    for (int r = 0; r < E.rows; r += 3)
    {
        Mat E1 = E.rowRange(r, r + 3);
        E1 /= E1.at<double>(2, 2);
        err = std::min(err, norm(E0, E1) / norm(E0) / norm(E1));
    }
}
