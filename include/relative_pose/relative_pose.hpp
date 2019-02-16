#ifndef RELATIVE_POSE_HPP
#define RELATIVE_POSE_HPP
#include <opencv2/opencv.hpp>

namespace cv
{

static void processInputArray(InputArray _points1, InputArray _points2,
        InputArray _camera_matrix, double _thresh,
        Mat& points1, Mat& points2, Mat camera_matrix, double& thresh)
{
    _points1.getMat().convertTo(points1, CV_64F);
    _points2.getMat().convertTo(points2, CV_64F);
    _camera_matrix.getMat().convertTo(camera_matrix, CV_64F);

    int npoints = points1.checkVector(2);
    CV_Assert( npoints >= 0 && points2.checkVector(2) == npoints &&
                              points1.type() == points2.type());

    CV_Assert(camera_matrix.rows == 3 && camera_matrix.cols == 3 && camera_matrix.channels() == 1);

    if (points1.channels() > 1)
    {
        points1 = points1.reshape(1, npoints);
        points2 = points2.reshape(1, npoints);
    }

    double fx = camera_matrix.at<double>(0,0);
    double fy = camera_matrix.at<double>(1,1);
    double cx = camera_matrix.at<double>(0,2);
    double cy = camera_matrix.at<double>(1,2);

    points1.col(0) = (points1.col(0) - cx) / fx;
    points2.col(0) = (points2.col(0) - cx) / fx;
    points1.col(1) = (points1.col(1) - cy) / fy;
    points2.col(1) = (points2.col(1) - cy) / fy;

    std::cerr << "points1 = " << points1 << std::endl;
    std::cerr << "points2 = " << points2 << std::endl;

    // Reshape data to fit opencv ransac function
    points1 = points1.reshape(2, npoints);
    points2 = points2.reshape(2, npoints);

    thresh = _thresh / (fx + fy) * 2;

}

Mat estimateRelativePose_PC4PRA(double angle, InputArray _points1, InputArray _points2,
        InputArray _camera_matrix, int method, double prob, double threshold,
        OutputArray _mask);
void estimateRelativePose_PC4PRA(double angle, InputArray _points1, InputArray _points2,
        InputArray _camera_matrix, int method, double prob, double threshold,
        OutputArray _rvecs, OutputArray _tvecs, OutputArray _mask);

Mat estimateRelativePose_PC5P_LiH(InputArray _points1, InputArray _points2,
        InputArray _cameraMatrix, int method, double prob, double threshold,
        OutputArray _mask);
void estimateRelativePose_PC5P_LiH(InputArray _points1, InputArray _points2,
        InputArray _cameraMatrix, int method, double prob, double threshold,
        OutputArray _rvecs, OutputArray _tvecs, OutputArray _mask);
Mat estimateRelativePose_PC4PST0_NullE(
        InputArray _points1, InputArray _points2,
        InputArray _cameraMatrix, int method, double prob, double threshold,
        OutputArray _mask);
}
#endif
