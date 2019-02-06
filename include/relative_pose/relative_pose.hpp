#include <opencv2/opencv.hpp>

namespace cv
{


void estimateRelativePose_PC4PRA(double angle, InputArray _points1, InputArray _points2,
        InputArray _cameraMatrix, int method, double prob, double threshold,
        OutputArray _rvecs, OutputArray _tvecs, OutputArray _mask);
}
