#include <opencv2/opencv.hpp>
#include "gtest/gtest.h"
#include "../src/relative_pose/precomp.hpp"
#include "relative_pose/relative_pose.hpp"
#include "data_sampler.hpp"

double const E_ERR_THRESH = 1e-5;

TEST(PC_4PRA, Minimal)
{
    using namespace cv;

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
    EXPECT_LT(err, E_ERR_THRESH);
    EXPECT_EQ(countNonZero(mask), mask.total());
}

TEST(PC_4PRA, RANSAC)
{
    using namespace cv;

    DataSampler sampler;
    sampler.outlier_rate = 0.6;

    Vec3d rvec, tvec;
    std::vector<Point2d> image_points1, image_points2;
    sampler.sample(1000, rvec, tvec, image_points1, image_points2);
    Mat rmat;
    Rodrigues(rvec, rmat);
    Mat E0 = skew(tvec) * rmat;
    std::cerr << "E0 = " << E0 << std::endl;

    auto camera_matrix = sampler.cameraMatrix();
    Mat rvecs, tvecs, mask;
    Mat E = estimateRelativePose_PC4PRA(cv::norm(rvec),
            image_points1, image_points2,
            camera_matrix, cv::RANSAC, 0.99, 1e-2, mask);

    E0 /= E0.at<double>(2, 2);
    E /= E.at<double>(2, 2);
    double err = norm(E0, E) / norm(E0) / norm(E);

    EXPECT_LT(err, E_ERR_THRESH);
}

TEST(PC_5P_LiH, Minimal)
{
    using namespace cv;

    DataSampler sampler;
    Vec3d rvec, tvec;
    std::vector<Point2d> image_points1, image_points2;
    sampler.sample(5, rvec, tvec, image_points1, image_points2);

    Mat rmat;
    Rodrigues(rvec, rmat);
    Mat E0 = skew(tvec) * rmat;
    std::cerr << E0 / E0.at<double>(2, 2) << std::endl;

    auto camera_matrix = sampler.cameraMatrix();
    Mat rvecs, tvecs, mask;
    Mat E = estimateRelativePose_PC5P_LiH(image_points1, image_points2,
            camera_matrix, cv::RANSAC, 0.99, 1e-2, mask);

    E0 /= E0.at<double>(2, 2);
    double err = std::numeric_limits<double>::max();
    for (int r = 0; r < E.rows; r += 3)
    {
        Mat E1 = E.rowRange(r, r + 3);
        E1 /= E1.at<double>(2, 2);
        err = std::min(err, norm(E0, E1) / norm(E0) / norm(E1));
    }
    EXPECT_LT(err, E_ERR_THRESH);
    EXPECT_EQ(countNonZero(mask), mask.total());

}

TEST(PC_5P_LiH, RANSAC)
{
    using namespace cv;

    DataSampler sampler;
    sampler.outlier_rate = 0.6; 

    Vec3d rvec, tvec;
    std::vector<Point2d> image_points1, image_points2;
    sampler.sample(1000, rvec, tvec, image_points1, image_points2);
    Mat rmat;
    Rodrigues(rvec, rmat);
    Mat E0 = skew(tvec) * rmat;

    auto camera_matrix = sampler.cameraMatrix();
    Mat rvecs, tvecs, mask;
    Mat E = estimateRelativePose_PC5P_LiH(image_points1, image_points2,
            camera_matrix, cv::RANSAC, 0.99, 1e-2, mask);

    E0 /= E0.at<double>(2, 2);
    E /= E.at<double>(2, 2);
    double err = norm(E0, E) / norm(E0) / norm(E);

    EXPECT_LT(err, E_ERR_THRESH);
}


TEST(PC_4PST0_NullE, Minimal)
{
    using namespace cv;

    DataSampler sampler;
    Vec3d rvec, tvec;
    std::vector<Point2d> image_points1, image_points2;
    sampler.sample(4, rvec, tvec, image_points1, image_points2);
    Mat rmat;
    Rodrigues(rvec, rmat);
    Mat1d E0 = skew(tvec) * rmat;

    auto camera_matrix = sampler.cameraMatrix();
    Mat rvecs, tvecs, mask;
    Mat E = estimateRelativePose_PC4PST0_NullE(image_points1, image_points2,
            camera_matrix, cv::RANSAC, 0.99, 1e-2, mask);

    E0 /= E0(2, 2);
    double err = std::numeric_limits<double>::max();
    for (int r = 0; r < E.rows; r += 3)
    {
        Mat1d E1 = E.rowRange(r, r + 3);
        E1 /= E1(2, 2);
        err = std::min(err, norm(E0, E1) / norm(E0) / norm(E1));
    }
    EXPECT_LT(err, E_ERR_THRESH);
    EXPECT_EQ(countNonZero(mask), mask.total());
}

TEST(PC_4PST0_NullE, RANSAC)
{
    using namespace cv;

    DataSampler sampler;
    sampler.outlier_rate = 0.6;

    Vec3d rvec, tvec;
    std::vector<Point2d> image_points1, image_points2;
    sampler.sample(1000, rvec, tvec, image_points1, image_points2);
    Mat rmat;
    Rodrigues(rvec, rmat);
    Mat E0 = skew(tvec) * rmat;

    auto camera_matrix = sampler.cameraMatrix();
    Mat rvecs, tvecs, mask;
    Mat E = estimateRelativePose_PC5P_LiH(image_points1, image_points2,
            camera_matrix, cv::RANSAC, 0.99, 1e-2, mask);

    E0 /= E0.at<double>(2, 2);
    E /= E.at<double>(2, 2);
    double err = norm(E0, E) / norm(E0) / norm(E);

    EXPECT_LT(err, E_ERR_THRESH);
}
