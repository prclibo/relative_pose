#include <opencv2/opencv.hpp>
#include "gtest/gtest.h"
#include "../src/relative_pose/precomp.hpp"
#include "relative_pose/relative_pose.hpp"
#include "data_sampler.hpp"

double const E_ERR_THRESH = 1e-5;
double const RAY_ERR_THRESH = 1e-5;

using namespace cv;

class ScrewPlanarMotionTest: public ::testing::Test
{
protected:
    void setup(int num_rays, float outlier_rate)
    {
        DataSampler sampler;
        sampler.outlier_rate = outlier_rate;
        sampler.zero_screw_transl = true;
        sampler.sampleRays(num_rays, rvec_, tvec_, image_rays1_, image_rays2_);

        Rodrigues(rvec_, rmat_);
        E0_ = skew(tvec_) * rmat_;

        camera_matrix = sampler.cameraMatrix();
    }
    std::vector<Point3d> image_rays1_, image_rays2_;
    Vec3d rvec_, tvec_;
    Mat rmat_, E0_, E_, camera_matrix, mask_;
};

void expectEqualE(Mat E, Mat E0)
{
    assert(E.type() == CV_64F && E0.type() == CV_64F);
    E0 /= E0.at<double>(2, 2);
    double err = std::numeric_limits<double>::max();
    for (int r = 0; r < E.rows; r += 3)
    {
        Mat E1 = E.rowRange(r, r + 3);
        E1 /= E1.at<double>(2, 2);
        err = std::min(err, norm(E0, E1) / norm(E0) / norm(E1));
    }
    EXPECT_LT(err, E_ERR_THRESH);
}

// ---------------------------------------------------------------------------

TEST_F(ScrewPlanarMotionTest, PC_4PRA_Minimal_StdVectorPoint3d)
{
    setup(4, 0);
    E_ = estimateRelativePose_PC4PRA(norm(rvec_),
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);

    expectEqualE(E_, E0_);
    EXPECT_EQ(countNonZero(mask_), mask_.total());
}

TEST_F(ScrewPlanarMotionTest, PC_4PRA_Minimal_Mat3d)
{
    setup(4, 0);
    Mat3d image_rays1 = Mat(image_rays1_),
          image_rays2 = Mat(image_rays2_);

    E_ = estimateRelativePose_PC4PRA(norm(rvec_),
            image_rays1, image_rays2, RANSAC, 0.99, RAY_ERR_THRESH, mask_);

    expectEqualE(E_, E0_);
    EXPECT_EQ(countNonZero(mask_), mask_.total());
}

TEST_F(ScrewPlanarMotionTest, PC_4PRA_RANSAC_StdVectorPoint3d)
{
    setup(1000, 0.6);
    E_ = estimateRelativePose_PC4PRA(norm(rvec_),
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);
    expectEqualE(E_, E0_);
}

TEST_F(ScrewPlanarMotionTest, PC_4PRA_RANSAC_Mat3d)
{
    setup(1000, 0.6);
    Mat3d image_rays1 = Mat(image_rays1_),
          image_rays2 = Mat(image_rays2_);

    E_ = estimateRelativePose_PC4PRA(norm(rvec_),
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);
    expectEqualE(E_, E0_);
}

// ---------------------------------------------------------------------------

TEST_F(ScrewPlanarMotionTest, PC_5P_LiH_Minimal_StdVectorPoint3d)
{
    setup(5, 0);
    E_ = estimateRelativePose_PC5P_LiH(
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);

    expectEqualE(E_, E0_);
    EXPECT_EQ(countNonZero(mask_), mask_.total());
}

TEST_F(ScrewPlanarMotionTest, PC_5P_LiH_Minimal_Mat3d)
{
    setup(5, 0);
    Mat3d image_rays1 = Mat(image_rays1_),
          image_rays2 = Mat(image_rays2_);

    E_ = estimateRelativePose_PC5P_LiH(
            image_rays1, image_rays2, RANSAC, 0.99, RAY_ERR_THRESH, mask_);

    expectEqualE(E_, E0_);
    EXPECT_EQ(countNonZero(mask_), mask_.total());
}

TEST_F(ScrewPlanarMotionTest, PC_5P_LiH_RANSAC_StdVectorPoint3d)
{
    setup(1000, 0.6);
    E_ = estimateRelativePose_PC5P_LiH(
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);
    expectEqualE(E_, E0_);
}

// ---------------------------------------------------------------------------

TEST_F(ScrewPlanarMotionTest, PC_4PST0_NullE_Minimal_StdVectorPoint3d)
{
    setup(4, 0);
    E_ = estimateRelativePose_PC4PST0_NullE(
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);

    expectEqualE(E_, E0_);
    EXPECT_EQ(countNonZero(mask_), mask_.total());
}

TEST_F(ScrewPlanarMotionTest, PC_4PST0_NullE_LiH_RANSAC_StdVectorPoint3d)
{
    setup(1000, 0.6);
    E_ = estimateRelativePose_PC4PST0_NullE(
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);
    expectEqualE(E_, E0_);
}

// ---------------------------------------------------------------------------

TEST_F(ScrewPlanarMotionTest, PC_3PRAST0_T2D_Minimal_StdVectorPoint3d)
{
    setup(3, 0);
    E_ = estimateRelativePose_PC3PRAST0_T2D(norm(rvec_),
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);

    expectEqualE(E_, E0_);
    EXPECT_EQ(countNonZero(mask_), mask_.total());
}

TEST_F(ScrewPlanarMotionTest, PC_3PRAST0_T2D_RANSAC_StdVectorPoint3d)
{
    setup(1000, 0.6);
    E_ = estimateRelativePose_PC3PRAST0_T2D(norm(rvec_),
            image_rays1_, image_rays2_, RANSAC, 0.99, RAY_ERR_THRESH, mask_);
    expectEqualE(E_, E0_);
}


// TEST(PC_4PST0_NullE, Minimal)
// {
//     using namespace cv;
// 
//     DataSampler sampler;
//     Vec3d rvec, tvec;
//     std::vector<Point2d> image_points1, image_points2;
//     sampler.sample(4, rvec, tvec, image_points1, image_points2);
//     Mat rmat;
//     Rodrigues(rvec, rmat);
//     Mat1d E0 = skew(tvec) * rmat;
// 
//     auto camera_matrix = sampler.cameraMatrix();
//     Mat rvecs, tvecs, mask;
//     Mat E = estimateRelativePose_PC4PST0_NullE(image_points1, image_points2,
//             camera_matrix, RANSAC, 0.99, RAY_ERR_THRESH, mask);
// 
//     E0 /= E0(2, 2);
//     double err = std::numeric_limits<double>::max();
//     for (int r = 0; r < E.rows; r += 3)
//     {
//         Mat1d E1 = E.rowRange(r, r + 3);
//         E1 /= E1(2, 2);
//         err = std::min(err, norm(E0, E1) / norm(E0) / norm(E1));
//     }
//     EXPECT_LT(err, E_ERR_THRESH);
//     EXPECT_EQ(countNonZero(mask), mask.total());
// }
// 
// TEST(PC_4PST0_NullE, RANSAC)
// {
//     using namespace cv;
// 
//     DataSampler sampler;
//     sampler.outlier_rate = 0.6;
// 
//     Vec3d rvec, tvec;
//     std::vector<Point2d> image_points1, image_points2;
//     sampler.sample(1000, rvec, tvec, image_points1, image_points2);
//     Mat rmat;
//     Rodrigues(rvec, rmat);
//     Mat E0 = skew(tvec) * rmat;
// 
//     auto camera_matrix = sampler.cameraMatrix();
//     Mat rvecs, tvecs, mask;
//     Mat E = estimateRelativePose_PC4PST0_NullE(image_points1, image_points2,
//             camera_matrix, RANSAC, 0.99, RAY_ERR_THRESH, mask);
// 
//     E0 /= E0.at<double>(2, 2);
//     E /= E.at<double>(2, 2);
//     double err = norm(E0, E) / norm(E0) / norm(E);
// 
//     EXPECT_LT(err, E_ERR_THRESH);
// }
