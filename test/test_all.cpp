#include <chrono>
#include <opencv2/opencv.hpp>
#include "relative_pose/relative_pose.hpp"
#include "data_sampler.hpp"

enum METHOD
{
    PC_5P = 0,
    PC_4PRA = 1,
    PC_4PST0 = 2
};

int const LOOPS = 1000;

int main()
{
    using namespace cv;

    std::chrono::duration<double> ttime_pc5p_lih(0),
                                  ttime_pc4pra(0),
                                  ttime_pc4pst0(0),
                                  ttime_pc3prast0(0);
    for (double stdd = 0.1; stdd <= 1.0; stdd += 0.1)
    {
        double sigma = stdd; 
        DataSampler sampler;

        for (int i = 0; i < LOOPS; i++)
        {    
            Vec3d rvec, tvec;
            std::vector<Point3d> image_rays1, image_rays2;
            sampler.sampleRays(5, rvec, tvec, image_rays1, image_rays2);

            {
                cv::Mat rvecs, tvecs, mask;
                auto start = std::chrono::system_clock::now();
                Mat E = estimateRelativePose_PC5P_LiH(image_rays1, image_rays2,
                        cv::RANSAC, 0.99, 1e-2, mask);
                auto end = std::chrono::system_clock::now();
                ttime_pc5p_lih += (end - start);
            }

            image_rays1.resize(4); image_rays2.resize(4);
            {
                cv::Mat rvecs, tvecs, mask;
                auto start = std::chrono::system_clock::now();
                Mat E = estimateRelativePose_PC4PRA(cv::norm(rvec),
                        image_rays1, image_rays2,
                        cv::RANSAC, 0.99, 1e-2, mask);
                auto end = std::chrono::system_clock::now();
                ttime_pc4pra += (end - start);
            }
            {
                cv::Mat rvecs, tvecs, mask;
                auto start = std::chrono::system_clock::now();
                Mat E = estimateRelativePose_PC4PST0_NullE_Poly(
                        image_rays1, image_rays2,
                        cv::RANSAC, 0.99, 1e-2, mask);
                auto end = std::chrono::system_clock::now();
                ttime_pc4pst0 += (end - start);
            }
            image_rays1.resize(3); image_rays2.resize(3);
            {
                cv::Mat rvecs, tvecs, mask;
                auto start = std::chrono::system_clock::now();
                Mat E = estimateRelativePose_PC3PRAST0_T2D(cv::norm(rvec),
                        image_rays1, image_rays2,
                        cv::RANSAC, 0.99, 1e-2, mask);
                auto end = std::chrono::system_clock::now();
                ttime_pc3prast0 += (end - start);
            }


            // estimateRelativePose_PC4PRA(angle,
            //         image_rays1.rowRange(0, 4), image_rays2.rowRange(0, 4),
            //         camera_matrix, cv::RANSAC, 0.99, 1, rvecs, tvecs, mask);
            // std::cerr << rvecs << std::endl;
            // std::cerr << rvec << std::endl;
            // auto end = std::chrono::system_clock::now();
            // total_time += (end - start);

            // std::cerr << rvecs << tvecs << std::endl;


        }
        break;

    }

    std::cout << ttime_pc5p_lih.count() / LOOPS << std::endl;
    std::cout << ttime_pc4pra.count() / LOOPS << std::endl;
    std::cout << ttime_pc4pst0.count() / LOOPS << std::endl;
    std::cout << ttime_pc3prast0.count() / LOOPS << std::endl;

}
