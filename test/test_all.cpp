#include <opencv2/opencv.hpp>
#include "relative_pose/relative_pose.hpp"

enum METHOD
{
    PC_5P = 0,
    PC_4PRA = 1,
    PC_4PST0 = 2
};

int main()
{
    double angle_bound = CV_PI / 18; 

    double nearest_dist = 10; 
    double baseline_dev = 0.05; 
    double baseline = -1; 
    double depth = 10; 

    double focal = 300; 
    double bound_2d = 175; 

    cv::RNG rng; 
    std::unordered_map<int, std::vector<cv::Mat>> t_angles;
    
    for (double stdd = 0.1; stdd <= 1.0; stdd += 0.1)
    {
        double sigma = stdd; 

        for (int i = 0; i < 1000; i++)
        {    
            cv::Mat rvec(3, 1, CV_64F), tvec(3, 1, CV_64F), cvec(3, 1, CV_64F); 
            rng.fill(rvec, cv::RNG::UNIFORM, -angle_bound, angle_bound); 
            rng.fill(cvec, cv::RNG::UNIFORM, -baseline_dev, baseline_dev); 
            
            cvec.at<double>(2) = baseline;  
            normalize(cvec, cvec); 

            cv::Mat rmat; 
            Rodrigues(rvec, rmat); 
            double angle = cv::norm(rvec);

            tvec = -rmat * cvec; 

            cv::Mat camera_matrix = (cv::Mat_<double>(3, 3) << focal, 0, 0, 0, focal, 0, 0, 0, 1); 
            
            cv::Mat object_points(5, 3, CV_64F); 
            rng.fill(object_points, cv::RNG::UNIFORM, -bound_2d, bound_2d); 
            object_points.col(2) = 1; 

            object_points = camera_matrix.inv() * object_points.t();
            object_points = object_points.t();

            cv::Mat ds_(5, 1, CV_64F), ds;
            rng.fill(ds_, cv::RNG::UNIFORM, nearest_dist, nearest_dist + depth); 
            cv::repeat(ds_, 1, 3, ds);
            cv::multiply(object_points, ds, object_points);

            cv::Mat image_points1, image_points2;
            cv::Mat zero3d = cv::Mat::zeros(1, 3, CV_64F);
            cv::projectPoints(object_points, zero3d, zero3d,
                    camera_matrix, cv::noArray(), image_points1);
            cv::projectPoints(object_points, rvec, tvec,
                    camera_matrix, cv::noArray(), image_points2);
            
            std::vector<cv::Point2f> vec = {{1, 2}, {2, 3}};

            cv::Mat u;
            cv::normalize(rvec, u);
            u *= std::sin(angle / 2);
            std::cerr << "u = " << u << std::endl;

            cv::Mat rvecs, tvecs, mask;
            estimateRelativePose_PC4PRA(angle, image_points1, image_points2, camera_matrix, 0, 0.99, 1,
                    rvecs, tvecs, mask);

//     std::cerr << "a4" << std::endl;
//             cv::Mat x1s = K * Xs.t(); 
//             cv::Mat x2s = rmat * Xs.t(); 
//             for (int j = 0; j < x2s.cols; j++) x2s.col(j) += tvec; 
//             x2s = K * x2s; 
// 
// 
//         
//             x1s.row(0) /= x1s.row(2); 
//             x1s.row(1) /= x1s.row(2); 
//             x1s.row(2) /= x1s.row(2); 
//         
//             x2s.row(0) /= x2s.row(2); 
//             x2s.row(1) /= x2s.row(2); 
//             x2s.row(2) /= x2s.row(2); 
//     std::cerr << "a5" << std::endl;
//     
//             x1s = x1s.t(); 
//             x2s = x2s.t(); 
//     
//             x1s = x1s.colRange(0, 2) * 1.0; 
//             x2s = x2s.colRange(0, 2) * 1.0; 
// 
//             x1s = x1s.reshape(2, 1);
//             x2s = x2s.reshape(2, 1);
//     std::cerr << "a6" << std::endl;
//     
//             cv::Mat noise1(x1s.size(), CV_64FC2), noise2(x2s.size(), CV_64FC2); 
//             rng.fill(noise1, cv::RNG::NORMAL, 0, sigma); 
//             cv::Mat x1s_noise = x1s + noise1; 
//             rng.fill(noise2, cv::RNG::NORMAL, 0, sigma); 
//             cv::Mat x2s_noise = x2s + noise2; 
// 
//     std::cerr << "a5" << std::endl;
//             tvec /= norm(tvec); 
//             std::cerr << "a1" << std::endl;
//             std::cerr << rvec << angle << std::endl;
//             std::cerr << tvec << std::endl;
// 
//             cv::Mat rvecs, tvecs, mask;
//             estimateRelativePose_PC4PRA(angle, x1s_noise, x2s_noise, K, 0, 0.99, 1,
//                     rvecs, tvecs, mask);

            exit(0);

        }

    }

}
