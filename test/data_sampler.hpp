#include <opencv2/opencv.hpp>

double ANGLE_EPSILON = 1e-5;

struct DataSampler
{
    bool zero_screw_transl;
    double angle_bound, nearest_dist, baseline_turb, depth,
           focal, half_size, noise, outlier_rate;
    cv::Vec3d baseline;
    cv::RNG rng;
    DataSampler();
    void sampleRays(int num_points, cv::OutputArray _rvec, cv::OutputArray _tvec,
            cv::OutputArray _image_rays1, cv::OutputArray _image_rays2);

    cv::Mat1d cameraMatrix() const;
    void samplePoints(int num_points, cv::OutputArray _rvec, cv::OutputArray _tvec,
            cv::OutputArray _image_points1, cv::OutputArray _image_points2);
};

DataSampler::DataSampler()
{
    rng = cv::RNG(211);
    angle_bound = CV_PI / 18; 
    nearest_dist = 10; 
    baseline_turb = 0.05; 
    baseline = {0, 0, -1};
    depth = 10; 
    focal = 300; 
    half_size = 175; 
    noise = 0;
    outlier_rate = 0;
    zero_screw_transl = true;
}

cv::Mat1d DataSampler::cameraMatrix() const
{
    cv::Mat1d camera_matrix(3, 3);
    camera_matrix << focal, 0, 0, 0, focal, 0, 0, 0, 1;
    return camera_matrix;
}

void assignOutputArray(cv::Mat mat, cv::OutputArray output)
{
    if (output.fixedType())
    {
        mat = mat.reshape(output.channels());
        mat.convertTo(mat, output.type());
    }

    output.create(mat.size(), mat.type());

    // FIXME(libo): Work-around for https://github.com/opencv/opencv/issues/5350
    if (output.kind() == cv::_InputArray::STD_VECTOR)
        mat = mat.t();

    mat.copyTo(output.getMat());

}

void DataSampler::sampleRays(int num_rays,
        cv::OutputArray _rvec, cv::OutputArray _tvec,
        cv::OutputArray _image_rays1, cv::OutputArray _image_rays2)
{
    using namespace cv;

    Mat2d image_points1, image_points2;
    samplePoints(num_rays, _rvec, _tvec, image_points1, image_points2);

    Mat3d image_rays1, image_rays2;
    convertPointsToHomogeneous(image_points1.reshape(1) / focal, image_rays1);
    convertPointsToHomogeneous(image_points2.reshape(1) / focal, image_rays2);

    assignOutputArray(image_rays1, _image_rays1);
    assignOutputArray(image_rays2, _image_rays2);
}


void DataSampler::samplePoints(int num_points,
        cv::OutputArray _rvec, cv::OutputArray _tvec,
        cv::OutputArray _image_points1, cv::OutputArray _image_points2)
{
    using namespace cv;
    Vec3d rvec, tvec, cvec;
    rng.fill(rvec, RNG::UNIFORM, -angle_bound, angle_bound); 
    rng.fill(cvec, RNG::UNIFORM, -baseline_turb, baseline_turb); 
    cvec += baseline;
    double angle = norm(rvec);

    Matx33d rmat; 
    Rodrigues(rvec, rmat); 
    tvec = -rmat * cvec;
    if (zero_screw_transl && angle > ANGLE_EPSILON)
        tvec -= (rvec / angle) * (tvec.dot(rvec) / angle);

    auto camera_matrix = cameraMatrix();

    Mat1d object_points(num_points, 3);
    rng.fill(object_points, RNG::UNIFORM, -half_size, half_size); 
    object_points.col(2) = 1; 
    object_points = camera_matrix.inv() * object_points.t();
    object_points = object_points.t();

    Mat1d ds(num_points, 1);
    rng.fill(ds, RNG::UNIFORM, nearest_dist, nearest_dist + depth); 
    multiply(object_points, repeat(ds, 1, 3), object_points);

    Mat image_points1, image_points2;
    projectPoints(object_points, Matx13d::zeros(), Matx13d::zeros(),
            camera_matrix, noArray(), image_points1);
    projectPoints(object_points, rvec, tvec,
            camera_matrix, noArray(), image_points2);

    int num_outliers = std::floor(num_points * outlier_rate);
    if (num_outliers > 0)
    {
    rng.fill(image_points1.rowRange(0, num_outliers), RNG::UNIFORM,
            -half_size, half_size);
    rng.fill(image_points2.rowRange(0, num_outliers), RNG::UNIFORM,
            -half_size, half_size);
    }
    if (_image_points1.fixedType())
    {
        assert(_image_points2.fixedType() &&
                image_points1.type() == image_points2.type());
        image_points1 = image_points1.reshape(_image_points1.channels());
        image_points1.convertTo(image_points1, _image_points1.type());
        image_points2 = image_points2.reshape(_image_points2.channels());
        image_points2.convertTo(image_points2, _image_points2.type());
    }

    assignOutputArray(image_points1, _image_points1);
    assignOutputArray(image_points2, _image_points2);

    _rvec.assign(Mat(rvec));
    _tvec.assign(Mat(tvec));
}

