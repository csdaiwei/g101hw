#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

float deg2rad(float angle){
    return angle / 180.0f * MY_PI;
}

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    Eigen::Matrix4f r;
    float rd = deg2rad(rotation_angle);
    r << std::cos(rd), -std::sin(rd), 0, 0, 
        std::sin(rd), std::cos(rd), 0, 0, 
        0, 0, 1, 0, 
        0, 0, 0, 1;

    model = r * model;

    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.
    float rd = deg2rad(eye_fov/2.0f);
    float t = zNear*std::tan(rd);
    float r = t * aspect_ratio;

    Eigen::Matrix4f m;
    // one way, works
    //m << zNear / r, 0, 0, 0, 0, zNear/t, 0, 0, 0, 0, -(zFar+zNear)/(zFar-zNear), -2*zFar*zNear/(zFar-zNear), 0, 0, -1, 0;
    
    //try another
    m << 1 / std::tan(rd)*aspect_ratio, 0, 0, 0, 0, 1/std::tan(rd), 0, 0, 0, 0, -(zFar+zNear)/(zFar-zNear), -2*zFar*zNear/(zFar-zNear), 0, 0, -1, 0;

    projection = m * projection;


    return projection;
}

Eigen::Matrix4f get_rotation(Vector3f axis, float angle)
{   
    float rd = deg2rad(angle);

    Eigen::Matrix3f M, N;
    Eigen::Vector3f n = axis.normalized();
    N << 0, -n.z(), n.y(), n.z(), 0, -n.x(), -n.y(), n.x(), 0;
    M = std::cos(rd) * Eigen::Matrix3f::Identity() + (1-std::cos(rd)) * (n * n.transpose()) + std::sin(rd) * N;

    Eigen::Matrix4f R = Eigen::Matrix4f::Identity();
    R.block(0, 0, 3, 3) = M;
    
    return R;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        //r.set_model(get_rotation({1.0f, 0.0f, 0.0f}, angle));
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, -0.1, -50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
