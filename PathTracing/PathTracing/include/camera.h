#ifndef CAMERA_H_
#define CAMERA_H_

#include <Eigen/Dense>

namespace util {

typedef Eigen::Matrix<double, 3, 4> CameraMatrix;

struct Camera {
  Camera(const Eigen::Vector3d &eye, const Eigen::Vector2d &bottom,
         const Eigen::Vector2d &top, const double width, const double height);
  ~Camera() {}

  CameraMatrix camera_matrix_;
};

}  // namespace util

#endif  // CAMERA_H_