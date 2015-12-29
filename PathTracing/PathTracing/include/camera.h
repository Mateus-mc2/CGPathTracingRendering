#ifndef CAMERA_H_
#define CAMERA_H_

#include <Eigen/Dense>

namespace util {

typedef Eigen::Matrix<double, 3, 4> CameraMatrix;

struct Camera {
  Camera(const Eigen::Vector3d &eye,
         const Eigen::Vector2d &bottom,
         const Eigen::Vector2d &top,
         const double width,
         const double height)
      : eye_(eye),
        bottom_(bottom),
        top_(top),
        width_(width),
        height_(height) {}

  ~Camera() {}

  const Eigen::Vector3d eye_;
  const Eigen::Vector2d bottom_;
  const Eigen::Vector2d top_;
  const double width_;
  const double height_;
};

}  // namespace util

#endif  // CAMERA_H_