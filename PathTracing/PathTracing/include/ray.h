#ifndef RAY_H_
#define RAY_H_

#include <Eigen\Dense>

namespace util {

struct Ray {
  Ray(const Eigen::Vector3d &origin, const Eigen::Vector3d &direction, int depth)
      : origin(origin),
        direction(direction),
        depth(depth) {}
  ~Ray() {}

  Eigen::Vector3d origin;
  Eigen::Vector3d direction;
  int depth;
};

}  // namespace util

#endif  // RAY_H_