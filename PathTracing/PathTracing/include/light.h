#ifndef LIGHT_H_
#define LIGHT_H_

#include <Eigen\Dense>

namespace util {

struct Light {
  Light(const Eigen::Vector3d &direction, double intensity)
      : direction(direction),
        intensity(intensity) {}
  ~Light() {}
      
  Eigen::Vector3d direction;
  double intensity;
};

}  // namespace util

#endif  // LIGHT_H_