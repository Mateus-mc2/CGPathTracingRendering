#ifndef LIGHT_H_
#define LIGHT_H_

#include <Eigen\Dense>

namespace util {

struct PointLight {
  PointLight(const Eigen::Vector3d &position, const double red, const double green,
             const double blue, const double intensity)
      : position(position),
        red(red),
        green(green),
        blue(blue),
        intensity(intensity) {}
  ~PointLight() {}
      
  Eigen::Vector3d position;

  double red;
  double green;
  double blue;

  double intensity;
};

}  // namespace util

#endif  // LIGHT_H_