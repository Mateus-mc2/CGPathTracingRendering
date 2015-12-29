#ifndef TRIANGULAR_OBJECT_H_
#define TRIANGULAR_OBJECT_H_

#include <Eigen/Dense>

#include <vector>

#include "math_lib.h"
#include "material.h"
#include "ray.h"

namespace util {

class TriangularObject {
  public:
    TriangularObject(const Material &material, const std::vector<Eigen::Vector3d> &vertices,
                     const std::vector<Eigen::Vector3i> &faces, const bool emissive_flag);
    ~TriangularObject() {}

    Eigen::Vector3d GetIntersectionParameters(const Ray &ray);
  private:
    static const double kEps;

    const Material kMaterial;
    const std::vector<Eigen::Vector3d> kVertices;
    const std::vector<Eigen::Vector3i> kFaces;
    const bool kEmissive;

    bool IsInnerPoint(const Eigen::Vector3d barycentric_coordinates) const;

    std::vector<Eigen::Vector4d> planes_coeffs_;
    std::vector<Eigen::Matrix3d> linear_systems_;  // Inverse system matrices.
};

}  // namespace util

#endif  // TRIANGULAR_OBJECT_H_