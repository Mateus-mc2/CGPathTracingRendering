#include "triangular_object.h"

using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix3d;

namespace util {

const double TriangularObject::kEps = 1.0e-03;

TriangularObject::TriangularObject(const Material &material,
                                   const std::vector<Eigen::Vector3d> &vertices,
                                   const std::vector<Eigen::Vector3i> &faces,
                                   const bool emissive_flag)
      :  kMaterial(material),
         kVertices(vertices),
         kFaces(faces),
         kEmissive(emissive_flag) {
  const int kNumFaces = faces.size();
  this->planes_coeffs_.resize(kNumFaces);
  this->linear_systems_.resize(kNumFaces);
    
  for (int i = 0; i < kNumFaces; ++i) {
    // Get plane equation (coefficients).
    const Vector3d kA = vertices[faces[i](0)];
    const Vector3d kB = vertices[faces[i](1)];
    const Vector3d kC = vertices[faces[i](2)];

    const Vector3d kAB = kB - kA;
    const Vector3d kAC = kC - kA;
    const Vector3d kNormal = kAB.cross(kAC);

    // The last coefficient is the additive inverse of the dot product of kA and kNormal.
    this->planes_coeffs_[i] << kNormal(0), kNormal(1), kNormal(2), -kNormal.dot(kA);

    // Get system's inverse matrix.
    Matrix3d linear_system;
    linear_system << kA(0), kB(0), kC(0),
                     kA(1), kB(1), kC(1), 
                     kA(2), kB(2), kC(2);
    this->linear_systems_[i] = linear_system.inverse();
  }
}

bool TriangularObject::IsInnerPoint(const Vector3d barycentric_coordinates) const {
  double alpha = barycentric_coordinates(0);
  double beta = barycentric_coordinates(1);
  double gamma = barycentric_coordinates(2);

  return (math::IsAlmostEqual(alpha + beta + gamma, 1.0, this->kEps) &&
          alpha >= 0 && beta >= 0 && gamma >= 0 &&
          alpha <= 1 && beta <= 1 && gamma <= 1);
}

double TriangularObject::GetIntersectionParameter(const Ray &ray) {
  Vector4d ray_origin(ray.origin(0), ray.origin(1), ray.origin(2), 1);

  // Parameters to return.
  double min_t = 0.0;
  Vector3d parameters;

  // Get nearest intersection point - need to check every single face of the object.
  for (int i = 0; i < this->planes_coeffs_.size(); ++i) {
    const double kNumerator = -(this->planes_coeffs_[i].dot(ray_origin));
    const double kDenominator = this->planes_coeffs_[i].dot(ray.direction);

    // Test if the ray and this plane are parallel (or if this plane contains the ray).
    // Returns a negative (dummy) parameter t if this happens.
    if (math::IsAlmostEqual(kDenominator, 0.0, this->kEps)) {
      return -1.0;
    }

    double curr_t = kNumerator / kDenominator;
    Vector3d intersection_point(ray.origin(0) + ray.direction(0)*curr_t,
                                ray.origin(1) + ray.direction(1)*curr_t,
                                ray.origin(2) + ray.direction(2)*curr_t);
    Vector3d barycentric_coords = this->linear_systems_[i]*intersection_point;  // x = A^(-1)*b.

    if (this->IsInnerPoint(barycentric_coords) && min_t > curr_t && curr_t > this->kEps) {
      min_t = curr_t;
      parameters = barycentric_coords;
    }
  }

  return min_t;
}

}  // namespace util