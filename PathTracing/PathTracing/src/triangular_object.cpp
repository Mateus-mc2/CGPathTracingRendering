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

// TODO(Mateus): avaliar se é melhor calcular a interseção de um raio e uma dada face,
// em vez de calcular a interseção de um raio com todas as faces do objeto.
Vector3d TriangularObject::GetIntersectionParameters(const Ray &ray) {
  // Ray parameters (origin and direction).
  double x_0 = ray.origin(0);
  double y_0 = ray.origin(1);
  double z_0 = ray.origin(2);

  double dx = ray.direction(0);
  double dy = ray.direction(1);
  double dz = ray.direction(2);

  // Parameters to return.
  double min_t = 0.0;
  Vector3d parameters;

  // Get nearest intersection point - need to check every single face of the object.
  for (int i = 0; i < this->planes_coeffs_.size(); ++i) {
    // Plane coefficients.
    const double kA = this->planes_coeffs_[i](0);
    const double kB = this->planes_coeffs_[i](1);
    const double kC = this->planes_coeffs_[i](2);
    const double kD = this->planes_coeffs_[i](3);

    double curr_t = - (kA*x_0 + kB*y_0 + kC*z_0 + kD) / (kA*dx + kB*dy + kC*dz);
    Vector3d intersection_point(x_0 + dx*curr_t, y_0 + dy*curr_t, z_0 + dz*curr_t);
    Vector3d barycentric_coords = this->linear_systems_[i]*intersection_point;

    if (this->IsInnerPoint(barycentric_coords) && min_t > curr_t && curr_t > this->kEps) {
      min_t = curr_t;
      parameters = barycentric_coords;
    }
  }

  return parameters;
}

}  // namespace util