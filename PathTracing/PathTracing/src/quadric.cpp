#include "quadric.h"

using Eigen::Vector3d;
using Eigen::VectorXd;

namespace util {
  
const double Quadric::kEps = 1.0e-03;

Quadric::Quadric(const VectorXd &coefficients, const Material &material)
      : coefficients_(coefficients),
        material_(material) {
  if (coefficients.size() != 10) {
    throw InvalidCoefficientsVectorException("Coefficients' vector doesn't have a valid size.");
  }
}

double Quadric::GetIntesectionParameter(const Ray &ray) {
  // Coefficients.
  double a = this->coefficients_(this->kCoeffA);
  double b = this->coefficients_(this->kCoeffB);
  double c = this->coefficients_(this->kCoeffC);

  double f = this->coefficients_(this->kCoeffF);
  double g = this->coefficients_(this->kCoeffG);
  double h = this->coefficients_(this->kCoeffH);
  
  double p = this->coefficients_(this->kCoeffP);
  double q = this->coefficients_(this->kCoeffQ);
  double r = this->coefficients_(this->kCoeffR);
  double d = this->coefficients_(this->kCoeffD);

  // Ray parameters (origin and direction).
  double x_0 = ray.origin(0);
  double y_0 = ray.origin(1);
  double z_0 = ray.origin(2);

  double dx = ray.direction(0);
  double dy = ray.direction(1);
  double dz = ray.direction(2);
  
  double t;  // Parameter to return.

  // Equation coefficients (degree 2: kA*t^2 + kB*t + kC = 0).
  const double kA = a*dx*dx + b*dy*dy + c*dz*dz + 2*(f*dy*dz + g*dx*dz + h*dx*dy);
  const double kB = 2*(a*x_0*dx + b*y_0*dy + c*z_0*dz + f*(dz*y_0 + dy*z_0) + g*(dz*x_0 + dx*z_0)
                    + h*(dy*x_0 + dx*y_0) + p*dx + q*dy+ r*dz);
  const double kC = a* x_0*x_0 + b*y_0*y_0 + c*z_0*z_0 + d + 2*(f*y_0*z_0 + g*x_0*z_0 + h*x_0*y_0
                    + p*x_0 + q*y_0 + r*z_0);

  if (math::IsAlmostEqual(kA, 0.0, this->kEps)) {
    // The equation has degree 1.
    if (math::IsAlmostEqual(kB, 0.0, this->kEps)) {
      // The equation has degree 0, thus it's degenerate (it has infinite - or even zero - roots).
      return -1.0;
    }

    t = (-kC) / kB;
  } else {
    double discriminant = kB*kB - 4*kA*kC;

    if (discriminant < this->kEps) {
      // No real roots.
      return -1.0;
    }

    double sqrt_delta = std::sqrt(discriminant);
    // Gets the nearest point in front of the ray center.
    t = (-kB - sqrt_delta) / 2*kA;

    if (t < 0.0) {  // It is behind the ray center.
      t = (-kB + sqrt_delta) / 2*kA;
    }
  }

  if (t < this->kEps) {  // If it's negative or almost zero.
    return -1.0;
  }

  return t;
}

}  // namespace util