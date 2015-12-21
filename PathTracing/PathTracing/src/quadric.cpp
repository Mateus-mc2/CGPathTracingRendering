#include "quadric.h"

using Eigen::Vector3d;
using Eigen::VectorXd;

namespace util {

Quadric::Quadric(const VectorXd &coefficients, const Material &material) : material_(material) {
  if (coefficients.size() != 10) {
    throw InvalidCoefficientsVectorException("Coefficients' vector doesn't have a valid size.");
  }
}

}  // namespace util