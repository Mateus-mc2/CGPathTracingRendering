#ifndef MATH_LIB_H_
#define MATH_LIB_H_

#include <cassert>
#include <cmath>

namespace math {

inline bool IsAlmostEqual(const double x, const double y, const double eps) {
  assert(eps < 1.0);
  return std::abs(x - y) <= std::abs(x)*eps;
}

}  // namespace math

#endif  // MATH_LIB_H_