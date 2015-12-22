#include "material.h"

namespace util {

Material::Material(const double &r, const double &g, const double &b, const double &k_a,
                   const double &k_d, const double &k_s, const int &n)
                   :  red(r), green(g), blue(b), k_a(k_a), k_d(k_d), k_s(k_s), n(n) {
  if (r < 0 || g < 0 || b < 0 || k_a < 0 || k_d < 0 || k_s < 0 || n <= 0 ||
      r > 1 || g > 1 || b > 1 || k_a > 1 || k_d > 1 || k_s > 1) {
    throw InvalidMaterialCoefficientsException("Material coefficients out of range.");
  }
}

}  // namespace util