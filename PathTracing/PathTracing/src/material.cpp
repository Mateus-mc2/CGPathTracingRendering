#include "material.h"

namespace util {

Material::Material(const double &r, const double &g, const double &b, const double &refraction_coeff,
                   const double &k_a, const double &k_d, const double &k_s, const double &k_t,
                   const int &n, const double &lp)
                   :  red(r),
                      green(g),
                      blue(b),
                      refraction_coeff(refraction_coeff),
                      k_a(k_a),
                      k_d(k_d),
                      k_s(k_s),
                      k_t(k_t),
                      n(n),
                      lp(lp) {
  if (r < 0 || g < 0 || b < 0 || k_a < 0 || k_d < 0 || k_s < 0 || k_t < 0 || n <= 0 ||
      r > 1 || g > 1 || b > 1 || k_a > 1 || k_d > 1 || k_s > 1 || k_t > 1 || lp < 0 ||
      lp > 1 || refraction_coeff < 1) {  // TODO: colocar teste para lp
    throw InvalidMaterialCoefficientsException("Material coefficients out of range.");
  }
}

}  // namespace util